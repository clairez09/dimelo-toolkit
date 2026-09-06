#!/usr/bin/env python3
"""
Optimized equalization + pileup + enrichment plots for multiple BED region sets.

Changes included:
  (Option A) Filter out chrM reads from BAM streams used for counting/subsampling.
  (Robust fix) Add bedtools '-g <bam.genome>' everywhere we intersect against a BAM
               stream (stdin) so bedtools does not error if the BAM contains contigs
               (e.g., chrY) not present in the BED.

NEW:
  - Add 30 bp smoothing ONLY (plot smoothing). Does NOT change pileup generation
    or read counting/equalization. Smoothing is applied to plotted profiles via
    a rolling mean over y-values.

Requires:
  samtools, bedtools, awk, python3, dimelo, matplotlib
"""

import argparse
import subprocess
import shlex
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np  # NEW

from dimelo import parse_bam
from dimelo.plot_enrichment_profile import plot_enrichment_profile as plot_profile

FLAG_FILTER = "0x904"
EXCLUDE_RNAME = "chrM"   # Option A: drop mitochondrial contig from BAM stream


def run(cmd, quiet=False, capture_stdout=False):
    if not quiet:
        print(f"[cmd] {cmd}")
    if capture_stdout:
        return subprocess.check_output(cmd, shell=True, text=True)
    subprocess.check_call(cmd, shell=True)
    return None


def ensure_bam_index(bam: Path, quiet=False, threads=4):
    bai1 = bam.with_suffix(bam.suffix + ".bai")  # foo.bam.bai
    bai2 = bam.with_suffix(".bai")               # foo.bai
    if bai1.exists() or bai2.exists():
        return
    run(f"samtools index -@ {threads} {shlex.quote(str(bam))}", quiet=quiet)


def bam_genome_file(bam: Path, cache_dir: Path, quiet: bool) -> Path:
    """
    Create a bedtools 'genome' file (contig<TAB>length) IN THE SAME ORDER as the BAM header (@SQ).
    This lets us sort BEDs to match the BAM contig order so bedtools intersect -sorted is fast and safe.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{bam.stem}.genome"
    if out.exists():
        return out

    cmd = (
        f"samtools view -H {shlex.quote(str(bam))} "
        f"| awk -F'\\t' '$1==\"@SQ\"{{"
        f"sn=\"\"; ln=\"\";"
        f"for(i=1;i<=NF;i++){{"
        f" if($i ~ /^SN:/){{sn=$i; sub(/^SN:/,\"\",sn)}}"
        f" if($i ~ /^LN:/){{ln=$i; sub(/^LN:/,\"\",ln)}}"
        f"}}"
        f"if(sn!=\"\" && ln!=\"\") print sn\"\\t\"ln"
        f"}}' > {shlex.quote(str(out))}"
    )
    run(cmd, quiet=quiet)
    return out


def sort_bed_once_bamorder(bed: Path, bam_genome: Path, cache_dir: Path, quiet: bool) -> Path:
    """
    Sort BED in BAM contig order using bedtools sort -g <bam.genome>.
    Cached by mtime.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / f"{bed.stem}.bamorder.sorted.bed"
    if out.exists() and out.stat().st_mtime >= bed.stat().st_mtime:
        return out

    cmd = (
        f"bedtools sort -g {shlex.quote(str(bam_genome))} "
        f"-i {shlex.quote(str(bed))} > {shlex.quote(str(out))}"
    )
    run(cmd, quiet=quiet)
    return out


def count_overlaps_sum_per_interval(
    bam: Path,
    bed_bamorder_sorted: Path,
    bam_genome: Path,
    quiet: bool,
    threads: int,
) -> int:
    """
    For each BED interval, count #overlapping reads; then sum counts across intervals.
    A read overlapping multiple intervals contributes multiple times.

    Fixes:
      - Excludes chrM reads (Option A)
      - Uses bedtools intersect -sorted -g <bam.genome> to avoid contig mismatch errors
        when BAM contains contigs not present in BED (e.g., chrY).
    """
    cmd = (
        f"samtools view -@ {threads} -h -F {FLAG_FILTER} {shlex.quote(str(bam))} "
        f"| awk 'BEGIN{{OFS=\"\\t\"}} "
        f"     /^@/{{print; next}} "
        f"     $3!=\"{EXCLUDE_RNAME}\" {{print}}' "
        f"| samtools view -@ {threads} -b - "
        f"| bedtools intersect -sorted -g {shlex.quote(str(bam_genome))} "
        f"    -a {shlex.quote(str(bed_bamorder_sorted))} "
        f"    -b - "
        f"    -c "
        f"| awk '{{s+=$NF}} END{{print s+0}}'"
    )
    out = run(cmd, quiet=quiet, capture_stdout=True).strip()
    return int(out) if out else 0


def make_subsampled_bam_from_regions(
    bam: Path,
    bed_for_L: Path,
    out_bam: Path,
    fraction: float,
    seed: int,
    threads: int,
    quiet: bool,
):
    """
    Build ONE coordinate-sorted, indexed BAM containing only reads overlapping bed_for_L,
    filtered by FLAG_FILTER, excluding chrM (Option A), and downsampled to fraction.
    """
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if fraction <= 0:
        raise ValueError(f"Downsample fraction <= 0 ({fraction}); cannot proceed.")

    base_stream = (
        f"samtools view -@ {threads} -h -F {FLAG_FILTER} "
        f"-L {shlex.quote(str(bed_for_L))} {shlex.quote(str(bam))} "
        f"| awk 'BEGIN{{OFS=\"\\t\"}} "
        f"     /^@/{{print; next}} "
        f"     $3!=\"{EXCLUDE_RNAME}\" {{print}}' "
        f"| samtools view -@ {threads} -b - "
    )

    if fraction >= 0.999999:
        cmd = base_stream + f"| samtools sort -@ {threads} -o {shlex.quote(str(out_bam))} -"
    else:
        frac_str = f"{seed}.{int(fraction * 10_000_000):07d}"
        cmd = (
            base_stream
            + f"| samtools view -@ {threads} -b -s {frac_str} - "
            + f"| samtools sort -@ {threads} -o {shlex.quote(str(out_bam))} -"
        )

    run(cmd, quiet=quiet)
    ensure_bam_index(out_bam, quiet=quiet, threads=threads)


def run_pileup(bam, fasta, regions_bed, outdir, pileup_prefix, motifs, window_size, cores, quiet, cleanup):
    outdir.mkdir(parents=True, exist_ok=True)
    parse_bam.pileup(
        str(bam),
        pileup_prefix,
        str(fasta),
        str(outdir),
        str(regions_bed),
        motifs=motifs,
        window_size=window_size,
        cores=cores,
        quiet=quiet,
        cleanup=cleanup,
    )
    mod_file = outdir / pileup_prefix / "pileup.sorted.bed.gz"
    if not mod_file.exists():
        raise FileNotFoundError(f"Expected pileup output not found: {mod_file}")
    return mod_file


# NEW: smoothing helper
def _smooth_line_xy(x, y, smooth_bp: int):
    """
    Rolling-mean smoothing by ~smooth_bp along the x-axis.
    Converts bp window to number of points using median dx.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if smooth_bp is None or smooth_bp <= 0:
        return x, y
    if len(x) < 3:
        return x, y

    dx = np.diff(x)
    dx_med = np.median(dx[np.isfinite(dx)]) if np.any(np.isfinite(dx)) else 1.0
    if dx_med <= 0 or not np.isfinite(dx_med):
        dx_med = 1.0

    win_pts = int(round(smooth_bp / dx_med))
    if win_pts < 2:
        return x, y
    if win_pts % 2 == 0:
        win_pts += 1  # odd window is nicer visually

    # Edge-padded convolution (prevents shrinkage)
    pad = win_pts // 2
    ypad = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(win_pts, dtype=float) / win_pts
    ys = np.convolve(ypad, kernel, mode="valid")
    return x, ys


def plot_profile_fixed_ylim(
    mod_file, regions_bed, motifs, sample_name, window_size, cores,
    regions_5to3prime, ymin, ymax, out_png, dpi,
    corner_text=None, corner_loc="upper left",
    smooth_bp: int = 30,   # NEW (default 30 bp smoothing)
):
    plot_profile(
        mod_file_names=[mod_file],
        regions_list=[regions_bed],
        motifs=motifs,
        sample_names=[sample_name],
        window_size=window_size,
        cores=cores,
        regions_5to3prime=regions_5to3prime,
    )
    ax = plt.gca()

    # NEW: apply smoothing ONLY to plotted lines
    if smooth_bp and smooth_bp > 0:
        for line in ax.lines:
            x = line.get_xdata()
            y = line.get_ydata()
            xs, ys = _smooth_line_xy(x, y, smooth_bp=smooth_bp)
            line.set_xdata(xs)
            line.set_ydata(ys)

        # Update limits after editing data (keep your fixed y-lim below)
        ax.relim()
        ax.autoscale_view(scalex=True, scaley=False)

    ax.set_ylim(ymin, ymax)

    if corner_text:
        locs = {
            "upper left":  (0.02, 0.98, "left",  "top"),
            "upper right": (0.98, 0.98, "right", "top"),
            "lower left":  (0.02, 0.02, "left",  "bottom"),
            "lower right": (0.98, 0.02, "right", "bottom"),
        }
        x, y, ha, va = locs[corner_loc]
        ax.text(
            x, y, corner_text,
            transform=ax.transAxes,
            ha=ha, va=va,
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="black", alpha=0.85)
        )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def parse_named_beds(items):
    out = []
    for it in items:
        if "=" not in it:
            raise ValueError(f"Bad --beds entry (expected name=path): {it}")
        name, path = it.split("=", 1)
        name = name.strip()
        path = Path(path)
        if not name:
            raise ValueError(f"Empty name in --beds entry: {it}")
        out.append((name, path))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--CTCF_targetting", required=True, type=Path)
    ap.add_argument("--control_non_specific", required=True, type=Path)
    ap.add_argument("--fasta", required=True, type=Path)

    ap.add_argument(
        "--beds", required=True, nargs="+",
        help="Region sets as name=path entries (e.g. HH=/path/HH.bed HL=/path/HL.bed ...)"
    )

    ap.add_argument("--outdir", required=True, type=Path)

    ap.add_argument("--motifs", default="A,0")
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--cores", type=int, default=16)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=0.35)
    ap.add_argument("--dpi", type=int, default=200)

    ap.add_argument("--quiet", action="store_true", default=False)
    ap.add_argument("--no_cleanup", action="store_true")

    ap.add_argument(
        "--skip_existing_pileups", action="store_true", default=False,
        help="If pileup.sorted.bed.gz exists, do not rerun parse_bam.pileup()"
    )

    ap.add_argument("--regions_5to3prime", dest="regions_5to3prime", action="store_true", default=True)
    ap.add_argument("--no_regions_5to3prime", dest="regions_5to3prime", action="store_false")

    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument(
        "--corner_loc", default="upper left",
        choices=["upper left", "upper right", "lower left", "lower right"]
    )
    ap.add_argument("--report_original_counts", action="store_true", default=False)

    # smoothing control
    ap.add_argument(
        "--smooth_bp", type=int, default=30,
        help="Smoothing window in bp applied ONLY to plotted enrichment profiles (0 disables). Default: 30"
    )

    args = ap.parse_args()

    motifs = [args.motifs]
    cleanup = (not args.no_cleanup)
    beds = parse_named_beds(args.beds)

    samples = [
        ("CTCF_targetting", args.CTCF_targetting),
        ("control_non_specific", args.control_non_specific),
    ]

    if not args.fasta.exists():
        raise FileNotFoundError(f"Missing fasta: {args.fasta}")
    for set_name, bed in beds:
        if not bed.exists():
            raise FileNotFoundError(f"Missing bed for set {set_name}: {bed}")
    for _, bam in samples:
        if not bam.exists():
            raise FileNotFoundError(f"Missing bam: {bam}")
        ensure_bam_index(bam, quiet=args.quiet, threads=args.threads)

    outdir = args.outdir
    pileup_root = outdir / "pileups"
    plot_root = outdir / "plots"
    tmp_root = outdir / "subsampled_bams"

    genome_cache = outdir / "bam_genomes"
    bed_cache = outdir / "beds_bamorder_cache"

    pileup_root.mkdir(parents=True, exist_ok=True)
    plot_root.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)
    genome_cache.mkdir(parents=True, exist_ok=True)
    bed_cache.mkdir(parents=True, exist_ok=True)

    for condition_name, bam in samples:
        if not args.quiet:
            print(f"\n=== {condition_name} ===")

        bam_genome = bam_genome_file(bam, cache_dir=genome_cache, quiet=args.quiet)

        beds_prepped = []
        per_bam_bed_cache = bed_cache / condition_name
        per_bam_bed_cache.mkdir(parents=True, exist_ok=True)

        for set_name, bed in beds:
            bed_bamorder_sorted = sort_bed_once_bamorder(
                bed=bed,
                bam_genome=bam_genome,
                cache_dir=per_bam_bed_cache,
                quiet=args.quiet,
            )
            beds_prepped.append((set_name, bed, bed_bamorder_sorted))

        counts = {}
        for set_name, bed, bed_bamorder_sorted in beds_prepped:
            c = count_overlaps_sum_per_interval(
                bam=bam,
                bed_bamorder_sorted=bed_bamorder_sorted,
                bam_genome=bam_genome,
                quiet=args.quiet,
                threads=args.threads,
            )
            counts[set_name] = c
            if not args.quiet:
                print(f"[count] {condition_name} {set_name}: {c}")

        target = min(counts.values()) if counts else 0
        if target <= 0:
            raise RuntimeError(f"{condition_name}: target count is {target}. Are beds empty / no overlaps?")

        if not args.quiet:
            print(f"[target] {condition_name}: equalizing all sets to {target}")

        for set_name, bed, bed_bamorder_sorted in beds_prepped:
            set_pileup_dir = pileup_root / set_name
            prefix = f"{condition_name}_{set_name}_pileup"
            expected_mod_file = set_pileup_dir / prefix / "pileup.sorted.bed.gz"

            subsampled_bam = tmp_root / condition_name / set_name / f"{condition_name}.{set_name}.subsampled.bam"

            orig = counts[set_name]
            frac = (target / orig) if orig > 0 else 0.0

            if not subsampled_bam.exists():
                make_subsampled_bam_from_regions(
                    bam=bam,
                    bed_for_L=bed_bamorder_sorted,
                    out_bam=subsampled_bam,
                    fraction=frac,
                    seed=args.seed,
                    threads=args.threads,
                    quiet=args.quiet,
                )

            if args.skip_existing_pileups and expected_mod_file.exists():
                mod_file = expected_mod_file
                if not args.quiet:
                    print(f"[skip] Using existing pileup: {mod_file}")
            else:
                mod_file = run_pileup(
                    bam=subsampled_bam,
                    fasta=args.fasta,
                    regions_bed=bed,
                    outdir=set_pileup_dir,
                    pileup_prefix=prefix,
                    motifs=motifs,
                    window_size=args.window,
                    cores=args.cores,
                    quiet=args.quiet,
                    cleanup=cleanup,
                )

            out_png = plot_root / set_name / condition_name / (
                f"{condition_name}_{set_name}_enrichment_equalized_y{args.ymin}-{args.ymax}.png"
            )
            if args.report_original_counts:
                corner_text = f"equalized={target}\norig={orig}\nfrac={frac:.4f}"
            else:
                corner_text = f"equalized={target}"

            plot_profile_fixed_ylim(
                mod_file=mod_file,
                regions_bed=bed,
                motifs=motifs,
                sample_name=f"{condition_name}_{set_name}",
                window_size=args.window,
                cores=args.cores,
                regions_5to3prime=args.regions_5to3prime,
                ymin=args.ymin,
                ymax=args.ymax,
                out_png=out_png,
                dpi=args.dpi,
                corner_text=corner_text,
                corner_loc=args.corner_loc,
                smooth_bp=args.smooth_bp,
            )

    print(f"\nDone.\nSubsampled BAMs: {tmp_root}\nPileups: {pileup_root}\nPlots:  {plot_root}")


if __name__ == "__main__":
    main()
