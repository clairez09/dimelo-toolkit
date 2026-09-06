#!/usr/bin/env python3
"""
Given N BED region sets (e.g. CTCF/RAD21 high/low intersections),
generate pileups + enrichment plots for CTCF targetting and control non specific.

Adds:
  --skip_existing_pileups : if pileup.sorted.bed.gz already exists, do not rerun parse_bam.pileup()

Defaults:
  --motifs "A,0"
  --window 2000
  --cores 16
  --ymin 0
  --ymax 0.35
  --dpi 200
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt

from dimelo import parse_bam
from dimelo.plot_enrichment_profile import plot_enrichment_profile as plot_profile


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


def plot_profile_fixed_ylim(mod_file, regions_bed, motifs, sample_name, window_size, cores,
                            regions_5to3prime, ymin, ymax, out_png, dpi):
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
    ax.set_ylim(ymin, ymax)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def parse_named_beds(items):
    """Parse: --beds name=path name=path ..."""
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
        "--beds",
        required=True,
        nargs="+",
        help="Region sets as name=path entries (e.g. HH=/path/CTCFhigh_RAD21high.bed)",
    )

    ap.add_argument("--outdir", required=True, type=Path)

    ap.add_argument("--motifs", default="A,0")
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--cores", type=int, default=16)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=0.35)
    ap.add_argument("--dpi", type=int, default=200)

    ap.add_argument("--quiet", action="store_true", default=False)
    ap.add_argument("--no_cleanup", action="store_true")
    ap.add_argument(
        "--skip_existing_pileups",
        action="store_true",
        default=False,
        help="If pileup.sorted.bed.gz exists, do not rerun parse_bam.pileup()",
    )

    ap.add_argument("--regions_5to3prime", dest="regions_5to3prime", action="store_true", default=True)
    ap.add_argument("--no_regions_5to3prime", dest="regions_5to3prime", action="store_false")

    args = ap.parse_args()

    motifs = [args.motifs]
    cleanup = (not args.no_cleanup)
    beds = parse_named_beds(args.beds)

    samples = [
        ("CTCF_targetting", args.CTCF_targetting),
        ("control_non_specific", args.control_non_specific),
    ]

    pileup_root = args.outdir / "pileups"
    plot_root = args.outdir / "plots"
    pileup_root.mkdir(parents=True, exist_ok=True)
    plot_root.mkdir(parents=True, exist_ok=True)

    for set_name, bed in beds:
        if not bed.exists():
            raise FileNotFoundError(f"Missing bed for set {set_name}: {bed}")

        for condition_name, bam in samples:
            set_pileup_dir = pileup_root / set_name
            prefix = f"{condition_name}_{set_name}_pileup"
            expected_mod_file = set_pileup_dir / prefix / "pileup.sorted.bed.gz"

            if args.skip_existing_pileups and expected_mod_file.exists():
                mod_file = expected_mod_file
                if not args.quiet:
                    print(f"[skip] Using existing pileup: {mod_file}")
            else:
                mod_file = run_pileup(
                    bam=bam,
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

            out_png = (
                plot_root
                / set_name
                / condition_name
                / f"{condition_name}_{set_name}_enrichment_y{args.ymin}-{args.ymax}.png"
            )
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
            )

    print(f"Done.\nPileups: {pileup_root}\nPlots:  {plot_root}")


if __name__ == "__main__":
    main()
