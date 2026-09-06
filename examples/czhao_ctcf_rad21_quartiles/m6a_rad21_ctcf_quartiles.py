#!/usr/bin/env python3
"""
Generate all 16 m6A enrichment profile plots:
- 2 conditions (CTCF targetting, control non specific)
- 2 factor sets (RAD21 and CTCF), each with 4 quartile BEDs
Total = 2 * 2 * 4 = 16 plots (and 16 pileups).

Defaults (override via CLI if desired):
  --motifs "A,0"
  --window 2000
  --cores 16
  --ymin 0
  --ymax 0.45
  --dpi 200
  --regions_5to3prime True

Outputs:
  {outdir}/pileups/{factor}/{condition}_{factor}_Q{1-4}_pileup/pileup.sorted.bed.gz
  {outdir}/plots/{factor}/{condition}/{condition}_{factor}_Q{1-4}_enrichment_y{ymin}-{ymax}.png
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from dimelo import parse_bam
from dimelo.plot_enrichment_profile import plot_enrichment_profile as plot_profile


def run_pileup(
    bam: Path,
    fasta: Path,
    regions_bed: Path,
    outdir: Path,
    pileup_prefix: str,
    motifs,
    window_size: int,
    cores: int,
    quiet: bool,
    cleanup: bool,
) -> Path:
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


def plot_profile_fixed_ylim(
    mod_file: Path,
    regions_bed: Path,
    motifs,
    sample_name: str,
    window_size: int,
    cores: int,
    regions_5to3prime: bool,
    ymin: float,
    ymax: float,
    out_png: Path,
    dpi: int,
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

    # Force identical y-scale for all plots
    ax = plt.gca()
    ax.set_ylim(ymin, ymax)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description="Generate m6A pileups + enrichment profiles for RAD21/CTCF quartiles (CTCF targetting vs control non specific)."
    )

    # Inputs
    ap.add_argument(
        "--CTCF_targetting",
        required=True,
        type=Path,
        help="CTCF targetting BAM",
    )
    ap.add_argument(
        "--control_non_specific",
        required=True,
        type=Path,
        help="control non specific BAM",
    )
    ap.add_argument("--fasta", required=True, type=Path, help="Reference FASTA")

    # Quartile BEDs for each factor
    ap.add_argument("--rad21_q1", required=True, type=Path)
    ap.add_argument("--rad21_q2", required=True, type=Path)
    ap.add_argument("--rad21_q3", required=True, type=Path)
    ap.add_argument("--rad21_q4", required=True, type=Path)

    ap.add_argument("--ctcf_q1", required=True, type=Path)
    ap.add_argument("--ctcf_q2", required=True, type=Path)
    ap.add_argument("--ctcf_q3", required=True, type=Path)
    ap.add_argument("--ctcf_q4", required=True, type=Path)

    # dimelo / plotting params (defaults set here)
    ap.add_argument("--motifs", default="A,0", help='Motif string like "A,0"')
    ap.add_argument("--window", type=int, default=2000, help="Window size around region center")
    ap.add_argument("--cores", type=int, default=16)
    ap.add_argument("--ymin", type=float, default=0.0)
    ap.add_argument("--ymax", type=float, default=0.45)
    ap.add_argument("--dpi", type=int, default=200, help="PNG DPI for saved figures")

    # regions_5to3prime: default True, allow disabling with --no_regions_5to3prime
    ap.add_argument(
        "--regions_5to3prime",
        dest="regions_5to3prime",
        action="store_true",
        default=True,
        help="Orient regions 5'->3' before plotting (default: True)",
    )
    ap.add_argument(
        "--no_regions_5to3prime",
        dest="regions_5to3prime",
        action="store_false",
        help="Do not orient regions 5'->3'",
    )

    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--quiet", action="store_true", default=False)
    ap.add_argument("--no_cleanup", action="store_true", help="Do not delete intermediates")

    args = ap.parse_args()

    motifs = [args.motifs]
    cleanup = (not args.no_cleanup)

    samples = {
        "CTCF targetting": args.CTCF_targetting,
        "control non specific": args.control_non_specific,
    }

    factors = {
        "RAD21": {
            "Q1": args.rad21_q1,
            "Q2": args.rad21_q2,
            "Q3": args.rad21_q3,
            "Q4": args.rad21_q4,
        },
        "CTCF": {
            "Q1": args.ctcf_q1,
            "Q2": args.ctcf_q2,
            "Q3": args.ctcf_q3,
            "Q4": args.ctcf_q4,
        },
    }

    pileup_root = args.outdir / "pileups"
    plot_root = args.outdir / "plots"
    pileup_root.mkdir(parents=True, exist_ok=True)
    plot_root.mkdir(parents=True, exist_ok=True)

    # 2 conditions × 2 factors × 4 quartiles = 16
    for factor_name, quartiles in factors.items():
        factor_pileup_dir = pileup_root / factor_name

        for condition_name, bam in samples.items():
            for qname, bed in quartiles.items():
                pileup_prefix = f"{condition_name}_{factor_name}_{qname}_pileup"

                mod_file = run_pileup(
                    bam=bam,
                    fasta=args.fasta,
                    regions_bed=bed,
                    outdir=factor_pileup_dir,
                    pileup_prefix=pileup_prefix,
                    motifs=motifs,
                    window_size=args.window,
                    cores=args.cores,
                    quiet=args.quiet,
                    cleanup=cleanup,
                )

                out_png = (
                    plot_root
                    / factor_name
                    / condition_name
                    / f"{condition_name}_{factor_name}_{qname}_enrichment_y{args.ymin}-{args.ymax}.png"
                )

                plot_profile_fixed_ylim(
                    mod_file=mod_file,
                    regions_bed=bed,
                    motifs=motifs,
                    sample_name=f"{condition_name}_{factor_name}_{qname}",
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
