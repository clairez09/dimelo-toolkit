#!/usr/bin/env python3
# liftover_both.py: liftOver CTCF.median.bed and RAD21.median.bed hg38->CHM13

import argparse
import subprocess
from pathlib import Path

DEFAULT_FACTORS = ["CTCF", "RAD21"]

def run(cmd):
    subprocess.check_call(cmd)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--chain",
        default="/oak/stanford/groups/altemose/references/chains/hg38-chm13v2.over.chain",
        help="hg38->chm13 chain file",
    )
    ap.add_argument(
        "--indir",
        default=".",
        help="Directory containing <FACTOR>.median.bed (default: current dir)",
    )
    ap.add_argument(
        "--outdir",
        default=".",
        help="Where to write lifted beds (default: current dir)",
    )
    ap.add_argument(
        "--factors",
        nargs="*",
        default=DEFAULT_FACTORS,
        help="Factors to process (default: CTCF RAD21)",
    )
    ap.add_argument(
        "--liftover-bin",
        default="liftOver",
        help="liftOver executable name/path (default: liftOver)",
    )
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for factor in args.factors:
        in_bed = indir / f"{factor}.median.bed"
        out_bed = outdir / f"{factor}.chm13.median.bed"
        unmapped = outdir / f"{factor}.chm13.unmapped.bed"

        if not in_bed.exists():
            print(f"[SKIP] missing: {in_bed}")
            continue

        print(f"[RUN] {factor}: {in_bed} -> {out_bed}")
        run([args.liftover_bin, str(in_bed), str(args.chain), str(out_bed), str(unmapped)])

    print("[DONE]")

if __name__ == "__main__":
    main()
