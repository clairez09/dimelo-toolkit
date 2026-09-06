#!/usr/bin/env python3
"""Keep only summits that have no neighbouring summit within --dist bp on the
same chromosome, then emit BED6+1: chrom start end name score strand summit.

Two input flavours:
  --format motif       BED6 of FIMO motif matches; summit = centre of the match
  --format narrowpeak  ENCODE narrowPeak;          summit = start + col10, score forced to 0

This is the score/strand-preserving generalisation of
drop_nearby_summits_incl_interval.py, which emitted only
5 columns.  Verified to reproduce both CTCF.filteredPeaks.withSummit.bed and
RAD21.filteredPeaks.withSummit.bed byte-for-byte.
"""
import argparse
import pandas as pd

MOTIF_COLS = ["chrom", "start", "end", "name", "score", "strand"]
NP_COLS = ["chrom", "start", "end", "name", "score", "strand",
           "signal", "pval", "qval", "peak"]
OUT_COLS = ["chrom", "start", "end", "name", "score", "strand", "summit"]


def natural_chrom_key(c):
    """chr1 < chr2 < ... < chr22 < chrM < chrX < chrY (karyotypic order)."""
    s = str(c)
    s = s[3:] if s.startswith("chr") else s
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--format", choices=["motif", "narrowpeak"], required=True)
    ap.add_argument("--dist", type=int, default=1000,
                    help="drop a summit if another summit is within this many bp")
    ap.add_argument("--prefix", default="CTCF_peak",
                    help="name prefix for surviving records (default: CTCF_peak)")
    ap.add_argument("--chrom-order", choices=["natural", "lexicographic"],
                    default="natural",
                    help="natural = chr1,chr2,...,chr22,chrX (karyotypic); "
                         "lexicographic = chr1,chr10,chr11,...,chr2,... "
                         "The two delivered files used different orders: CTCF "
                         "is natural, RAD21 is lexicographic.")
    a = ap.parse_args()

    if a.format == "motif":
        # score kept as a string so FIMO's "11" does not become "11.0"
        df = pd.read_csv(a.inp, sep="\t", header=None, names=MOTIF_COLS,
                         dtype={"score": str})
        df["summit"] = df["start"] + (df["end"] - df["start"]) // 2
    else:
        df = pd.read_csv(a.inp, sep="\t", header=None, names=NP_COLS, comment="#")
        df["summit"] = df["start"] + df["peak"].astype(int)
        df["score"] = "0"

    if a.chrom_order == "natural":
        order = sorted(df["chrom"].unique(), key=natural_chrom_key)
        df["chrom"] = pd.Categorical(df["chrom"], categories=order, ordered=True)

    df = df.sort_values(["chrom", "summit"]).reset_index(drop=True)

    prev_summit = df.groupby("chrom", observed=True)["summit"].shift(1)
    next_summit = df.groupby("chrom", observed=True)["summit"].shift(-1)
    keep = ~(((df["summit"] - prev_summit) <= a.dist).fillna(False) |
             ((next_summit - df["summit"]) <= a.dist).fillna(False))

    out = df.loc[keep].copy()
    out["name"] = ["%s_%d" % (a.prefix, i + 1) for i in range(len(out))]
    out[OUT_COLS].to_csv(a.out, sep="\t", header=False, index=False)

    print("[filter_isolated_summits] %s: %d -> %d kept (dist=%d)"
          % (a.inp, len(df), len(out), a.dist))


if __name__ == "__main__":
    main()
