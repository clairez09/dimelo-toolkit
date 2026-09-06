#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print("Usage: plot_heatmap_from_tsv.py <matrix.tsv> <out_prefix>", file=sys.stderr)
        sys.exit(1)

    tsv, out_prefix = sys.argv[1], sys.argv[2]
    df = pd.read_csv(tsv, sep="\t", index_col=0)

    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(6, 5))
    ax = sns.heatmap(df, annot=True, fmt="d", cmap="viridis",
                     cbar_kws={"label": "# overlaps"})
    ax.set_xlabel("CTCF quartile")
    ax.set_ylabel("RAD21 quartile")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.heatmap.png", dpi=300)
    plt.savefig(f"{out_prefix}.heatmap.pdf")

if __name__ == "__main__":
    main()
