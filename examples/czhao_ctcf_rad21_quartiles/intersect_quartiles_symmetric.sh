#!/usr/bin/env bash
set -euo pipefail

OUT="${OUT:-./intersections_150bp}"
W="${W:-150}"
mkdir -p "$OUT"

ml biology bedtools

CTCF_H="${CTCF_H:-./CTCF.chm13.median.Q4.bed}"
CTCF_L="${CTCF_L:-./CTCF.chm13.median.Q1.bed}"
RAD21_H="${RAD21_H:-./RAD21.chm13.median.Q4.bed}"
RAD21_L="${RAD21_L:-./RAD21.chm13.median.Q1.bed}"

# Symmetric "within W bp" using bedtools window
bedtools window -u -w "$W" -a "$CTCF_H" -b "$RAD21_H" > "$OUT/CTCFhigh_RAD21high.w${W}.bed"
bedtools window -u -w "$W" -a "$CTCF_H" -b "$RAD21_L" > "$OUT/CTCFhigh_RAD21low.w${W}.bed"
bedtools window -u -w "$W" -a "$CTCF_L" -b "$RAD21_H" > "$OUT/CTCFlow_RAD21high.w${W}.bed"
bedtools window -u -w "$W" -a "$CTCF_L" -b "$RAD21_L" > "$OUT/CTCFlow_RAD21low.w${W}.bed"

echo "Wrote outputs to: $OUT"
