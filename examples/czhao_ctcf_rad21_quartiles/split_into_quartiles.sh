#!/usr/bin/env bash
set -euo pipefail

inbed="$1"                 # e.g. FACTOR.lifted.chm13v2.bed
outprefix="${2:-$(basename "$inbed" .bed)}"  # prefix for outputs

# Sort by score (col 5). If score is different column, change here.
tmp_sorted=$(mktemp)
sort -k5,5n "$inbed" > "$tmp_sorted"

# Count rows
n=$(wc -l < "$tmp_sorted")

# Compute cutoffs by row index (rank-based quartiles)
# Q1: [1 .. floor(n/4)]
# Q2: [floor(n/4)+1 .. floor(n/2)]
# Q3: [floor(n/2)+1 .. floor(3n/4)]
# Q4: [floor(3n/4)+1 .. n]
q1=$(( n/4 ))
q2=$(( n/2 ))
q3=$(( (3*n)/4 ))

awk -v q1="$q1" -v q2="$q2" -v q3="$q3" -v p="$outprefix" 'BEGIN{
  OFS="\t"
}
{
  if (NR <= q1)       print > (p ".Q1.bed");
  else if (NR <= q2)  print > (p ".Q2.bed");
  else if (NR <= q3)  print > (p ".Q3.bed");
  else                print > (p ".Q4.bed");
}' "$tmp_sorted"

rm -f "$tmp_sorted"

echo "Wrote:"
ls -1 "${outprefix}.Q"{1..4}".bed"
