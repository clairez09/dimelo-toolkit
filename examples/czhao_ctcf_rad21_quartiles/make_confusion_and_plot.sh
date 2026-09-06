#!/usr/bin/env bash
set -euo pipefail

# Sibling scripts and data live next to this file (symlinks resolved).
SELFDIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)

ml biology bedtools || true

usage () {
  cat <<'EOF'
Compute a 4x4 matrix with:
  rows = RAD21 quartiles (Q1..Q4)
  cols = CTCF quartiles (Q1..Q4)

Cell(i,j):
  unique mode = # unique CTCF Qj intervals that overlap >=1 RAD21 Qi interval
  pairs  mode = total # overlap pairs from CTCF intervals to RAD21 intervals

Expected BED naming:
  <RAD_PREFIX>.Q1.bed ... <RAD_PREFIX>.Q4.bed
  <CTCF_PREFIX>.Q1.bed ... <CTCF_PREFIX>.Q4.bed

Tolerance (bedtools v2.30.0 compatible):
  "within T bp counts as overlap" is implemented by expanding RAD21 intervals
  by +/-T using bedtools slop (-b T). Requires a genome sizes file (-g) if T>0.

Usage:
  ./make_confusion_and_plot.sh [-r RAD_PREFIX] [-c CTCF_PREFIX]
                               [-m unique|pairs] [-f MIN_FRAC] [-R]
                               [-t TOL_BP] [-g genome.sizes]

Options:
  -r   RAD prefix (default: RAD21.chm13.median)
  -c   CTCF prefix (default: CTCF.chm13.median)
  -m   mode [unique|pairs] (default: unique)
  -f   min fraction overlap of A (CTCF) required (bedtools -f). Example: 0.5
  -R   reciprocal overlap (bedtools -r); only useful with -f
  -t   tolerance in bp (slop RAD21 by +/-t). Default: 150. Use 0 for exact.
  -g   genome sizes file (default: chm13v2.chrom.sizes next to this script)
  -h   help
EOF
}

RAD_PREFIX="RAD21.chm13.median"
CTCF_PREFIX="CTCF.chm13.median"
MODE="unique"
MIN_FRAC=""
RECIP=""
TOL_BP="150"
GENOME="$SELFDIR/chm13v2.chrom.sizes"

while getopts ":r:c:m:f:Rt:g:h" opt; do
  case $opt in
    r) RAD_PREFIX="$OPTARG" ;;
    c) CTCF_PREFIX="$OPTARG" ;;
    m) MODE="$OPTARG" ;;
    f) MIN_FRAC="$OPTARG" ;;
    R) RECIP="yes" ;;
    t) TOL_BP="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: invalid option -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "ERROR: option -$OPTARG requires an argument" >&2; usage; exit 1 ;;
  esac
done

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found in PATH"; exit 1; }
command -v python3  >/dev/null 2>&1 || { echo "ERROR: python3 not found in PATH"; exit 1; }

# Optional intersect args
declare -a BT_ARGS
BT_ARGS=()
if [[ -n "${MIN_FRAC:-}" ]]; then BT_ARGS+=("-f" "$MIN_FRAC"); fi
if [[ -n "${RECIP:-}" ]]; then BT_ARGS+=("-r"); fi

rad=("${RAD_PREFIX}.Q1.bed" "${RAD_PREFIX}.Q2.bed" "${RAD_PREFIX}.Q3.bed" "${RAD_PREFIX}.Q4.bed")
ctc=("${CTCF_PREFIX}.Q1.bed" "${CTCF_PREFIX}.Q2.bed" "${CTCF_PREFIX}.Q3.bed" "${CTCF_PREFIX}.Q4.bed")

for f in "${rad[@]}" "${ctc[@]}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty file: $f" >&2; exit 1; }
done

if [[ "$TOL_BP" != "0" ]]; then
  [[ -s "$GENOME" ]] || { echo "ERROR: genome sizes missing/empty: $GENOME (required for -t > 0)" >&2; exit 1; }
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# Slop RAD21 if tolerance enabled
rad_use=()
if [[ "$TOL_BP" == "0" ]]; then
  rad_use=("${rad[@]}")
else
  for r in "${rad[@]}"; do
    out="$tmpdir/$(basename "$r" .bed).slop${TOL_BP}.bed"
    bedtools slop -b "$TOL_BP" -i "$r" -g "$GENOME" > "$out"
    rad_use+=("$out")
  done
fi

OUT_TSV="${RAD_PREFIX}__x__${CTCF_PREFIX}__CTCFcount__byRAD__${MODE}__tol${TOL_BP}bp.tsv"
OUT_PREFIX="${RAD_PREFIX}__x__${CTCF_PREFIX}__CTCFcount__byRAD__${MODE}__tol${TOL_BP}bp"

printf "RAD21\\\\CTCF\tQ1\tQ2\tQ3\tQ4\n" > "$OUT_TSV"

for i in 0 1 2 3; do
  r="${rad_use[$i]}"
  rq="Q$((i+1))"
  printf "%s" "$rq" >> "$OUT_TSV"

  for c in "${ctc[@]}"; do
    if [[ "$MODE" == "unique" ]]; then
      n=$(
        bedtools intersect ${BT_ARGS[@]+"${BT_ARGS[@]}"} \
          -u -a "$c" -b "$r" | wc -l | awk '{print $1}'
      )
    elif [[ "$MODE" == "pairs" ]]; then
      n=$(
        bedtools intersect ${BT_ARGS[@]+"${BT_ARGS[@]}"} \
          -c -a "$c" -b "$r" | awk '{s+=$NF} END{print s+0}'
      )
    else
      echo "ERROR: MODE must be 'unique' or 'pairs' (got: $MODE)" >&2
      exit 1
    fi
    printf "\t%s" "$n" >> "$OUT_TSV"
  done
  printf "\n" >> "$OUT_TSV"
done

echo "Wrote matrix: $OUT_TSV"
python3 "$SELFDIR/plot_heatmap_from_tsv.py" "$OUT_TSV" "$OUT_PREFIX"
echo "Wrote heatmaps: ${OUT_PREFIX}.heatmap.png  ${OUT_PREFIX}.heatmap.pdf"
