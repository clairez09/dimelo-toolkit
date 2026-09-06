#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Robust logging + "where it failed"
# -----------------------------
LOGDIR="${SLURM_SUBMIT_DIR:-$PWD}/logs"
mkdir -p "$LOGDIR"
LOG="${LOGDIR}/${SLURM_JOB_NAME:-ctcf_rad21_confusion}.${SLURM_JOB_ID:-$$}.log"

exec > >(tee -a "$LOG") 2>&1
export PS4='+ [$(date "+%F %T")] ${BASH_SOURCE}:${LINENO}: '
set -x
trap 'rc=$?; echo; echo "ERROR (exit=$rc) at ${BASH_SOURCE}:${LINENO} while running:"; echo " $BASH_COMMAND"; echo "Log: $LOG"; exit $rc' ERR

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "PWD: $PWD"
echo "Log: $LOG"

# -----------------------------
# Modules
# -----------------------------
module purge
module load viz gcc/12.4.0
module load biology samtools bedtools
module load python/3.14.2 py-numpy/2.3.5_py314 py-pandas/2.3.3_py314 py-matplotlib/3.10.8_py314
# seaborn via: python3 -m pip install --user seaborn

# -----------------------------
# Defaults
# -----------------------------
RAD_PREFIX="RAD21.chm13.median"
CTCF_PREFIX="CTCF.chm13.median"

CTCF_TARGETTING="/oak/stanford/groups/altemose/data/20250110_NG_one_pot/barcode17.merged.sorted.bam"
CONTROL_NON_SPECIFIC="/oak/stanford/groups/altemose/data/20250110_NG_one_pot/barcode18.merged.sorted.bam"

OUTDIR="ctcf_vs_rad21_confusion_w150"
THREADS=4

# Proximity window for bedtools window (like your intersections_150bp script)
WINDOW="${WINDOW:-150}"

# Filtering matches your equalization script
FLAG_FILTER="0x904"
EXCLUDE_RNAME="chrM"

# If user supplies --beds, we use those instead of defaults
BEDS=()

usage() {
cat <<EOF
Usage: $0 [options]

Build a rectangular "confusion" matrix:

  rows = CTCF_* BED sets
  cols = RAD21_* BED sets

For each cell (CTCF_i, RAD21_j):
  regions = CTCF_i intervals that are within WINDOW bp of any RAD21_j interval
            (bedtools window -u -w WINDOW -a CTCF_i -b RAD21_j)
  count   = sum over those regions of (# overlapping reads), using cached filtered BAM:
            bedtools intersect -c per interval then sum

BAM filtering (cached once per BAM):
  samtools view -h -F 0x904
  remove alignments with RNAME == chrM
  sort + index

Options:
  --rad-prefix STR        (default: ${RAD_PREFIX})
  --ctcf-prefix STR       (default: ${CTCF_PREFIX})
  --beds NAME=PATH ...    (override all default beds)
  --CTCF_targetting PATH
  --control_non_specific PATH
  --outdir DIR
  --threads INT
  --window INT            (default: 150; can also set env WINDOW)
  -h, --help

Defaults (if no --beds):
  CTCF_Q1..Q4 = ${CTCF_PREFIX}.Q1-4.bed
  RAD21_Q1..Q4 = ${RAD_PREFIX}.Q1-4.bed

Outputs:
  counts.CTCF_targetting.tsv, matrix.CTCF_targetting.tsv, heatmap.CTCF_targetting.png
  counts.control_non_specific.tsv, matrix.control_non_specific.tsv, heatmap.control_non_specific.png
EOF
}

# -----------------------------
# Parse args
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rad-prefix) RAD_PREFIX="$2"; shift 2 ;;
    --ctcf-prefix) CTCF_PREFIX="$2"; shift 2 ;;
    --beds)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do
        BEDS+=("$1")
        shift
      done
      ;;
    --CTCF_targetting) CTCF_TARGETTING="$2"; shift 2 ;;
    --control_non_specific) CONTROL_NON_SPECIFIC="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --window) WINDOW="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$OUTDIR"
export OUTDIR

# -----------------------------
# Default BED sets
# -----------------------------
DEFAULT_BEDS=(
  "CTCF_Q1=${CTCF_PREFIX}.Q1.bed"
  "CTCF_Q2=${CTCF_PREFIX}.Q2.bed"
  "CTCF_Q3=${CTCF_PREFIX}.Q3.bed"
  "CTCF_Q4=${CTCF_PREFIX}.Q4.bed"
  "RAD21_Q1=${RAD_PREFIX}.Q1.bed"
  "RAD21_Q2=${RAD_PREFIX}.Q2.bed"
  "RAD21_Q3=${RAD_PREFIX}.Q3.bed"
  "RAD21_Q4=${RAD_PREFIX}.Q4.bed"
)

if [[ ${#BEDS[@]} -eq 0 ]]; then
  BEDS=("${DEFAULT_BEDS[@]}")
fi

# -----------------------------
# Validate inputs
# -----------------------------
for item in "${BEDS[@]}"; do
  [[ "$item" == *"="* ]] || { echo "Bad bed entry (expected NAME=PATH): $item" >&2; exit 1; }
  name="${item%%=*}"
  path="${item#*=}"
  [[ -n "$name" ]] || { echo "Empty bed name in: $item" >&2; exit 1; }
  [[ -s "$path" ]] || { echo "Missing/empty bed for $name: $path" >&2; exit 1; }
done

[[ -s "$CTCF_TARGETTING" ]] || { echo "Missing CTCF_TARGETTING BAM: $CTCF_TARGETTING" >&2; exit 1; }
[[ -s "$CONTROL_NON_SPECIFIC" ]] || { echo "Missing CONTROL_NON_SPECIFIC BAM: $CONTROL_NON_SPECIFIC" >&2; exit 1; }

ensure_index() {
  local bam="$1"
  if [[ -s "${bam}.bai" ]]; then
    :
  elif [[ -s "${bam%.bam}.bai" ]]; then
    :
  else
    samtools index -@ "$THREADS" "$bam"
  fi
}
ensure_index "$CTCF_TARGETTING"
ensure_index "$CONTROL_NON_SPECIFIC"

# -----------------------------
# Genome file in BAM header order for bedtools intersect -sorted
# -----------------------------
bam_genome_for() {
  local bam="$1"
  local out="$2"
  samtools view -H "$bam" \
  | awk -F'\t' '$1=="@SQ"{
      sn=""; ln="";
      for(i=1;i<=NF;i++){
        if($i ~ /^SN:/){sn=$i; sub(/^SN:/,"",sn)}
        if($i ~ /^LN:/){ln=$i; sub(/^LN:/,"",ln)}
      }
      if(sn!="" && ln!="") print sn"\t"ln
    }' > "$out"
}

GENOME_TARGET="${OUTDIR}/$(basename "${CTCF_TARGETTING%.bam}").genome"
GENOME_CONTROL="${OUTDIR}/$(basename "${CONTROL_NON_SPECIFIC%.bam}").genome"
[[ -s "$GENOME_TARGET" ]] || bam_genome_for "$CTCF_TARGETTING" "$GENOME_TARGET"
[[ -s "$GENOME_CONTROL" ]] || bam_genome_for "$CONTROL_NON_SPECIFIC" "$GENOME_CONTROL"
[[ -s "$GENOME_TARGET" ]] || { echo "Failed to create genome file: $GENOME_TARGET" >&2; exit 1; }
[[ -s "$GENOME_CONTROL" ]] || { echo "Failed to create genome file: $GENOME_CONTROL" >&2; exit 1; }

# -----------------------------
# Parse bed entries into maps
# -----------------------------
NAMES=()
PATHS=()
declare -A PATH_BY_NAME=()
for item in "${BEDS[@]}"; do
  n="${item%%=*}"
  p="${item#*=}"
  NAMES+=("$n")
  PATHS+=("$p")
  PATH_BY_NAME["$n"]="$p"
done

# Define row/col sets
CTCF_NAMES=()
RAD_NAMES=()
for n in "${NAMES[@]}"; do
  if [[ "$n" == CTCF_* ]]; then
    CTCF_NAMES+=("$n")
  elif [[ "$n" == RAD21_* || "$n" == RAD_* ]]; then
    RAD_NAMES+=("$n")
  fi
done

[[ ${#CTCF_NAMES[@]} -gt 0 ]] || { echo "No CTCF_* beds found (names must start with CTCF_)" >&2; exit 1; }
[[ ${#RAD_NAMES[@]}  -gt 0 ]] || { echo "No RAD21_* (or RAD_*) beds found (names must start with RAD21_ or RAD_)" >&2; exit 1; }

# -----------------------------
# Working directories
# -----------------------------
SORTDIR="${OUTDIR}/beds_sorted_bamorder"
WINDIR="${OUTDIR}/bed_windows_w${WINDOW}"
FBAMDIR="${OUTDIR}/filtered_bams"
mkdir -p "$SORTDIR" "$WINDIR" "$FBAMDIR"

sorted_bed_for_name_genome() {
  local name="$1"
  local path="$2"
  local genome="$3"
  local tag="$4"
  local out="${SORTDIR}/${tag}.${name}.sorted.bed"
  if [[ ! -s "$out" || "$out" -ot "$path" || "$out" -ot "$genome" ]]; then
    bedtools sort -g "$genome" -i "$path" > "$out"
  fi
  echo "$out"
}

# CTCF_i within WINDOW bp of RAD21_j (returns CTCF intervals)
window_regions() {
  local ctcf_name="$1"
  local ctcf_bed="$2"
  local rad_name="$3"
  local rad_bed="$4"
  local tag="$5"
  local out="${WINDIR}/${tag}.${ctcf_name}__within${WINDOW}bp__${rad_name}.bed"

  if [[ ! -s "$out" || "$out" -ot "$ctcf_bed" || "$out" -ot "$rad_bed" ]]; then
    bedtools window -u -w "$WINDOW" -a "$ctcf_bed" -b "$rad_bed" > "$out"
  fi
  echo "$out"
}

make_filtered_bam() {
  local inbam="$1"
  local tag="$2"
  local outbam="${FBAMDIR}/${tag}.F${FLAG_FILTER}.no${EXCLUDE_RNAME}.sorted.bam"

  if [[ -s "$outbam" && ( -s "${outbam}.bai" || -s "${outbam%.bam}.bai" ) ]]; then
    echo "$outbam"
    return
  fi

  samtools view -@ "$THREADS" -h -F "$FLAG_FILTER" "$inbam" \
  | awk -v ex="$EXCLUDE_RNAME" 'BEGIN{OFS="\t"} /^@/{print; next} $3!=ex {print}' \
  | samtools view -@ "$THREADS" -b - \
  | samtools sort -@ "$THREADS" -o "$outbam" -

  samtools index -@ "$THREADS" "$outbam"
  echo "$outbam"
}

count_sum_per_interval_filteredbam() {
  local filtered_bam="$1"
  local genome="$2"
  local regions_bed="$3"

  bedtools intersect -sorted -g "$genome" -a "$regions_bed" -b "$filtered_bam" -c \
  | awk '{s+=$NF} END{print s+0}'
}

count_matrix_for_bam() {
  local bam="$1"
  local bamname="$2"
  local genome="$3"
  local out_tsv="${OUTDIR}/counts.${bamname}.tsv"

  echo -e "bam\tctcf_set\trad21_set\toverlaps" > "$out_tsv"

  local filtered_bam
  filtered_bam="$(make_filtered_bam "$bam" "$bamname")"

  # Sort all beds once per BAM contig order
  declare -A SORTED=()
  for n in "${NAMES[@]}"; do
    SORTED["$n"]="$(sorted_bed_for_name_genome "$n" "${PATH_BY_NAME[$n]}" "$genome" "$bamname")"
  done

  # Only compute CTCF(rows) x RAD21(cols)
  for ctcf in "${CTCF_NAMES[@]}"; do
    for rad in "${RAD_NAMES[@]}"; do
      echo "Counting ${bamname}: ${ctcf} within ${WINDOW}bp of ${rad}"
      local regions
      regions="$(window_regions "$ctcf" "${SORTED[$ctcf]}" "$rad" "${SORTED[$rad]}" "$bamname")"

      # Optional warning if regions are empty (helps catch contig mismatches)
      if [[ ! -s "$regions" ]]; then
        echo "WARNING: empty regions for ${ctcf} vs ${rad} (WINDOW=${WINDOW})" >&2
      fi

      local c
      c="$(count_sum_per_interval_filteredbam "$filtered_bam" "$genome" "$regions")"
      echo -e "${bamname}\t${ctcf}\t${rad}\t${c}" >> "$out_tsv"
    done
  done
}

count_matrix_for_bam "$CTCF_TARGETTING" "CTCF_targetting" "$GENOME_TARGET"
count_matrix_for_bam "$CONTROL_NON_SPECIFIC" "control_non_specific" "$GENOME_CONTROL"

# -----------------------------
# Build matrix TSVs (rectangular)
# -----------------------------
python3 - <<'PY'
import os
import pandas as pd

outdir = os.environ["OUTDIR"]

def make_matrix(in_tsv, out_matrix_tsv):
    df = pd.read_csv(in_tsv, sep="\t")
    # preserve encountered orders
    row_order = list(dict.fromkeys(df["ctcf_set"].tolist()))
    col_order = list(dict.fromkeys(df["rad21_set"].tolist()))
    mat = df.pivot(index="ctcf_set", columns="rad21_set", values="overlaps")
    mat = mat.reindex(index=row_order, columns=col_order)
    mat.to_csv(out_matrix_tsv, sep="\t", index=True)

make_matrix(f"{outdir}/counts.CTCF_targetting.tsv", f"{outdir}/matrix.CTCF_targetting.tsv")
make_matrix(f"{outdir}/counts.control_non_specific.tsv", f"{outdir}/matrix.control_non_specific.tsv")
PY

# -----------------------------
# Plot heatmaps
# -----------------------------
python3 - <<'PY'
import os
import pandas as pd
import matplotlib.pyplot as plt

outdir = os.environ["OUTDIR"]

try:
    import seaborn as sns
except ImportError:
    raise SystemExit(
        "Missing seaborn. Install with:\n"
        "  python3 -m pip install --user seaborn\n"
        "Then rerun."
    )

sns.set_context("talk")

def plot_matrix(matrix_tsv, out_png, title):
    mat = pd.read_csv(matrix_tsv, sep="\t", index_col=0)
    plt.figure(figsize=(1.2*mat.shape[1] + 3, 1.0*mat.shape[0] + 2.5))
    ax = sns.heatmap(
        mat, annot=True, fmt="g", cmap="viridis",
        linewidths=0.5, linecolor="white",
        cbar_kws={"label":"Sum over CTCF intervals (within window) of (# overlapping reads), chrM excluded, -F 0x904"}
    )
    ax.set_xlabel("RAD21 set")
    ax.set_ylabel("CTCF set")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

plot_matrix(f"{outdir}/matrix.CTCF_targetting.tsv", f"{outdir}/heatmap.CTCF_targetting.png",
            "CTCF targetting: CTCF vs RAD21 proximity matrix")
plot_matrix(f"{outdir}/matrix.control_non_specific.tsv", f"{outdir}/heatmap.control_non_specific.png",
            "control non specific: CTCF vs RAD21 proximity matrix")
PY

echo "Done. Outputs in: $OUTDIR"
echo "Finished: $(date)"
