#!/bin/bash
# =============================================================================
# CTCF / RAD21 DiMeLo pipeline.
#
# Self-contained: every helper script it calls lives next to this file, and no
# path outside this directory is hard-coded. Reference data (genome, chain,
# BAMs) and the run directory are all options -- see below.
#
# Typical use:
#   mkdir -p myrun/logs && cd myrun
#   sbatch /path/to/pipeline_scripts/example_run.sbatch
#
# Peak-set prep (ENCODE bigBed -> FIMO -> CTCF/RAD21 filteredPeaks) lives in
# prepare_peaksets.sh and is the first step below.
#
#   --outdir DIR         working/output directory        (default: $PWD)
#   --scan-dir DIR       genome / FIMO cache             (default: OUTDIR/scan)
#   --scripts-dir DIR    where the helper scripts live   (default: this file's
#                        directory, resolved through symlinks)
#
#   ENCODE inputs
#   --ctcf-bigbed FILE   CTCF IDR peaks   (ENCSR000DRZ / ENCFF485TGR)
#   --ctcf-bw FILE       CTCF fold-change over control   (ENCFF644EEX)
#   --rad21-bigbed FILE  RAD21 IDR peaks  (ENCSR000BMY / ENCFF101UQZ)
#   --rad21-bw FILE      RAD21 fold-change over control  (ENCFF496LTS)
#
#   references
#   --genome-fa FILE     NCBI GRCh38.p14 genomic FASTA
#   --genome-ucsc FILE   renamed primary-chromosome genome, built if absent
#   --chrom-map FILE     RefSeq -> UCSC chromosome name map
#   --motif-meme FILE    CTCF motif, MEME format (JASPAR MA0139.1)
#   --meme-sif FILE      MEME Suite 5.5.5 container
#   --chain FILE         hg38 -> chm13v2 liftOver chain
#   --chm13-fasta FILE   chm13v2.0 FASTA for the pileups
#
#   DiMeLo reads
#   --target-bam FILE    CTCF-targeting sample     (barcode17)
#   --control-bam FILE   non-specific control      (barcode18)
#
#   tuning
#   --cores N            (default: 16)
#   --window N           pileup half-window, bp     (default: 2000)
#   --isolation-dist N   summit isolation distance  (default: 1000)
#   --fimo-thresh X      FIMO p-value threshold     (default: 1e-4)
#   --motif-shift N      0 = the historical 2 bp-left motif coordinates,
#                        2 = corrected                (default: 0)
# =============================================================================

set -euo pipefail   # NOTE: the original had no error checking; drop this if a
                    # step you expect to fail is now aborting the run.

# Helper scripts live next to this file (symlinks resolved), overridable with
# --scripts-dir.
SCRIPTS=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)

OAK=/oak/stanford/groups/altemose

OUTDIR="$PWD"
SCANDIR=""
CTCF_BIGBED=""
RAD21_BIGBED=""
CTCF_BW=""
RAD21_BW=""
GENOME_FA=""
GENOME_UCSC=""
CHROM_MAP=""
MOTIF_MEME=""
MEME_SIF=""
CHAIN=$OAK/references/chains/hg38-chm13v2.over.chain
CHM13_FASTA=$OAK/references/fastas/chm13v2.0.fasta
TARGET_BAM=$OAK/data/20250110_NG_one_pot/barcode17.merged.sorted.bam
CONTROL_BAM=$OAK/data/20250110_NG_one_pot/barcode18.merged.sorted.bam
CORES=16
WINDOW=2000
ISOLATION_DIST=1000
FIMO_THRESH=1e-4
MOTIF_SHIFT=0

usage() { sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --outdir)         OUTDIR="$2";         shift 2 ;;
    --scan-dir)       SCANDIR="$2";        shift 2 ;;
    --scripts-dir)    SCRIPTS="$2";        shift 2 ;;
    --ctcf-bigbed)    CTCF_BIGBED="$2";    shift 2 ;;
    --ctcf-bw)        CTCF_BW="$2";        shift 2 ;;
    --rad21-bigbed)   RAD21_BIGBED="$2";   shift 2 ;;
    --rad21-bw)       RAD21_BW="$2";       shift 2 ;;
    --genome-fa)      GENOME_FA="$2";      shift 2 ;;
    --genome-ucsc)    GENOME_UCSC="$2";    shift 2 ;;
    --chrom-map)      CHROM_MAP="$2";      shift 2 ;;
    --motif-meme)     MOTIF_MEME="$2";     shift 2 ;;
    --meme-sif)       MEME_SIF="$2";       shift 2 ;;
    --chain)          CHAIN="$2";          shift 2 ;;
    --chm13-fasta)    CHM13_FASTA="$2";    shift 2 ;;
    --target-bam)     TARGET_BAM="$2";     shift 2 ;;
    --control-bam)    CONTROL_BAM="$2";    shift 2 ;;
    --cores)          CORES="$2";          shift 2 ;;
    --window)         WINDOW="$2";         shift 2 ;;
    --isolation-dist) ISOLATION_DIST="$2"; shift 2 ;;
    --fimo-thresh)    FIMO_THRESH="$2";    shift 2 ;;
    --motif-shift)    MOTIF_SHIFT="$2";    shift 2 ;;
    -h|--help)        usage 0 ;;
    *) echo "ERROR: unknown argument: $1  (try --help)" >&2; exit 1 ;;
  esac
done

OUTDIR=$(cd "$OUTDIR" && pwd)
: "${SCANDIR:=$OUTDIR/scan}"
: "${CTCF_BIGBED:=$OUTDIR/ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed}"
: "${RAD21_BIGBED:=$OUTDIR/ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed}"
: "${CTCF_BW:=$OUTDIR/ENCSR000DRZ/fc_over_control.ENCFF644EEX.bigWig}"
: "${RAD21_BW:=$OUTDIR/ENCSR000BMY/fc_over_control.ENCFF496LTS.bigWig}"
: "${GENOME_FA:=$SCANDIR/hg38.p14.fa}"
: "${GENOME_UCSC:=$SCANDIR/hg38.primary.ucsc.fa}"
: "${CHROM_MAP:=$SCANDIR/refseq_to_ucsc_primary.tsv}"
: "${MOTIF_MEME:=$SCANDIR/MA0139.1.meme}"
: "${MEME_SIF:=$SCANDIR/meme.sif}"

ml biology ucsc-utils
ml biology samtools
ml biology bedtools

cd "$OUTDIR"

#Build CTCF/RAD21 peak sets from the ENCODE downloads
bash $SCRIPTS/prepare_peaksets.sh \
  --outdir "$OUTDIR" --scan-dir "$SCANDIR" \
  --ctcf-bigbed "$CTCF_BIGBED" --rad21-bigbed "$RAD21_BIGBED" \
  --genome-fa "$GENOME_FA" --genome-ucsc "$GENOME_UCSC" --chrom-map "$CHROM_MAP" \
  --motif-meme "$MOTIF_MEME" --meme-sif "$MEME_SIF" \
  --isolation-dist "$ISOLATION_DIST" --fimo-thresh "$FIMO_THRESH" \
  --motif-shift "$MOTIF_SHIFT"
python3 $SCRIPTS/score_summits_by_peakmean_new.py \
  --ctcf-peaks-bed CTCF.filteredPeaks.withSummit.bed --ctcf-bw "$CTCF_BW" \
  --rad21-peaks-bed RAD21.filteredPeaks.withSummit.bed --rad21-bw "$RAD21_BW" \
  --outdir ./
#Convert via liftOver
python3 $SCRIPTS/liftover_both.py --chain "$CHAIN"
#Split into quartiles
bash $SCRIPTS/split_into_quartiles.sh CTCF.chm13.median.bed CTCF.chm13.median
bash $SCRIPTS/split_into_quartiles.sh RAD21.chm13.median.bed RAD21.chm13.median
#Pileup & Plot
python3 $SCRIPTS/m6a_rad21_ctcf_quartiles.py \
  --CTCF_targetting "$TARGET_BAM" \
  --control_non_specific "$CONTROL_BAM" \
  --fasta "$CHM13_FASTA" \
  --rad21_q1 RAD21.chm13.median.Q1.bed \
  --rad21_q2 RAD21.chm13.median.Q2.bed \
  --rad21_q3 RAD21.chm13.median.Q3.bed \
  --rad21_q4 RAD21.chm13.median.Q4.bed \
  --ctcf_q1  CTCF.chm13.median.Q1.bed \
  --ctcf_q2  CTCF.chm13.median.Q2.bed \
  --ctcf_q3  CTCF.chm13.median.Q3.bed \
  --ctcf_q4  CTCF.chm13.median.Q4.bed \
  --outdir ./m6A_profiles_chm13 \
  --cores "$CORES"
python3 $SCRIPTS/m6a_bgsub_from_plots.py \
  --pileup_root ./m6A_profiles_chm13/pileups \
  --outdir ./m6A_profiles_chm13/bgsub_plots \
  --rad21_q1 RAD21.chm13.median.Q1.bed \
  --rad21_q2 RAD21.chm13.median.Q2.bed \
  --rad21_q3 RAD21.chm13.median.Q3.bed \
  --rad21_q4 RAD21.chm13.median.Q4.bed \
  --ctcf_q1  CTCF.chm13.median.Q1.bed \
  --ctcf_q2  CTCF.chm13.median.Q2.bed \
  --ctcf_q3  CTCF.chm13.median.Q3.bed \
  --ctcf_q4  CTCF.chm13.median.Q4.bed \
  --window "$WINDOW" \
  --cores "$CORES" \
  --ymin -0.1 --ymax 0.45
#Intersect RAD21 with CTCF
bash $SCRIPTS/intersect_quartiles_symmetric.sh
#Pilup & Plot
python3 $SCRIPTS/m6a_pileup_4sets.py \
  --CTCF_targetting "$TARGET_BAM" \
  --control_non_specific "$CONTROL_BAM" \
  --fasta "$CHM13_FASTA" \
  --beds \
    HH=./intersections_150bp/CTCFhigh_RAD21high.w150.bed \
    HL=./intersections_150bp/CTCFhigh_RAD21low.w150.bed \
    LH=./intersections_150bp/CTCFlow_RAD21high.w150.bed \
    LL=./intersections_150bp/CTCFlow_RAD21low.w150.bed \
  --outdir ./intersection_profiles_w150 \
  --cores "$CORES" \
  --window "$WINDOW" \
  --skip_existing_pileups \
  --ymax 0.45
python3 $SCRIPTS/m6a_bgsub_4sets_from_plots.py \
  --pileup_root ./intersection_profiles_w150/pileups \
  --outdir ./intersection_profiles_w150/bgsub \
  --beds \
    HH=./intersections_150bp/CTCFhigh_RAD21high.w150.bed \
    HL=./intersections_150bp/CTCFhigh_RAD21low.w150.bed \
    LH=./intersections_150bp/CTCFlow_RAD21high.w150.bed \
    LL=./intersections_150bp/CTCFlow_RAD21low.w150.bed \
  --window "$WINDOW" \
  --cores "$CORES" \
  --ymin -0.1 --ymax 0.35
#Same number of reads pileup & plot
python3 $SCRIPTS/m6a_pileup_4sets_subsampled.py \
  --CTCF_targetting "$TARGET_BAM" \
  --control_non_specific "$CONTROL_BAM" \
  --fasta "$CHM13_FASTA" \
  --beds \
    HH=./intersections_150bp/CTCFhigh_RAD21high.w150.bed \
    HL=./intersections_150bp/CTCFhigh_RAD21low.w150.bed \
    LH=./intersections_150bp/CTCFlow_RAD21high.w150.bed \
    LL=./intersections_150bp/CTCFlow_RAD21low.w150.bed \
  --outdir ./intersection_profiles_w150_equalized \
  --ymax 0.45
python3 $SCRIPTS/m6a_bgsub_4sets_from_plots_subsampled.py \
  --outdir ./intersection_profiles_w150_subsampled/bgsub \
  --CTCF_targetting "$TARGET_BAM" \
  --control_non_specific "$CONTROL_BAM" \
  --fasta "$CHM13_FASTA" \
  --beds \
    HH=./intersections_150bp/CTCFhigh_RAD21high.w150.bed \
    HL=./intersections_150bp/CTCFhigh_RAD21low.w150.bed \
    LH=./intersections_150bp/CTCFlow_RAD21high.w150.bed \
    LL=./intersections_150bp/CTCFlow_RAD21low.w150.bed
#Confusion Matricies
bash $SCRIPTS/make_confusion_and_plot.sh
bash $SCRIPTS/rad_ctcf_quartile_confmat.sh \
  --CTCF_targetting "$TARGET_BAM" \
  --control_non_specific "$CONTROL_BAM"
