#!/bin/bash
#SBATCH --job-name=example_run
#SBATCH --time=24:00:00
#SBATCH --partition=altemose,owners,normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --output=logs/example_run.%j.out
#SBATCH --error=logs/example_run.%j.err

# =============================================================================
# Sample run of the full pipeline, WITH the motif coordinate fix.
#
# Run it from the directory you want the results in:
#
#   mkdir -p myrun/logs && cd myrun
#   sbatch /path/to/pipeline_scripts/example_run.sbatch
#
# The log paths above are relative, so they land in myrun/logs/ -- create that
# directory first or Slurm will refuse the job.
#
# myrun/ needs the ENCODE downloads (real files or symlinks):
#   ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed
#   ENCSR000DRZ/fc_over_control.ENCFF644EEX.bigWig
#   ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed
#   ENCSR000BMY/fc_over_control.ENCFF496LTS.bigWig
# and a scan/ directory holding (or symlinking) the genome and motif files.
# Anything missing from scan/ is rebuilt, which is slow: the genome rename is
# ~3 GB and the FIMO scan takes a few minutes. Symlink an existing scan/ to
# skip both -- the cached FIMO output is coordinate-independent, since the
# shift is applied afterwards when fimo.tsv is converted to BED.
#
# Override anything on the sbatch command line, e.g.
#   sbatch --export=ALL,CORES=32,MOTIF_SHIFT=0 .../example_run.sbatch
# =============================================================================

set -euo pipefail

# Where this bundle lives. Derived from the path sbatch was given; override
# with --export=ALL,SCRIPTS=/path/to/pipeline_scripts if that ever fails.
if [ -z "${SCRIPTS:-}" ]; then
  SCRIPTS=$(dirname "$(scontrol show job "$SLURM_JOB_ID" | sed -n 's/^ *Command=//p' | awk '{print $1}')")
fi
SCRIPTS=$(cd "$SCRIPTS" && pwd)

# Run directory = wherever you submitted from.
RUN="${RUN:-${SLURM_SUBMIT_DIR:-$PWD}}"

OAK=/oak/stanford/groups/altemose

# ENCODE downloads (download_encode_gm12878_ctcf_rad21_hg38.py)
CTCF_BIGBED="${CTCF_BIGBED:-$RUN/ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed}"    # ENCSR000DRZ
CTCF_BW="${CTCF_BW:-$RUN/ENCSR000DRZ/fc_over_control.ENCFF644EEX.bigWig}"
RAD21_BIGBED="${RAD21_BIGBED:-$RUN/ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed}"  # ENCSR000BMY
RAD21_BW="${RAD21_BW:-$RUN/ENCSR000BMY/fc_over_control.ENCFF496LTS.bigWig}"

# references
SCANDIR="${SCANDIR:-$RUN/scan}"
GENOME_FA="${GENOME_FA:-$SCANDIR/hg38.p14.fa}"                 # NCBI GRCh38.p14
GENOME_UCSC="${GENOME_UCSC:-$SCANDIR/hg38.primary.ucsc.fa}"    # built from the above if absent
CHROM_MAP="${CHROM_MAP:-$SCANDIR/refseq_to_ucsc_primary.tsv}"
MOTIF_MEME="${MOTIF_MEME:-$SCANDIR/MA0139.1.meme}"             # JASPAR CTCF PWM, 19 bp
MEME_SIF="${MEME_SIF:-$SCANDIR/meme.sif}"                      # MEME Suite 5.5.5
CHAIN="${CHAIN:-$OAK/references/chains/hg38-chm13v2.over.chain}"
CHM13_FASTA="${CHM13_FASTA:-$OAK/references/fastas/chm13v2.0.fasta}"

# DiMeLo reads
TARGET_BAM="${TARGET_BAM:-$OAK/data/20250110_NG_one_pot/barcode17.merged.sorted.bam}"
CONTROL_BAM="${CONTROL_BAM:-$OAK/data/20250110_NG_one_pot/barcode18.merged.sorted.bam}"

# tuning
CORES="${CORES:-${SLURM_CPUS_PER_TASK:-16}}"
WINDOW="${WINDOW:-2000}"
ISOLATION_DIST="${ISOLATION_DIST:-1000}"
FIMO_THRESH="${FIMO_THRESH:-1e-4}"
MOTIF_SHIFT="${MOTIF_SHIFT:-2}"   # <-- the fix: motif intervals at their true position

bash "$SCRIPTS/pipeline.sh" \
  --outdir       "$RUN" \
  --scan-dir     "$SCANDIR" \
  --scripts-dir  "$SCRIPTS" \
  --ctcf-bigbed  "$CTCF_BIGBED" \
  --ctcf-bw      "$CTCF_BW" \
  --rad21-bigbed "$RAD21_BIGBED" \
  --rad21-bw     "$RAD21_BW" \
  --genome-fa    "$GENOME_FA" \
  --genome-ucsc  "$GENOME_UCSC" \
  --chrom-map    "$CHROM_MAP" \
  --motif-meme   "$MOTIF_MEME" \
  --meme-sif     "$MEME_SIF" \
  --chain        "$CHAIN" \
  --chm13-fasta  "$CHM13_FASTA" \
  --target-bam   "$TARGET_BAM" \
  --control-bam  "$CONTROL_BAM" \
  --cores        "$CORES" \
  --window       "$WINDOW" \
  --isolation-dist "$ISOLATION_DIST" \
  --fimo-thresh  "$FIMO_THRESH" \
  --motif-shift  "$MOTIF_SHIFT"
