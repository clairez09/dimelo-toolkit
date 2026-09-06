#!/bin/bash
# =============================================================================
# prepare_peaksets.sh -- everything between "downloaded from ENCODE" and the
# CTCF/RAD21 peak sets that pipeline.sh consumes.
#
# Reconstruction of prep that was originally done interactively and never
# written down. Verified to regenerate both peak files byte-for-byte:
#   CTCF.filteredPeaks.withSummit.bed   md5 ba105e0a2b915aa2a7e4a446b99bf0da
#   RAD21.filteredPeaks.withSummit.bed  md5 856eaede8aae8e9b395c0783d3bedc89
#
# Not covered (run before this): download_encode_gm12878_ctcf_rad21_hg38.py
#   ENCSR000DRZ (CTCF)  idr_peaks.ENCFF485TGR.bigBed / fc_over_control.ENCFF644EEX.bigWig
#   ENCSR000BMY (RAD21) idr_peaks.ENCFF101UQZ.bigBed / fc_over_control.ENCFF496LTS.bigWig
# That script also writes a metadata file per experiment recording the ENCODE
# accessions, URLs and md5s it selected.
#
# Usage:
#   bash prepare_peaksets.sh [--outdir DIR] [--scan-dir DIR] [options]
#
#   --outdir DIR         where the peak sets are written        (default: $PWD)
#   --scan-dir DIR       genome / FIMO cache                    (default: OUTDIR/scan)
#   --ctcf-bigbed FILE   ENCODE CTCF IDR peaks                  (default: OUTDIR/ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed)
#   --rad21-bigbed FILE  ENCODE RAD21 IDR peaks                 (default: OUTDIR/ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed)
#   --genome-fa FILE     NCBI GRCh38.p14 genomic FASTA          (default: SCANDIR/hg38.p14.fa)
#   --genome-ucsc FILE   renamed genome, built if absent        (default: SCANDIR/hg38.primary.ucsc.fa)
#   --chrom-map FILE     RefSeq -> UCSC name map, 2 cols        (default: SCANDIR/refseq_to_ucsc_primary.tsv)
#   --motif-meme FILE    motif in MEME format                   (default: SCANDIR/MA0139.1.meme)
#   --meme-sif FILE      MEME Suite container                   (default: SCANDIR/meme.sif)
#   --filter-script FILE summit isolation filter                (default: next to this script)
#   --isolation-dist N   drop summits with a neighbour this close (default: 1000)
#   --fimo-thresh X      FIMO p-value threshold                 (default: 1e-4)
#   --motif-shift N      0 = the historical 2 bp-left coordinates,
#                        2 = corrected coordinates              (default: 0)
#   --force              rebuild cached intermediates (genome, regions.fa, FIMO)
#
# !! MOTIF_SHIFT: the original motif BED sat 2 bp left of the true match.
# !! Default reproduces that so existing downstream results stay consistent.
# !! Changing it invalidates every downstream file -- delete them and rerun.
# =============================================================================

set -euo pipefail

# Sibling scripts live next to this file (symlinks resolved).
SELFDIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)

OUTDIR=""
SCANDIR=""
CTCF_BIGBED=""
RAD21_BIGBED=""
GENOME_FA=""
GENOME_UCSC=""
CHROM_MAP=""
MOTIF_MEME=""
MEME_SIF=""
FILTER_SCRIPT=""
ISOLATION_DIST=1000
FIMO_THRESH=1e-4
MOTIF_SHIFT=0
FORCE=0

die() { echo "ERROR: $*" >&2; exit 1; }
usage() { sed -n '2,/^# ====/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [ $# -gt 0 ]; do
  case "$1" in
    --outdir)         OUTDIR="$2";        shift 2 ;;
    --scan-dir)       SCANDIR="$2";       shift 2 ;;
    --ctcf-bigbed)    CTCF_BIGBED="$2";   shift 2 ;;
    --rad21-bigbed)   RAD21_BIGBED="$2";  shift 2 ;;
    --genome-fa)      GENOME_FA="$2";     shift 2 ;;
    --genome-ucsc)    GENOME_UCSC="$2";   shift 2 ;;
    --chrom-map)      CHROM_MAP="$2";     shift 2 ;;
    --motif-meme)     MOTIF_MEME="$2";    shift 2 ;;
    --meme-sif)       MEME_SIF="$2";      shift 2 ;;
    --filter-script)  FILTER_SCRIPT="$2"; shift 2 ;;
    --isolation-dist) ISOLATION_DIST="$2";shift 2 ;;
    --fimo-thresh)    FIMO_THRESH="$2";   shift 2 ;;
    --motif-shift)    MOTIF_SHIFT="$2";   shift 2 ;;
    --force)          FORCE=1;            shift   ;;
    -h|--help)        usage 0 ;;
    *) die "unknown argument: $1  (try --help)" ;;
  esac
done

# defaults that depend on other options
: "${OUTDIR:=$PWD}"
: "${SCANDIR:=$OUTDIR/scan}"
: "${CTCF_BIGBED:=$OUTDIR/ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed}"
: "${RAD21_BIGBED:=$OUTDIR/ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed}"
: "${GENOME_FA:=$SCANDIR/hg38.p14.fa}"
: "${GENOME_UCSC:=$SCANDIR/hg38.primary.ucsc.fa}"
: "${CHROM_MAP:=$SCANDIR/refseq_to_ucsc_primary.tsv}"
: "${MOTIF_MEME:=$SCANDIR/MA0139.1.meme}"
: "${MEME_SIF:=$SCANDIR/meme.sif}"
: "${FILTER_SCRIPT:=$SELFDIR/filter_isolated_summits.py}"

MOTIF_TAG=$(basename "$MOTIF_MEME" .meme)     # e.g. MA0139.1
FIMO_DIR=$SCANDIR/fimo_${MOTIF_TAG}

for f in "$CTCF_BIGBED" "$RAD21_BIGBED" "$CHROM_MAP" "$FILTER_SCRIPT"; do
  [ -s "$f" ] || die "missing input: $f"
done

ml biology ucsc-utils
ml biology samtools
ml biology bedtools

mkdir -p "$SCANDIR" "$OUTDIR"

# -----------------------------------------------------------------------------
# STEP 1  bigBed -> narrowPeak   (body of batch_convert_bigbed_narrowpeak.sh)
# -----------------------------------------------------------------------------
CTCF_NP="${CTCF_BIGBED%.bigBed}.narrowPeak"
RAD21_NP="${RAD21_BIGBED%.bigBed}.narrowPeak"
for pair in "$CTCF_BIGBED:$CTCF_NP" "$RAD21_BIGBED:$RAD21_NP"; do
  bb="${pair%:*}"; np="${pair#*:}"
  if [ "$FORCE" = 1 ] || [ ! -s "$np" ]; then
    echo "[prep] bigBedToBed $bb -> $np"
    bigBedToBed "$bb" "$np"
  fi
done

# -----------------------------------------------------------------------------
# STEP 2  reference: NCBI GRCh38.p14 -> primary chromosomes with UCSC names
#         GENOME_FA is the RefSeq assembly, i.e.
#         https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/\
#           GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
#         CHROM_MAP keeps chr1..chr22,chrX,chrY,chrM and drops the scaffolds.
# -----------------------------------------------------------------------------
if [ "$FORCE" = 1 ] || [ ! -s "$GENOME_UCSC" ]; then
  [ -s "$GENOME_FA" ] || die "missing genome FASTA: $GENOME_FA (see URL above)"
  [ -s "$GENOME_FA.fai" ] || samtools faidx "$GENOME_FA"
  echo "[prep] building $GENOME_UCSC"
  cut -f1 "$CHROM_MAP" \
    | xargs samtools faidx "$GENOME_FA" \
    | awk 'NR==FNR { map[$1]=$2; next }
           /^>/     { acc=substr($1,2); print ">" (acc in map ? map[acc] : acc); next }
                    { print }' "$CHROM_MAP" - \
    > "$GENOME_UCSC"
  samtools faidx "$GENOME_UCSC"
fi

# -----------------------------------------------------------------------------
# STEP 3  scan regions = the CTCF IDR peaks verbatim, plus their sequence
#         Nothing is recentered: these are ENCODE's native IDR intervals. Most
#         are 250 bp with the point source at offset 125, but ~11.5k are not.
#         -name+ writes ".::chr:start-end" headers, which is what lets FIMO
#         report genomic coordinates.
# -----------------------------------------------------------------------------
cp -f "$CTCF_NP" "$SCANDIR/regions.bed"
if [ "$FORCE" = 1 ] || [ ! -s "$SCANDIR/regions.fa" ]; then
  echo "[prep] bedtools getfasta -> $SCANDIR/regions.fa"
  bedtools getfasta -fi "$GENOME_UCSC" -bed "$SCANDIR/regions.bed" -name+ -fo "$SCANDIR/regions.fa"
fi

# -----------------------------------------------------------------------------
# STEP 4  FIMO scan for the CTCF motif (JASPAR MA0139.1, 19 bp)
#         MA0139.1.meme: https://jaspar.elixir.no/api/v1/matrix/MA0139.1.meme
#         meme.sif is MEME Suite 5.5.5 (version recorded in fimo.xml); the
#         exact `apptainer pull` that produced it was never recorded.
#         Only the PWM comes from JASPAR -- sites are called de novo here, no
#         precomputed TFBS track is used anywhere in this pipeline.
# -----------------------------------------------------------------------------
if [ "$FORCE" = 1 ] || [ ! -s "$FIMO_DIR/fimo.tsv" ]; then
  [ -s "$MEME_SIF" ]   || die "missing MEME container: $MEME_SIF (MEME Suite 5.5.5)"
  [ -s "$MOTIF_MEME" ] || die "missing motif file: $MOTIF_MEME (see URL above)"
  echo "[prep] FIMO $MOTIF_TAG (thresh $FIMO_THRESH)"
  apptainer exec --bind "$SCANDIR":/work "$MEME_SIF" \
    fimo --oc "/work/$(basename "$FIMO_DIR")" \
         --thresh "$FIMO_THRESH" \
         --parse-genomic-coord \
         "/work/$(basename "$MOTIF_MEME")" /work/regions.fa
fi

# -----------------------------------------------------------------------------
# STEP 5  fimo.tsv -> BED6 of motif hits
#         fimo.tsv cols: seq(3) start(4) stop(5) strand(6) score(7)
#         The -2/-1 is the coordinate bug; --motif-shift 2 cancels it.
#         NF>=10 skips the header and FIMO's trailing comment/blank lines. The
#         original awk did not, so the historical file carries one junk record
#         (43403 lines vs 43402 here). STEP 6 drops it either way.
# -----------------------------------------------------------------------------
HITS=$SCANDIR/${MOTIF_TAG}_hits.bed
awk -F'\t' -v s="$MOTIF_SHIFT" '
  NF>=10 && $1!="motif_id" && $1 !~ /^#/ {
    printf "%s\t%d\t%d\t%s\t%s\t%s\n", $3, $4-2+s, $5-1+s, $2, $7, $6
  }' "$FIMO_DIR/fimo.tsv" > "$HITS"
cp -f "$HITS" "$OUTDIR/$(basename "$HITS")"

# -----------------------------------------------------------------------------
# STEP 6  CTCF peak set = isolated motif hits, summit at the motif centre
#         43402 hits -> 28162. Karyotypic chrom order (chr1,chr2,...,chr22,chrX).
#         Note this drops BOTH members of any pair/cluster closer than
#         ISOLATION_DIST -- it does not keep the better-scoring one.
# -----------------------------------------------------------------------------
python3 "$FILTER_SCRIPT" \
  --in  "$HITS" \
  --out "$OUTDIR/CTCF.filteredPeaks.withSummit.bed" \
  --format motif \
  --dist "$ISOLATION_DIST" \
  --chrom-order natural

# -----------------------------------------------------------------------------
# STEP 7  RAD21 peak set = isolated IDR peaks (no motif filter)
#         44515 peaks -> 41373. Summit from narrowPeak col10, score forced to 0.
#         Lexicographic chrom order and the "CTCF_peak_" name prefix are quirks
#         of the original run, kept so the file reproduces exactly.
# -----------------------------------------------------------------------------
python3 "$FILTER_SCRIPT" \
  --in  "$RAD21_NP" \
  --out "$OUTDIR/RAD21.filteredPeaks.withSummit.bed" \
  --format narrowpeak \
  --dist "$ISOLATION_DIST" \
  --chrom-order lexicographic \
  --prefix CTCF_peak

echo "[prep] done -> $OUTDIR/{CTCF,RAD21}.filteredPeaks.withSummit.bed"
