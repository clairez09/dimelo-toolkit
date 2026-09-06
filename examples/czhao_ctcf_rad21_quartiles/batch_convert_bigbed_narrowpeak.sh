#!/usr/bin/env bash
set -euo pipefail

OUTDIR="encode_gm12878_hg38_tfchip"   # must match what you use in the python downloader

# Convert any downloaded .bigBed -> .narrowPeak (BED text)
find "$OUTDIR/downloads" -type f -name "*.bigBed" -print0 |
while IFS= read -r -d '' bb; do
  # write next to the bigBed file
  out="${bb%.bigBed}.narrowPeak"
  bigBedToBed "$bb" "$out"
done
