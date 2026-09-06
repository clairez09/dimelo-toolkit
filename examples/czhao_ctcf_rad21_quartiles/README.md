# CTCF / RAD21 DiMeLo pipeline

Takes ENCODE GM12878 CTCF and RAD21 ChIP-seq, builds isolated motif-anchored
site sets, lifts them to CHM13v2, splits them into ChIP-signal quartiles, and
produces m6A pileup profiles from DiMeLo-seq reads around each set.

One command runs all of it. This directory is self-contained — every script it
calls lives here, and every path it touches is an option.

---

## 1. Requirements

Sherlock, with these modules (the scripts load them themselves):
`biology/ucsc-utils`, `biology/samtools`, `biology/bedtools`, `apptainer`.

Run it through Slurm. The genome build, FIMO scan and pileups are all far too
heavy for a login node.

## 2. Set up a run directory

```bash
mkdir -p $SCRATCH/myrun/logs && cd $SCRATCH/myrun
```

`logs/` must exist before you submit — the sbatch log paths are relative and
Slurm refuses the job if the directory is missing.

Then put the four ENCODE files in place (real files or symlinks):

```
ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed        CTCF IDR peaks
ENCSR000DRZ/fc_over_control.ENCFF644EEX.bigWig  CTCF fold-change
ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed        RAD21 IDR peaks
ENCSR000BMY/fc_over_control.ENCFF496LTS.bigWig  RAD21 fold-change
```

To fetch them fresh:

```bash
python3 /path/to/pipeline_scripts/download_encode_gm12878_ctcf_rad21_hg38.py \
  --outdir ./encode
```

(`--targets CTCF RAD21`, `--overwrite`, `--no-md5` are the other options. It
writes a metadata file per experiment with the accessions, URLs and md5s.)

And a `scan/` directory with the reference/motif files:

```
scan/hg38.p14.fa                     NCBI GCF_000001405.40_GRCh38.p14_genomic.fna.gz
scan/refseq_to_ucsc_primary.tsv      RefSeq -> UCSC chromosome name map
scan/MA0139.1.meme                   https://jaspar.elixir.no/api/v1/matrix/MA0139.1.meme
scan/meme.sif                        MEME Suite 5.5.5 apptainer container
```

**Symlink an existing `scan/` if you have one.** It is a cache: the renamed
genome (`hg38.primary.ucsc.fa`, ~3 GB), `regions.fa`, and the FIMO output are
only rebuilt when missing. Reusing it saves the genome build and the scan. The
cached FIMO output is valid for either `--motif-shift` value, since the shift
is applied afterwards.

Those four files are not downloaded automatically. If one is missing the
pipeline stops and prints where to get it.

## 3. Run it

```bash
sbatch /path/to/pipeline_scripts/example_run.sbatch
```

Submit from inside the run directory — results go wherever you submitted from.

Two wrappers, identical apart from the motif coordinates:

| | `--motif-shift` | use when |
|---|---|---|
| `example_run.sbatch` | 2 | normal use — motif intervals at their true position |
| `pipeline.sbatch` | 0 | matching results produced by an earlier run |

Both request 16 CPUs, 128 GB, 24 h, on `altemose,owners,normal`.

Check on it:

```bash
squeue --me
tail -f logs/example_run.*.out
seff <jobid>          # after it finishes
```

## 4. Change settings

Override any variable on the submit line — no file editing:

```bash
sbatch --export=ALL,CORES=32,WINDOW=1000 .../example_run.sbatch
```

Every option below is settable this way (uppercase, e.g. `--motif-shift` →
`MOTIF_SHIFT`). Bump `--cpus-per-task` alongside `CORES` if you raise it.

Or call `pipeline.sh` directly from an interactive allocation:

```bash
sh_dev -c 16
bash /path/to/pipeline_scripts/pipeline.sh --help
bash /path/to/pipeline_scripts/pipeline.sh --outdir . --motif-shift 2
```

### Options

| Option | Env var | Default |
|---|---|---|
| `--outdir` | `RUN` | submit directory |
| `--scan-dir` | `SCANDIR` | `OUTDIR/scan` |
| `--scripts-dir` | `SCRIPTS` | this directory |
| `--ctcf-bigbed` | `CTCF_BIGBED` | `OUTDIR/ENCSR000DRZ/idr_peaks.ENCFF485TGR.bigBed` |
| `--ctcf-bw` | `CTCF_BW` | `OUTDIR/ENCSR000DRZ/fc_over_control.ENCFF644EEX.bigWig` |
| `--rad21-bigbed` | `RAD21_BIGBED` | `OUTDIR/ENCSR000BMY/idr_peaks.ENCFF101UQZ.bigBed` |
| `--rad21-bw` | `RAD21_BW` | `OUTDIR/ENCSR000BMY/fc_over_control.ENCFF496LTS.bigWig` |
| `--genome-fa` | `GENOME_FA` | `SCANDIR/hg38.p14.fa` |
| `--genome-ucsc` | `GENOME_UCSC` | `SCANDIR/hg38.primary.ucsc.fa` (built if absent) |
| `--chrom-map` | `CHROM_MAP` | `SCANDIR/refseq_to_ucsc_primary.tsv` |
| `--motif-meme` | `MOTIF_MEME` | `SCANDIR/MA0139.1.meme` |
| `--meme-sif` | `MEME_SIF` | `SCANDIR/meme.sif` |
| `--chain` | `CHAIN` | `$OAK/references/chains/hg38-chm13v2.over.chain` |
| `--chm13-fasta` | `CHM13_FASTA` | `$OAK/references/fastas/chm13v2.0.fasta` |
| `--target-bam` | `TARGET_BAM` | `$OAK/data/20250110_NG_one_pot/barcode17.merged.sorted.bam` |
| `--control-bam` | `CONTROL_BAM` | `$OAK/data/20250110_NG_one_pot/barcode18.merged.sorted.bam` |
| `--cores` | `CORES` | `$SLURM_CPUS_PER_TASK`, else 16 |
| `--window` | `WINDOW` | 2000 (pileup half-window, bp) |
| `--isolation-dist` | `ISOLATION_DIST` | 1000 (summit isolation distance, bp) |
| `--fimo-thresh` | `FIMO_THRESH` | 1e-4 |
| `--motif-shift` | `MOTIF_SHIFT` | 0 in `pipeline.sbatch`, 2 in `example_run.sbatch` |

To run on your own DiMeLo data, the two you almost always change are
`--target-bam` and `--control-bam`.

---

## What it does

```
ENCODE GM12878 ChIP-seq
  ENCSR000DRZ (CTCF)   idr_peaks.bigBed + fc_over_control.bigWig
  ENCSR000BMY (RAD21)  idr_peaks.bigBed + fc_over_control.bigWig
         |
         |  [1] bigBedToBed
         v
  narrowPeak  (CTCF 42,381 IDR peaks / RAD21 44,515)
         |
   CTCF  |                                      RAD21
         v                                        |
  [2] hg38.p14 -> primary chroms, UCSC names      |
  [3] bedtools getfasta over the peaks            |
  [4] FIMO scan, JASPAR MA0139.1, p<1e-4          |
         v                                        |
  43,402 motif hits                               |
         |                                        |
         |  [5] fimo.tsv -> BED6                  |
         |  [6] isolation filter (1 kb)           |  [7] isolation filter (1 kb)
         v                                        v
  CTCF.filteredPeaks.withSummit.bed        RAD21.filteredPeaks.withSummit.bed
  28,162 sites, summit = motif centre      41,373 sites, summit = IDR point source
         \                                        /
          \        [8] mean ChIP signal over each peak interval
           v      (bigWigAverageOverBed) -> 1 bp summit BED
        CTCF.median.bed / RAD21.median.bed
                          |
                          |  [9] liftOver hg38 -> chm13v2
                          v
        CTCF.chm13.median.bed / RAD21.chm13.median.bed
                          |
                          |  [10] rank quartiles by score
                          v
                    *.chm13.median.Q1..Q4.bed
                          |
        [11] m6A pileups (target vs control BAM) per quartile
        [12] background-subtracted profiles
        [13] CTCF x RAD21 intersections (±150 bp) -> HH/HL/LH/LL
        [14] pileups + bgsub for the 4 sets, plain and read-count-equalised
        [15] 4x4 quartile confusion matrices + heatmaps
```

The site counts above are what a default run on the stock ENCODE inputs
produces — use them to sanity-check yours.

Details worth knowing when reading the outputs:

- **Isolation filter (steps 6/7).** A summit is dropped if *another summit lies
  within `--isolation-dist` on either side*. It drops **both** members of a
  close pair, not the weaker one; a cluster of three loses all three. At the
  1 kb default, 15,240 of 43,402 CTCF motif hits go.
- **CTCF summit** is the motif centre (`start + 9`) and keeps the FIMO score.
  **RAD21 summit** is `start + narrowPeak col10` and its score is forced to 0.
- **Scoring (step 8) is a mean, not a median**, despite the `*.median.bed`
  filenames. It is the mean fold-change-over-control across the full peak
  interval. ENCODE's own `score` column is capped at 1000 for 34,378 peaks, so
  it can't be used for ranking.
- **Quartiles (step 10) split by rank, not by score value.** Q1 lowest,
  Q4 highest.
- **Only the PWM comes from JASPAR.** Sites are called de novo by FIMO; no
  precomputed TFBS track is used.
- **RAD21 sites are named `CTCF_peak_N`.** Cosmetic, kept so files stay
  comparable across runs.

## Outputs

All written into the run directory:

```
CTCF.filteredPeaks.withSummit.bed    isolated CTCF motif sites (BED6+summit)
RAD21.filteredPeaks.withSummit.bed   isolated RAD21 peaks
MA0139.1_hits.bed                    all FIMO motif hits
CTCF.median.bed / RAD21.median.bed   1 bp summits scored by mean ChIP signal
*.chm13.median.bed                   lifted to CHM13v2
*.chm13.median.Q1..Q4.bed            signal quartiles
*.chm13.unmapped.bed                 liftOver failures (~1-2 kB, expected)
m6A_profiles_chm13/                  per-quartile pileups + plots
  └── bgsub_plots/                   background-subtracted
intersections_150bp/                 HH / HL / LH / LL site sets
intersection_profiles_w150/          pileups for the four sets
intersection_profiles_w150_equalized/  same, read-count matched
*__CTCFcount__byRAD__*.tsv/.png/.pdf   quartile confusion matrix + heatmaps
scan/                                genome + FIMO cache (reusable)
```

## Running one step

Each script is standalone and takes `--help` or reads environment variables.
Useful when you only want to redo part of a run — run them from inside the run
directory, in pipeline order:

```bash
bash prepare_peaksets.sh --outdir . --scan-dir ./scan --motif-shift 2
bash prepare_peaksets.sh --force ...          # rebuild the genome/FIMO cache
bash split_into_quartiles.sh CTCF.chm13.median.bed CTCF.chm13.median
W=300 bash intersect_quartiles_symmetric.sh   # different intersection window
bash make_confusion_and_plot.sh --help
```

Changing `--motif-shift` or `--isolation-dist` invalidates everything
downstream — delete the affected outputs and rerun from `prepare_peaksets.sh`.

## Troubleshooting

- **Job rejected immediately** — `logs/` doesn't exist in the submit directory.
- **Stops with a URL** — one of `hg38.p14.fa`, `MA0139.1.meme`, `meme.sif` is
  missing from `scan/`. Fetch it from the printed source; the pipeline won't
  guess.
- **Run aborts at a step you expected to fail** — `set -euo pipefail` is on in
  `pipeline.sh` and `prepare_peaksets.sh`. Drop that line to continue past
  errors.
- **A shared `scan/` got clobbered** — a real run overwrites `regions.bed` and
  `<motif>_hits.bed` inside `SCANDIR`. Point `--scan-dir` at a private
  directory and symlink the cached genome and FIMO output into it.
- **Half the CPUs idle** — `--cores` and `--cpus-per-task` are independent. The
  wrappers default `CORES` to `$SLURM_CPUS_PER_TASK`, but check
  `m6a_pileup_4sets_subsampled.py` before raising it, in case its subsampling
  draws depend on worker count.
- **Run went wrong after you edited a script mid-job** — bash reads scripts
  incrementally, so editing in place corrupts a running job. Write a new file
  and `mv` it into place; the rename is atomic and the running job keeps its
  original inode.
- **`SCRIPTS` resolved wrong** — pass
  `--export=ALL,SCRIPTS=/path/to/pipeline_scripts`.

## Files

**Pipeline**
| | |
|---|---|
| `pipeline.sh` | the pipeline; every path is an option (`--help`) |
| `prepare_peaksets.sh` | steps 1–7, ENCODE downloads → peak sets |
| `filter_isolated_summits.py` | the summit isolation filter |
| `example_run.sbatch` | Slurm wrapper, `--motif-shift 2` |
| `pipeline.sbatch` | Slurm wrapper, `--motif-shift 0` |

**Analysis** — `score_summits_by_peakmean_new.py`, `liftover_both.py`,
`split_into_quartiles.sh`, `m6a_rad21_ctcf_quartiles.py`,
`m6a_bgsub_from_plots.py`, `intersect_quartiles_symmetric.sh`,
`m6a_pileup_4sets.py`, `m6a_bgsub_4sets_from_plots.py`,
`m6a_pileup_4sets_subsampled.py`, `m6a_bgsub_4sets_from_plots_subsampled.py`,
`make_confusion_and_plot.sh`, `rad_ctcf_quartile_confmat.sh`,
`plot_heatmap_from_tsv.py`

**Upstream / reference** — `download_encode_gm12878_ctcf_rad21_hg38.py`,
`batch_convert_bigbed_narrowpeak.sh`, `chm13v2.chrom.sizes`

Scripts resolve their own directory through `readlink -f`, so the bundle works
from any location and through symlinks.
