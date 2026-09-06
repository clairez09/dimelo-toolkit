#!/usr/bin/env python3
import argparse
import math
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple


def read_bed_with_summit_and_strand(path: Path) -> List[Tuple[str, int, int, str, str, int]]:
    """
    Read a BED-like file that contains peak interval + strand + summit.

    Expected columns (>=7):
      chrom  start  end  name  score  strand  summit

    Notes:
    - "score" is read only to keep BED column positions; it is not used.
    - summit is assumed to be a 0-based coordinate (as typical in BED-style outputs).
    """
    out: List[Tuple[str, int, int, str, str, int]] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                raise SystemExit(
                    f"Expected >=7 columns (chrom,start,end,name,score,strand,summit) in {path}"
                )
            chrom = parts[0]
            peak_start = int(parts[1])
            peak_end = int(parts[2])
            name = parts[3]
            strand = parts[5]
            summit = int(parts[6])
            out.append((chrom, peak_start, peak_end, name, strand, summit))
    return out


def run_bigwig_average_over_bed(bigwig: Path, bed: Path, out_tab: Path) -> None:
    subprocess.check_call(["bigWigAverageOverBed", str(bigwig), str(bed), str(out_tab)])


def parse_ucsc_tab_mean(tab_path: Path) -> Dict[str, float]:
    """
    Parse bigWigAverageOverBed output.
    Returns dict: name -> mean (column 6; 0-based index 5).
    """
    out: Dict[str, float] = {}
    with open(tab_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            name = parts[0]
            try:
                out[name] = float(parts[5])
            except ValueError:
                out[name] = float("nan")
    return out


def write_summit_bed_scored_by_peak_mean_stranded(
    peaks_bed_with_summit: Path,
    bigwig: Path,
    out_bed: Path
) -> None:
    """
    Input (>=7 cols): chrom start end name score strand summit
    Computes mean bigWig signal over [start,end) for each peak.
    Output (BED6): chrom summit summit+1 name mean strand
    """
    peaks = read_bed_with_summit_and_strand(peaks_bed_with_summit)
    if not peaks:
        raise SystemExit(f"No peaks read from: {peaks_bed_with_summit}")

    # NOTE: If 'name' is not unique, values will be overwritten in mean_by_name.
    with tempfile.TemporaryDirectory(prefix="bwavg_") as td:
        td = Path(td)
        tmp_bed = td / "peaks_for_bwavg.bed"
        out_tab = td / "avg.tab"

        # Use full peak intervals for averaging; keep 'name' for joining.
        with open(tmp_bed, "w") as w:
            for chrom, s, e, name, strand, summit in peaks:
                w.write(f"{chrom}\t{s}\t{e}\t{name}\n")

        run_bigwig_average_over_bed(bigwig, tmp_bed, out_tab)
        mean_by_name = parse_ucsc_tab_mean(out_tab)

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(out_bed, "w") as w:
        for chrom, peak_s, peak_e, name, strand, summit in peaks:
            m = mean_by_name.get(name, float("nan"))
            score = 0.0 if (m is None or math.isnan(m)) else m
            w.write(f"{chrom}\t{summit}\t{summit+1}\t{name}\t{score}\t{strand}\n")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Score 1bp summit BED entries by mean bigWig signal over their peak intervals, "
            "while preserving strand from the input."
        )
    )
    ap.add_argument(
        "--ctcf-peaks-bed", required=True,
        help="Input peaks (>=7 cols): chrom start end name score strand summit"
    )
    ap.add_argument("--ctcf-bw", required=True, help="CTCF bigWig")
    ap.add_argument(
        "--rad21-peaks-bed", required=True,
        help="Input peaks (>=7 cols): chrom start end name score strand summit"
    )
    ap.add_argument("--rad21-bw", required=True, help="RAD21 bigWig")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    write_summit_bed_scored_by_peak_mean_stranded(
        Path(args.ctcf_peaks_bed), Path(args.ctcf_bw), outdir / "CTCF.median.bed"
    )
    write_summit_bed_scored_by_peak_mean_stranded(
        Path(args.rad21_peaks_bed), Path(args.rad21_bw), outdir / "RAD21.median.bed"
    )


if __name__ == "__main__":
    main()
