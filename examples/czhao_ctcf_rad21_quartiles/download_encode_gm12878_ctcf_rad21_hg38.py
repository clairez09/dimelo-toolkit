#!/usr/bin/env python3
"""
Download released ENCODE GM12878 TF ChIP-seq (hg38/GRCh38) analysis files for:
  - CTCF
  - RAD21

For each matching experiment, download:
  - IDR thresholded peaks (narrowPeak or bigBed)
  - fold-change-over-control bigWig

Selection policy (when multiple candidates exist for an experiment):
  1) Prefer ⭐ files (ENCODE "preferred_default", then "default")
  2) Prefer newer files by date_released/date_created
  3) Peaks: prefer narrowPeak over bigBed, then larger file_size
     FC:    prefer larger file_size

Outputs:
  - downloads/  (files)
  - metadata/   (per-experiment metadata text)
  - a combined TSV summary

Requirements: Python 3.9+ (standard library only)
Usage:
  python download_encode_gm12878_ctcf_rad21_hg38.py --outdir encode_out
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ENCODE = "https://www.encodeproject.org"
UA = "encode-downloader/1.2 (urllib; GM12878 CTCF/RAD21 GRCh38; prefer starred+newer)"


@dataclass
class EncodeFile:
    kind: str  # "peak" or "fc"
    accession: str
    url: str
    output_type: str
    file_format: str
    assembly: str
    status: str
    md5sum: Optional[str]
    file_size: Optional[int]
    date_released: Optional[str]
    date_created: Optional[str]
    preferred_default: bool  # ⭐ on ENCODE UI
    default: bool            # often also true for ⭐


@dataclass
class ExperimentSelection:
    target: str
    exp_accession: str
    exp_url: str
    exp_title: str
    peak: Optional[EncodeFile]
    fc: Optional[EncodeFile]


def http_get_json(url: str) -> Dict[str, Any]:
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": UA},
        method="GET",
    )
    with urllib.request.urlopen(req) as resp:
        data = resp.read()
    return json.loads(data.decode("utf-8"))


def safe_name(s: str) -> str:
    s = s.strip()
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s


def classify_file(f: Dict[str, Any]) -> Optional[str]:
    """Return kind 'peak' or 'fc' if file matches our criteria."""
    if f.get("status") != "released":
        return None
    if f.get("assembly") != "GRCh38":
        return None

    fmt = (f.get("file_format") or "").lower()
    otype = (f.get("output_type") or "").lower()

    # IDR thresholded peaks (narrowPeak or bigBed)
    if "idr thresholded peaks" in otype and fmt in ("narrowpeak", "bigbed"):
        return "peak"

    # fold-change over control bigWig
    if ("fold change over control" in otype or "fold-change over control" in otype) and fmt == "bigwig":
        return "fc"

    return None


def encode_file_url(f: Dict[str, Any]) -> str:
    href = f.get("href") or f.get("@id")
    if not href:
        raise ValueError("ENCODE file missing href/@id")
    return ENCODE + href


def _parse_encode_datetime(s: Optional[str]) -> int:
    """
    Convert ENCODE datetime/date strings to a sortable integer timestamp (seconds).
    Returns 0 if missing/unparseable.
    Common formats:
      - YYYY-MM-DD
      - YYYY-MM-DDTHH:MM:SS(.ffffff)?(Z or ±HH:MM)
    """
    if not s:
        return 0
    s = s.strip()
    try:
        if len(s) == 10 and s[4] == "-" and s[7] == "-":
            dt = datetime.strptime(s, "%Y-%m-%d").replace(tzinfo=timezone.utc)
            return int(dt.timestamp())
        s2 = s.replace("Z", "+00:00")
        dt = datetime.fromisoformat(s2)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return int(dt.timestamp())
    except Exception:
        return 0


def pick_best(files: List[EncodeFile]) -> Tuple[Optional[EncodeFile], Optional[EncodeFile]]:
    """
    Choose one peak and one fc file from candidates.

    Preference:
      1) ⭐ preferred_default, then default
      2) Newer date (date_released/date_created)
      3) Peaks: narrowPeak > bigBed, then size
         FC: size
    """
    peaks = [x for x in files if x.kind == "peak"]
    fcs = [x for x in files if x.kind == "fc"]

    def star_score(x: EncodeFile) -> Tuple[int, int]:
        return (1 if x.preferred_default else 0, 1 if x.default else 0)

    def date_score(x: EncodeFile) -> int:
        return max(_parse_encode_datetime(x.date_released), _parse_encode_datetime(x.date_created))

    def peak_score(x: EncodeFile):
        fmt_score = 2 if x.file_format.lower() == "narrowpeak" else (1 if x.file_format.lower() == "bigbed" else 0)
        size = x.file_size or 0
        return (star_score(x), date_score(x), fmt_score, size)

    def fc_score(x: EncodeFile):
        size = x.file_size or 0
        return (star_score(x), date_score(x), size)

    peak = max(peaks, key=peak_score) if peaks else None
    fc = max(fcs, key=fc_score) if fcs else None
    return peak, fc


def download_file(url: str, outpath: Path, expected_md5: Optional[str] = None, overwrite: bool = False) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists() and outpath.stat().st_size > 0 and not overwrite:
        if expected_md5:
            if md5sum_file(outpath).lower() == expected_md5.lower():
                return
        else:
            return

    req = urllib.request.Request(url, headers={"User-Agent": UA}, method="GET")
    with urllib.request.urlopen(req) as resp, open(outpath, "wb") as w:
        h = hashlib.md5() if expected_md5 else None
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            w.write(chunk)
            if h:
                h.update(chunk)

    if expected_md5:
        got = md5sum_file(outpath)
        if got.lower() != expected_md5.lower():
            raise RuntimeError(f"MD5 mismatch for {outpath.name}: expected {expected_md5} got {got}")


def md5sum_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as r:
        for chunk in iter(lambda: r.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def experiment_search_url(target_label: str) -> str:
    params = [
        ("type", "Experiment"),
        ("status", "released"),
        ("assay_title", "TF ChIP-seq"),
        ("biosample_ontology.term_name", "GM12878"),
        ("target.label", target_label),
        ("limit", "all"),
        ("format", "json"),
        ("frame", "object"),
    ]
    return ENCODE + "/search/?" + urllib.parse.urlencode(params)


def fetch_experiment_full(exp_obj: Dict[str, Any]) -> Dict[str, Any]:
    exp_id = exp_obj.get("@id")
    if not exp_id:
        raise ValueError("Experiment missing @id")
    return http_get_json(ENCODE + exp_id + "?format=json")


def select_from_experiment(target: str, exp_full: Dict[str, Any]) -> ExperimentSelection:
    exp_acc = exp_full.get("accession", "")
    exp_url = ENCODE + exp_full.get("@id", f"/experiments/{exp_acc}/")
    exp_title = exp_full.get("description") or exp_full.get("title") or ""

    candidates: List[EncodeFile] = []
    for f in exp_full.get("files", []):
        kind = classify_file(f)
        if not kind:
            continue

        candidates.append(
            EncodeFile(
                kind=kind,
                accession=f.get("accession", ""),
                url=encode_file_url(f),
                output_type=f.get("output_type", ""),
                file_format=f.get("file_format", ""),
                assembly=f.get("assembly", ""),
                status=f.get("status", ""),
                md5sum=f.get("md5sum"),
                file_size=f.get("file_size"),
                date_released=f.get("date_released"),
                date_created=f.get("date_created"),
                preferred_default=bool(f.get("preferred_default", False)),
                default=bool(f.get("default", False)),
            )
        )

    peak, fc = pick_best(candidates)
    return ExperimentSelection(
        target=target,
        exp_accession=exp_acc,
        exp_url=exp_url,
        exp_title=exp_title,
        peak=peak,
        fc=fc,
    )


def write_metadata(sel: ExperimentSelection, md_path: Path) -> None:
    md_path.parent.mkdir(parents=True, exist_ok=True)

    def block(label: str, f: Optional[EncodeFile]) -> str:
        if not f:
            return f"{label}:\tNOT_FOUND\n"
        return "\n".join(
            [
                f"{label}_accession:\t{f.accession}",
                f"{label}_url:\t{f.url}",
                f"{label}_output_type:\t{f.output_type}",
                f"{label}_file_format:\t{f.file_format}",
                f"{label}_assembly:\t{f.assembly}",
                f"{label}_status:\t{f.status}",
                f"{label}_preferred_default:\t{f.preferred_default}",
                f"{label}_default:\t{f.default}",
                f"{label}_date_released:\t{f.date_released or ''}",
                f"{label}_date_created:\t{f.date_created or ''}",
                f"{label}_md5:\t{f.md5sum or ''}",
                f"{label}_file_size:\t{f.file_size if f.file_size is not None else ''}",
                "",
            ]
        )

    txt = "\n".join(
        [
            f"target:\t{sel.target}",
            f"experiment_accession:\t{sel.exp_accession}",
            f"experiment_url:\t{sel.exp_url}",
            f"experiment_title:\t{sel.exp_title}",
            "",
            block("idr_thresholded_peaks", sel.peak),
            block("fold_change_over_control", sel.fc),
        ]
    )
    md_path.write_text(txt, encoding="utf-8")


def append_summary_tsv(sel: ExperimentSelection, tsv_path: Path) -> None:
    header = "\t".join(
        [
            "target",
            "experiment_accession",
            "experiment_title",
            "experiment_url",
            "peak_file_accession",
            "peak_file_format",
            "peak_output_type",
            "peak_url",
            "peak_preferred_default",
            "peak_default",
            "peak_date_released",
            "peak_date_created",
            "peak_md5",
            "peak_size",
            "fc_file_accession",
            "fc_output_type",
            "fc_url",
            "fc_preferred_default",
            "fc_default",
            "fc_date_released",
            "fc_date_created",
            "fc_md5",
            "fc_size",
        ]
    )

    if not tsv_path.exists():
        tsv_path.write_text(header + "\n", encoding="utf-8")

    def v(x: Any) -> str:
        if x is None:
            return ""
        s = str(x)
        return s.replace("\t", " ").replace("\n", " ")

    peak = sel.peak
    fc = sel.fc

    row = "\t".join(
        [
            v(sel.target),
            v(sel.exp_accession),
            v(sel.exp_title),
            v(sel.exp_url),
            v(peak.accession if peak else ""),
            v(peak.file_format if peak else ""),
            v(peak.output_type if peak else ""),
            v(peak.url if peak else ""),
            v(peak.preferred_default if peak else ""),
            v(peak.default if peak else ""),
            v(peak.date_released if peak else ""),
            v(peak.date_created if peak else ""),
            v(peak.md5sum if peak else ""),
            v(peak.file_size if peak else ""),
            v(fc.accession if fc else ""),
            v(fc.output_type if fc else ""),
            v(fc.url if fc else ""),
            v(fc.preferred_default if fc else ""),
            v(fc.default if fc else ""),
            v(fc.date_released if fc else ""),
            v(fc.date_created if fc else ""),
            v(fc.md5sum if fc else ""),
            v(fc.file_size if fc else ""),
        ]
    )

    with open(tsv_path, "a", encoding="utf-8") as w:
        w.write(row + "\n")


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    ap.add_argument("--outdir", default="encode_gm12878_hg38_tfchip", help="Output directory")
    ap.add_argument("--targets", nargs="+", default=["CTCF", "RAD21"], help="Targets (default: CTCF RAD21)")
    ap.add_argument("--overwrite", action="store_true", help="Re-download files even if present")
    ap.add_argument("--no-md5", action="store_true", help="Skip md5 verification")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    downloads = outdir / "downloads"
    metadata = outdir / "metadata"
    outdir.mkdir(parents=True, exist_ok=True)
    downloads.mkdir(parents=True, exist_ok=True)
    metadata.mkdir(parents=True, exist_ok=True)

    summary_tsv = outdir / "encode_gm12878_CTCF_RAD21_hg38_download_summary.tsv"

    for target in args.targets:
        search_url = experiment_search_url(target)
        search = http_get_json(search_url)
        exps = search.get("@graph", [])
        if not exps:
            print(f"[WARN] No experiments found for target={target}", file=sys.stderr)
            continue

        for exp_obj in exps:
            exp_full = fetch_experiment_full(exp_obj)
            sel = select_from_experiment(target, exp_full)

            md_path = metadata / f"{safe_name(sel.target)}_{safe_name(sel.exp_accession)}.metadata.txt"
            write_metadata(sel, md_path)
            append_summary_tsv(sel, summary_tsv)

            if sel.peak:
                peak_ext = sel.peak.file_format
                peak_out = downloads / sel.target / sel.exp_accession / f"idr_peaks.{sel.peak.accession}.{peak_ext}"
                download_file(
                    sel.peak.url,
                    peak_out,
                    None if args.no_md5 else sel.peak.md5sum,
                    overwrite=args.overwrite,
                )

            if sel.fc:
                fc_ext = sel.fc.file_format
                fc_out = downloads / sel.target / sel.exp_accession / f"fc_over_control.{sel.fc.accession}.{fc_ext}"
                download_file(
                    sel.fc.url,
                    fc_out,
                    None if args.no_md5 else sel.fc.md5sum,
                    overwrite=args.overwrite,
                )

            print(
                f"{sel.target} {sel.exp_accession}: "
                f"peak={'OK' if sel.peak else 'MISSING'} fc={'OK' if sel.fc else 'MISSING'} "
                f"(peak_star={sel.peak.preferred_default if sel.peak else ''} "
                f"fc_star={sel.fc.preferred_default if sel.fc else ''})"
            )

    print(f"\nSummary TSV: {summary_tsv}")
    print(f"Metadata dir: {metadata}")
    print(f"Downloads dir: {downloads}")


if __name__ == "__main__":
    main()
