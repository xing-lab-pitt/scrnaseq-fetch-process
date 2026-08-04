#!/usr/bin/env python3
"""Deterministic QC gate over STARsolo Gene/Summary.csv files.

Applies FIXED thresholds (passed in from config) to each sample and writes a
pass/fail table + a list of passing samples. It does NOT choose thresholds
itself, so the decision is reproducible and documentable. Always exits 0 (it
produces the report); enforcement (stop vs. warn) is a separate step so the
report is never deleted on failure.
"""
import argparse
import csv
import sys
from pathlib import Path

# Metric -> substrings matched (case-insensitive) against Summary.csv keys.
# STARsolo key names vary slightly by version, so we match loosely.
METRIC_KEYS = {
    "valid_barcodes":    ["reads with valid barcodes"],
    "reads_mapped_gene": ["reads mapped to gene: unique+multiple",
                          "reads mapped to genefull: unique+multiple",
                          "reads mapped to gene: unique gene"],
    "estimated_cells":   ["estimated number of cells"],
    "saturation":        ["sequencing saturation"],
}


def parse_summary(path):
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh):
            if len(row) >= 2:
                out[row[0].strip().lower()] = row[1].strip()
    return out


def get_metric(summary, names):
    for key, val in summary.items():
        if any(n in key for n in names):
            try:
                return float(val)
            except ValueError:
                return None
    return None


def evaluate(summary_path, thresholds):
    s = parse_summary(summary_path)
    metrics = {m: get_metric(s, keys) for m, keys in METRIC_KEYS.items()}
    reasons = []

    def check_min(name, thr):
        if thr is None:
            return
        val = metrics[name]
        if val is None:
            reasons.append(f"{name}:missing")
        elif val < thr:
            reasons.append(f"{name}={val:.4g}<{thr}")

    check_min("valid_barcodes", thresholds.get("min_valid_barcodes"))
    check_min("reads_mapped_gene", thresholds.get("min_reads_mapped_gene"))
    check_min("estimated_cells", thresholds.get("min_estimated_cells"))
    status = "PASS" if not reasons else "FAIL"
    return metrics, status, ";".join(reasons)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("entries", nargs="+", help="sample=<path/to/Gene/Summary.csv>")
    ap.add_argument("--report", required=True)
    ap.add_argument("--passlist", required=True)
    ap.add_argument("--min-valid-barcodes", type=float)
    ap.add_argument("--min-reads-mapped-gene", type=float)
    ap.add_argument("--min-estimated-cells", type=float)
    args = ap.parse_args()

    thresholds = {
        "min_valid_barcodes": args.min_valid_barcodes,
        "min_reads_mapped_gene": args.min_reads_mapped_gene,
        "min_estimated_cells": args.min_estimated_cells,
    }

    rows, passing = [], []
    for entry in args.entries:
        sample, _, path = entry.partition("=")
        metrics, status, reason = evaluate(path, thresholds)
        rows.append({
            "sample": sample,
            "valid_barcodes": metrics["valid_barcodes"],
            "reads_mapped_gene": metrics["reads_mapped_gene"],
            "estimated_cells": metrics["estimated_cells"],
            "saturation": metrics["saturation"],
            "status": status,
            "reasons": reason,
        })
        if status == "PASS":
            passing.append(sample)

    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    Path(args.passlist).write_text("\n".join(passing) + ("\n" if passing else ""))

    print(f"QC gate: {len(passing)}/{len(rows)} passed. Report: {args.report}")
    for r in rows:
        if r["status"] == "FAIL":
            print(f"  FAIL {r['sample']}: {r['reasons']}", file=sys.stderr)
    sys.exit(0)  # report always succeeds; enforcement is a separate rule


if __name__ == "__main__":
    main()
