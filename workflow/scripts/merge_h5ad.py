#!/usr/bin/env python3
"""Concatenate per-sample .h5ad files into one, preserving spliced/unspliced layers.

QC-aware: pass paired --qc {sample}.qc.json files and only samples with
status == PASS are included by default. Failing samples are left on disk
(nothing is deleted) but excluded from the merged object and recorded in
adata.uns['qc_excluded']. Use --include-failed to merge everything anyway.
"""
import argparse
import json
from pathlib import Path

import anndata as ad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--h5ad", nargs="+", required=True, help="per-sample .h5ad files")
    ap.add_argument("--qc", nargs="*", default=[], help="per-sample .qc.json files")
    ap.add_argument("--include-failed", action="store_true",
                    help="merge QC-failing samples too (still recorded in uns)")
    args = ap.parse_args()

    status = {}
    for f in args.qc:
        v = json.load(open(f))
        status[v["sample"]] = v["status"]

    kept, excluded = [], []
    for f in args.h5ad:
        sample = Path(f).stem
        st = status.get(sample, "UNKNOWN")
        if st == "FAIL" and not args.include_failed:
            excluded.append(sample)
            continue
        a = ad.read_h5ad(f)
        a.obs["qc_status"] = st
        kept.append(a)

    if not kept:
        raise SystemExit("No samples passed QC; nothing to merge "
                         "(re-run with --include-failed to override).")

    combined = ad.concat(kept, join="outer", merge="same")
    combined.uns["qc_excluded"] = excluded
    combined.write_h5ad(args.output)
    print(f"Merged {len(kept)} samples -> {args.output}: {combined.n_obs} cells, "
          f"layers={list(combined.layers)}; excluded (QC fail): {excluded or 'none'}")


if __name__ == "__main__":
    main()
