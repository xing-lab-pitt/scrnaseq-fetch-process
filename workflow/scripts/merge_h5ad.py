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
    ap.add_argument("--var-names", choices=["gene_ids", "symbols"], default="gene_ids",
                    help="what to use as var_names in the merged object. "
                         "'gene_ids' (default): Ensembl IDs, unique by construction, "
                         "so cross-dataset merges are collision-free; raw symbols are "
                         "kept in var['gene_symbols']. 'symbols': legacy behaviour "
                         "(scanpy-style), de-duplicated with -1/-2 suffixes.")
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
        # CHECKPOINT: Ensembl/GENCODE GTFs map several gene IDs to the SAME gene
        # symbol, so var_names (symbols) are not guaranteed unique. ad.concat(
        # merge="same") reindexes on the var axis and requires a unique index, so
        # a duplicate here would crash the merge with pandas InvalidIndexError.
        # This is intrinsic to the annotation — switching reference GTFs does NOT
        # avoid it; only Ensembl gene IDs are unique by construction.
        #
        # Default (--var-names gene_ids): index on the unique gene IDs so the
        # merge is collision-free, and keep the RAW symbols in var['gene_symbols']
        # (uniqueness is required on the index, not on metadata, so we do NOT
        # de-duplicate the symbol column). --var-names symbols restores the legacy
        # scanpy-style behaviour, de-duplicating symbols with -1/-2 suffixes.
        symbols = a.var_names  # per-sample .h5ad indexes on symbols (see starsolo_to_h5ad)
        dups = symbols[symbols.duplicated()].unique()
        if len(dups):
            preview = ", ".join(map(str, dups[:10]))
            more = f", … (+{len(dups) - 10} more)" if len(dups) > 10 else ""
            print(f"  [WARN] {sample}: {len(dups)} gene symbol(s) are non-unique "
                  f"({preview}{more}).")
        if args.var_names == "gene_ids":
            if "gene_ids" not in a.var:
                raise SystemExit(
                    f"{sample}: --var-names gene_ids requested but var['gene_ids'] "
                    f"is missing; re-run starsolo_to_h5ad or use --var-names symbols.")
            # Keep the RAW symbols as a column (uniqueness is only required on the
            # index, not on metadata — de-duplicating the symbols here would just
            # corrupt the real gene symbol for no benefit), then key the index on
            # the unique gene IDs so the concat is collision-free.
            a.var["gene_symbols"] = symbols.astype(str).values
            a.var_names = a.var["gene_ids"].astype(str)
            a.var_names_make_unique()  # gene_ids are already unique; cheap safety net
            print(f"  [{sample}] var_names -> gene_ids "
                  f"(raw symbols preserved in var['gene_symbols']).")
        elif len(dups):
            print(f"  [{sample}] var_names -> symbols; de-duplicating with -1/-2 suffixes.")
            a.var_names_make_unique()
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
