#!/usr/bin/env python3
"""Convert STARsolo output (Gene + Velocyto features) into a single .h5ad.

Layout produced:
    adata.X                 = Gene counts        (cells x genes, raw counts)
    adata.layers['spliced']   } aligned to X's
    adata.layers['unspliced'] } cells and genes, from the
    adata.layers['ambiguous'] } Velocyto feature

var_names are gene symbols (matching scanpy.read_10x_mtx), gene IDs kept in
var['gene_ids']. The result loads directly via run_pipeline.load_data().

The matrix logic lives in build_matrices() so it can be tested without anndata.
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp


def _find(directory: Path, names):
    """Return first existing file among candidates (plain or .gz)."""
    for name in names:
        for cand in (directory / name, directory / (name + ".gz")):
            if cand.exists():
                return cand
    raise FileNotFoundError(f"None of {names} found in {directory}")


def _read_lines(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


def _read_mtx(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            return scipy.io.mmread(fh)
    return scipy.io.mmread(str(path))


def load_matrix(mtx_dir, mtx_name):
    """Load one STARsolo matrix dir -> (barcodes, gene_ids, gene_names, M).

    STARsolo writes features x barcodes; we transpose to cells x genes.
    """
    mtx_dir = Path(mtx_dir)
    barcodes = _read_lines(_find(mtx_dir, ["barcodes.tsv"]))
    feats = [ln.split("\t") for ln in _read_lines(_find(mtx_dir, ["features.tsv", "genes.tsv"]))]
    gene_ids = [f[0] for f in feats]
    gene_names = [f[1] if len(f) > 1 else f[0] for f in feats]
    M = sp.csr_matrix(_read_mtx(_find(mtx_dir, [mtx_name])))  # features x barcodes
    M = M.T.tocsr()                                           # cells x genes
    return barcodes, gene_ids, gene_names, M


def reindex(M, src_rows, src_cols, tgt_rows, tgt_cols):
    """Reindex sparse M (src_rows x src_cols) onto (tgt_rows x tgt_cols).

    Missing rows/cols become zeros. Robust to different ordering or barcode
    sets between the Gene and Velocyto features.
    """
    M = sp.csr_matrix(M)
    n_r, n_c = M.shape
    # Pad one zero row and one zero col; missing keys map to the pad index.
    M = sp.vstack([M, sp.csr_matrix((1, n_c))]).tocsr()
    M = sp.hstack([M, sp.csr_matrix((n_r + 1, 1))]).tocsr()
    r_map = {b: i for i, b in enumerate(src_rows)}
    c_map = {g: i for i, g in enumerate(src_cols)}
    row_idx = np.fromiter((r_map.get(b, n_r) for b in tgt_rows), dtype=int, count=len(tgt_rows))
    col_idx = np.fromiter((c_map.get(g, n_c) for g in tgt_cols), dtype=int, count=len(tgt_cols))
    return sp.csr_matrix(M[row_idx][:, col_idx])


def build_matrices(solo_dir, use_filtered=True):
    """Assemble Gene counts + aligned Velocyto layers. No anndata dependency."""
    solo = Path(solo_dir)
    sub = "filtered" if use_filtered else "raw"
    barcodes, gene_ids, gene_names, X = load_matrix(solo / "Gene" / sub, "matrix.mtx")

    layers = {}
    velo = solo / "Velocyto" / sub
    if velo.exists():
        for layer, fname in [("spliced", "spliced.mtx"),
                             ("unspliced", "unspliced.mtx"),
                             ("ambiguous", "ambiguous.mtx")]:
            vbc, vgid, _, VM = load_matrix(velo, fname)
            layers[layer] = reindex(VM, vbc, vgid, barcodes, gene_ids)
    return barcodes, gene_ids, gene_names, X, layers


def _make_unique(names):
    seen, out = {}, []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}-{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("solo_dir", help="Path to '<prefix>Solo.out' directory")
    ap.add_argument("-o", "--output", required=True, help="Output .h5ad path")
    ap.add_argument("-s", "--sample", required=True, help="Sample name -> obs['sample']")
    ap.add_argument("--raw", action="store_true", help="Use raw/ instead of filtered/")
    args = ap.parse_args()

    import anndata as ad          # imported here so tests can skip it
    import pandas as pd

    barcodes, gene_ids, gene_names, X, layers = build_matrices(
        args.solo_dir, use_filtered=not args.raw)

    var = pd.DataFrame({"gene_ids": gene_ids, "feature_types": "Gene Expression"},
                       index=_make_unique(gene_names))
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=var)
    for name, mat in layers.items():
        adata.layers[name] = mat
    adata.obs["sample"] = args.sample
    adata.obs_names = [f"{args.sample}_{bc}" for bc in adata.obs_names]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f"Wrote {args.output}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"layers={list(adata.layers)}")


if __name__ == "__main__":
    main()
