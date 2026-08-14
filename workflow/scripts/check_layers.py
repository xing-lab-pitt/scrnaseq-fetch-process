#!/usr/bin/env python3
"""Check whether per-sample .h5ad files carry their RNA-velocity layers, and if
not, say which of the two causes it is.

An .h5ad with no spliced/unspliced/ambiguous layers looks identical whether
STARsolo never counted them or the converter failed to attach them — but the
fixes are in different files, so this reads BOTH ends and names the cause:

  OK              all expected layers present and non-empty
  ATTACHED_EMPTY  layers present but every one has nnz=0
  NOT_ATTACHED    Solo.out/Velocyto has counts, h5ad does not
                  -> converter bug: check which subdir starsolo_to_h5ad.py reads.
                     STARsolo writes Velocyto/raw/ only, never Velocyto/filtered/
  UPSTREAM_EMPTY  Solo.out/Velocyto/*.mtx exist but are all nnz=0
                  -> STARsolo counted nothing: `Gene` must be in --soloFeatures.
                     The Velocyto counter consumes the Gene pass, so
                     `--soloFeatures GeneFull Velocyto` silently yields zeros;
                     it needs `--soloFeatures Gene GeneFull Velocyto`
  NO_VELOCYTO     no Solo.out/Velocyto at all -> Velocyto not requested
  MISSING         no .h5ad for this sample yet

Uses h5py to read only the group structure, so matrices are never loaded into
memory — this stays fast and flat regardless of how big the h5ad is.

    python check_layers.py <workdir>
    python check_layers.py <workdir> --layers spliced,unspliced

Exit status 0 when every sample is OK, 1 otherwise, so it works as a gate before
deleting FASTQs or starting downstream velocity work.

Intended to be called through context-mode so the table is indexed, not dumped:

    ctx_execute(language="shell",
                code="python check_layers.py $W",
                intent="which samples are missing velocity layers and why")
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

EXPECTED = ["spliced", "unspliced", "ambiguous"]

HINTS = {
    "ATTACHED_EMPTY": "layers attached but all zero — same root cause as UPSTREAM_EMPTY.",
    "NOT_ATTACHED": "Solo.out/Velocyto has counts the h5ad lacks. Converter bug: "
                    "starsolo_to_h5ad.py reads Velocyto/<sub> where sub follows the "
                    "Gene matrix ('filtered'), but STARsolo only ever writes "
                    "Velocyto/raw/. Point it at raw/ and let reindex() realign.",
    "UPSTREAM_EMPTY": "STARsolo counted nothing. Its Velocyto counter consumes the "
                      "Gene pass, so '--soloFeatures GeneFull Velocyto' yields all "
                      "zeros silently. Use '--soloFeatures Gene GeneFull Velocyto'.",
    "NO_VELOCYTO": "no Solo.out/Velocyto — Velocyto was not requested in --soloFeatures.",
    "MISSING": "no .h5ad yet — the sample has not reached to_h5ad.",
}


def h5ad_layers(path):
    """Return (n_obs, n_var, {layer: nnz}) reading only structure, never matrices."""
    import h5py

    with h5py.File(path, "r") as f:
        def _n(name):
            g = f[name]
            idx = g.attrs.get("_index", "_index")
            if isinstance(idx, bytes):
                idx = idx.decode()
            return g[idx].shape[0]

        n_obs, n_var = _n("obs"), _n("var")
        layers = {}
        if "layers" in f:
            for name, node in f["layers"].items():
                # sparse layers are a group with a `data` dataset; dense is a dataset
                if hasattr(node, "keys") and "data" in node:
                    layers[name] = int(node["data"].shape[0])
                else:
                    layers[name] = int(getattr(node, "size", 0))
        return n_obs, n_var, layers


def velocyto_nnz(solo_dir: Path):
    """Return {layer: nnz} from the Velocyto .mtx headers, or None if absent.

    Reads the third line of each MatrixMarket file (rows cols nnz) — no parsing
    of the body, so this is O(1) per matrix no matter how large it is.
    """
    velo = solo_dir / "Velocyto"
    if not velo.is_dir():
        return None
    out = {}
    for sub in ("raw", "filtered"):
        d = velo / sub
        if not d.is_dir():
            continue
        for layer in EXPECTED:
            mtx = d / f"{layer}.mtx"
            if not mtx.is_file():
                continue
            with mtx.open() as fh:
                for i, line in enumerate(fh):
                    if i == 2:                       # rows cols nnz
                        try:
                            out[layer] = int(line.split()[2])
                        except (IndexError, ValueError):
                            out[layer] = -1
                        break
        if out:
            break                                    # prefer raw/, which is what STAR writes
    return out


def classify(h5ad: Path, solo: Path, expected):
    if not h5ad.is_file():
        return "MISSING", "", {}
    n_obs, n_var, layers = h5ad_layers(h5ad)
    shape = f"{n_obs}x{n_var}"
    have = {k: v for k, v in layers.items() if k in expected}

    if have and all(k in have for k in expected):
        return ("OK" if any(v > 0 for v in have.values()) else "ATTACHED_EMPTY"), shape, have

    up = velocyto_nnz(solo)
    if up is None:
        return "NO_VELOCYTO", shape, have
    if up and any(v > 0 for v in up.values()):
        return "NOT_ATTACHED", shape, have
    return "UPSTREAM_EMPTY", shape, have


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("workdir", help="pipeline workdir (contains star/, h5ad/)")
    ap.add_argument("--layers", default=",".join(EXPECTED),
                    help=f"comma-separated layers to require (default: {','.join(EXPECTED)})")
    a = ap.parse_args()

    expected = [s.strip() for s in a.layers.split(",") if s.strip()]
    workdir = Path(a.workdir)
    h5dir, stardir = workdir / "h5ad", workdir / "star"
    if not stardir.is_dir():
        sys.exit(f"no star/ directory under {workdir}")

    samples = sorted(p.name for p in stardir.iterdir() if p.is_dir())
    if not samples:
        print(f"no samples under {stardir}")
        return 1

    rows = []
    for s in samples:
        status, shape, have = classify(h5dir / f"{s}.h5ad", stardir / s / "Solo.out", expected)
        rows.append((s, status, shape, have))

    width = max(len(r[0]) for r in rows)
    print(f"workdir: {workdir}")
    print(f"expecting layers: {', '.join(expected)}")
    print()
    print(f"{'sample':<{width}}  {'status':<15}  {'shape':>14}  layer nnz")
    print("-" * (width + 50))
    bad = 0
    for s, status, shape, have in rows:
        detail = ", ".join(f"{k}={v}" for k, v in sorted(have.items())) or "-"
        print(f"{s:<{width}}  {status:<15}  {shape:>14}  {detail}")
        if status != "OK":
            bad += 1

    print()
    print(f"{len(rows) - bad}/{len(rows)} OK")
    if bad:
        print()
        for status in sorted({r[1] for r in rows if r[1] != "OK"}):
            print(f"{status}: {HINTS.get(status, 'no hint available')}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
