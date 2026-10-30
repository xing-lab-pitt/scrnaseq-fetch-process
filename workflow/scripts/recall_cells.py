#!/usr/bin/env python3
"""Re-call cells on existing STARsolo raw matrices, without re-aligning anything.

Wraps `STAR --runMode soloCellFiltering`, which re-applies a --soloCellFilter to a raw
count matrix that has already been produced by the `starsolo` rule. No FASTQ is read
and nothing is re-mapped, so this runs in minutes per sample instead of hours.

Generalized from the reference implementation that first proved this out on GSE310923
(analysis/scripts/recall_cells.sh in the aged_mouse_blood project), which found two
things worth carrying over to any future study:

1. STAR's soloCellFiltering mode hard-requires a file literally named `matrix.mtx` in
   the input dir. The Gene/GeneFull raw dir has this; the Velocyto raw dir instead
   holds spliced.mtx / unspliced.mtx / ambiguous.mtx, so soloCellFiltering cannot be
   run on it directly. Pass --recall-velocyto to slice those three matrices by the
   barcode list the Gene/GeneFull pass just called, instead of re-filtering them.
2. Every sample must get the IDENTICAL --soloCellFilter string. This script takes one
   --filter value for the whole run (not a per-sample option), so that is true by
   construction; do not call this script per-sample with different --filter values
   for the same study.

Usage:
    python workflow/scripts/recall_cells.py \\
        --raw-root  <workdir>/star \\
        --out-root  <workdir>/star_relaxed \\
        --samples sample1 sample2 ... \\
        --filter "EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000" \\
        --recall-velocyto

Not a Snakemake rule and not in `rule all`: this is a deliberate, on-demand re-call of
already-aligned data, run once a study's default cell filter has been reconsidered —
not something every pipeline run should redo automatically.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from pathlib import Path

try:
    from scipy.io import mmread, mmwrite
except ImportError:  # pragma: no cover - only needed for --recall-velocyto
    mmread = mmwrite = None


def _open_maybe_gz(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def _find(dirpath: Path, stem: str) -> Path:
    """Return `<dirpath>/<stem>` or `<dirpath>/<stem>.gz`, whichever exists."""
    plain, gz = dirpath / stem, dirpath / f"{stem}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"neither {plain} nor {gz} exists")


def recall_gene(sample: str, raw_dir: Path, out_dir: Path, run_dir: Path,
                 filter_args: str, star_bin: str, threads: int) -> Path:
    """Run STAR --runMode soloCellFiltering on one sample's Gene/GeneFull raw dir.

    Returns the path to the new barcodes.tsv (needed by recall_velocyto below).
    """
    matrix = _find(raw_dir, "matrix.mtx")
    if matrix.suffix == ".gz":
        raise SystemExit(
            f"{sample}: {matrix} is gzipped; STAR's soloCellFiltering mode reads plain "
            "matrix.mtx only. Decompress the raw dir before re-calling."
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)
    # Run with cwd=run_dir (a private, per-sample scratch dir) so STAR's own Log.out
    # never collides between samples running one after another.
    subprocess.run(
        [star_bin, "--runMode", "soloCellFiltering", str(raw_dir) + "/", str(out_dir) + "/",
         "--soloCellFilter", *filter_args.split(), "--runThreadN", str(threads)],
        cwd=run_dir, check=True,
    )
    barcodes = out_dir / "barcodes.tsv"
    if not barcodes.exists():
        raise SystemExit(f"{sample}: soloCellFiltering did not write {barcodes}")
    return barcodes


def recall_velocyto(sample: str, velocyto_raw_dir: Path, new_barcodes: Path,
                     out_dir: Path) -> None:
    """Slice the Velocyto raw spliced/unspliced/ambiguous matrices by a new barcode list.

    STAR cannot re-call Velocyto directly (see module docstring), so this reproduces the
    working fallback: read the ORIGINAL barcode order from Velocyto/raw, look up each
    barcode the Gene/GeneFull pass just called, and keep only those columns, in that order.
    """
    if mmread is None:
        raise SystemExit("scipy is required for --recall-velocyto (scipy.io.mmread/mmwrite)")

    with _open_maybe_gz(_find(velocyto_raw_dir, "barcodes.tsv")) as fh:
        raw_barcodes = [line.strip() for line in fh]
    raw_index = {bc: i for i, bc in enumerate(raw_barcodes)}

    with _open_maybe_gz(new_barcodes) as fh:
        wanted = [line.strip() for line in fh]
    missing = [bc for bc in wanted if bc not in raw_index]
    if missing:
        raise SystemExit(
            f"{sample}: {len(missing)} of {len(wanted)} newly-called barcodes are not "
            f"present in {velocyto_raw_dir}/barcodes.tsv (e.g. {missing[0]!r}). "
            "The Gene/GeneFull raw dir and the Velocyto raw dir must come from the same "
            "starsolo run."
        )
    cols = [raw_index[bc] for bc in wanted]

    out_dir.mkdir(parents=True, exist_ok=True)
    for name in ("spliced", "unspliced", "ambiguous"):
        mat = mmread(str(_find(velocyto_raw_dir, f"{name}.mtx"))).tocsc()
        sliced = mat[:, cols]
        mmwrite(str(out_dir / f"{name}.mtx"), sliced)

    with open(out_dir / "barcodes.tsv", "w") as fh:
        fh.write("\n".join(wanted) + "\n")
    features_src = _find(velocyto_raw_dir, "features.tsv")
    with _open_maybe_gz(features_src) as src, open(out_dir / "features.tsv", "w") as dst:
        dst.write(src.read())


def samples_from_tsv(path: str) -> list[str]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    seen = dict.fromkeys(r["sample"] for r in rows)
    return list(seen)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-root", required=True,
                     help="dir containing <sample>/Solo.out/<feature>/raw for every sample "
                          "(e.g. the pipeline's <workdir>/star)")
    ap.add_argument("--out-root", required=True,
                     help="dir to write <sample>/Solo.out/<feature>/filtered into; "
                          "the raw-root tree is opened READ-ONLY")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--samples", nargs="+", help="sample names, space-separated")
    group.add_argument("--samples-tsv", help="samples.tsv to read the sample list from")
    ap.add_argument("--filter", required=True,
                     help="the --soloCellFilter value, e.g. "
                          '"EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000". '
                          "Applied IDENTICALLY to every sample in this run.")
    ap.add_argument("--feature", default="Gene", choices=("Gene", "GeneFull"))
    ap.add_argument("--recall-velocyto", action="store_true",
                     help="also slice Velocyto/raw spliced/unspliced/ambiguous by the "
                          "barcode list this run just called for --feature")
    ap.add_argument("--star-bin", default="STAR")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    samples = args.samples if args.samples else samples_from_tsv(args.samples_tsv)
    raw_root, out_root = Path(args.raw_root), Path(args.out_root)

    for sample in samples:
        raw_dir = raw_root / sample / "Solo.out" / args.feature / "raw"
        out_dir = out_root / sample / "Solo.out" / args.feature / "filtered"
        run_dir = out_root / sample / "_recall_run"
        if not raw_dir.exists():
            raise SystemExit(f"{sample}: missing {raw_dir} — check --raw-root and --feature")

        barcodes = recall_gene(sample, raw_dir, out_dir, run_dir,
                                args.filter, args.star_bin, args.threads)
        n_cells = sum(1 for _ in open(barcodes))
        print(f"{sample}: {n_cells} cells called ({args.feature})")

        if args.recall_velocyto:
            velo_raw = raw_root / sample / "Solo.out" / "Velocyto" / "raw"
            velo_out = out_root / sample / "Solo.out" / "Velocyto" / "filtered"
            recall_velocyto(sample, velo_raw, barcodes, velo_out)
            print(f"{sample}: Velocyto layers sliced to the same {n_cells} barcodes")

    print(f"\n{len(samples)} sample(s) re-called with --soloCellFilter {args.filter!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
