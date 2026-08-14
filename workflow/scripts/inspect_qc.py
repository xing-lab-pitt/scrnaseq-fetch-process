#!/usr/bin/env python3
"""Report STARsolo QC metrics for every sample in a workdir, as one compact table.

This is the read-only inspection counterpart to `qc_gate.py`. The gate runs inside
the DAG and stops the workflow; this answers "how did the run actually do?" at any
point, including part-way through, without launching anything.

It imports `qc_gate.evaluate`, so the numbers here and the numbers the pipeline
enforces are the same code path and cannot drift apart.

    python inspect_qc.py <workdir> [--config config/config.yaml]
    python inspect_qc.py <workdir> --feature GeneFull --min-valid-barcodes 0.7

Thresholds resolve in this order: explicit --min-* flag, then the `qc:` block of
--config, then the built-in defaults below. --feature likewise falls back to the
config's `feature:` key, since the Summary.csv lives under Solo.out/<feature>/.

Exit status is 0 when every discovered sample passes, 1 when any fails or is
missing — so it doubles as a check in a shell pipeline.

Intended to be called through context-mode so the table is indexed, not dumped:

    ctx_execute(language="shell",
                code="python inspect_qc.py $W --config config/config.yaml",
                intent="which GSE219015 samples passed the QC gate")
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_gate import evaluate  # noqa: E402  (needs the path insert above)

DEFAULTS = {
    "min_valid_barcodes": 0.70,
    "min_reads_mapped_gene": 0.50,
    "min_estimated_cells": 100,
}
COLS = ["valid_barcodes", "reads_mapped_gene", "estimated_cells", "saturation"]


def load_config(path):
    """Pull the `qc:` thresholds and `feature:` out of the pipeline config, if given."""
    if not path:
        return {}, None
    try:
        import yaml
    except ImportError:
        sys.exit("--config needs PyYAML; omit it and pass --min-* flags instead")
    cfg = yaml.safe_load(Path(path).read_text()) or {}
    return (cfg.get("qc") or {}), cfg.get("feature")


def discover(workdir: Path, feature: str):
    """Yield (sample, summary_path_or_None) for every aligned sample, sorted."""
    star = workdir / "star"
    if not star.is_dir():
        sys.exit(f"no star/ directory under {workdir} — has anything aligned yet?")
    for sample_dir in sorted(p for p in star.iterdir() if p.is_dir()):
        summary = sample_dir / "Solo.out" / feature / "Summary.csv"
        yield sample_dir.name, (summary if summary.is_file() else None)


def fmt(val):
    if val is None:
        return "-"
    return f"{val:.4g}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("workdir", help="pipeline workdir (contains star/, h5ad/, qc/)")
    ap.add_argument("--config", help="pipeline config.yaml, for qc: thresholds and feature:")
    ap.add_argument("--feature", help="Gene or GeneFull (default: from --config, else Gene)")
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", type=float, default=None,
                        help=f"default {v}")
    a = ap.parse_args()

    cfg_qc, cfg_feature = load_config(a.config)
    feature = a.feature or cfg_feature or "Gene"
    thresholds = {
        k: (getattr(a, k) if getattr(a, k) is not None else cfg_qc.get(k, v))
        for k, v in DEFAULTS.items()
    }

    workdir = Path(a.workdir)
    rows, failed = [], 0
    for sample, summary in discover(workdir, feature):
        if summary is None:
            rows.append((sample, "MISSING", {}, f"no Solo.out/{feature}/Summary.csv"))
            failed += 1
            continue
        metrics, status, reasons = evaluate(summary, thresholds)
        rows.append((sample, status, metrics, reasons))
        if status != "PASS":
            failed += 1

    if not rows:
        print(f"no samples found under {workdir}/star/")
        return 1

    width = max(len(r[0]) for r in rows)
    thr_desc = ", ".join(f"{k.replace('min_', '')}>={v}" for k, v in thresholds.items())
    print(f"workdir: {workdir}")
    print(f"feature: {feature}    thresholds: {thr_desc}")
    print()
    cw = {c: max(10, len(c)) for c in COLS}          # header may be wider than the value
    header = f"{'sample':<{width}}  {'status':<7}  " + "  ".join(f"{c:>{cw[c]}}" for c in COLS)
    print(header)
    print("-" * len(header))
    for sample, status, metrics, reasons in rows:
        cells = "  ".join(f"{fmt(metrics.get(c)):>{cw[c]}}" for c in COLS)
        print(f"{sample:<{width}}  {status:<7}  {cells}")
        if reasons:
            print(f"{'':<{width}}  └─ {reasons}")
    print()
    print(f"{len(rows) - failed}/{len(rows)} passed")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
