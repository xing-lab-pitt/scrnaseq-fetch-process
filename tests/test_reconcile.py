#!/usr/bin/env python3
"""Build a synthetic workdir hitting every reconciler category and assert.

Covers: done, missing, corrupt, read_qc_fail, qc_fail, no_layers — plus the
`action` each carries (rerun vs flag), and the append-only success ledger
(only DONE samples recorded, idempotent on re-run).
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Repo root = parent of this tests/ dir, so the suite runs wherever the repo lives.
PIPE = Path(__file__).resolve().parents[1]
RECON = PIPE / "workflow" / "scripts" / "reconcile.py"


def make_h5ad(path, with_layers=True):
    X = sp.csr_matrix(np.array([[1, 0], [0, 2]], dtype="float32"))
    a = ad.AnnData(X=X,
                   obs=pd.DataFrame(index=["c1", "c2"]),
                   var=pd.DataFrame(index=["g1", "g2"]))
    if with_layers:
        for l in ("spliced", "unspliced", "ambiguous"):
            a.layers[l] = X.copy()
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


def main():
    tmp = Path(tempfile.mkdtemp(prefix="recon_fx_"))
    wd = tmp / "wd"
    (wd / "h5ad").mkdir(parents=True)
    (wd / "qc").mkdir(parents=True)

    # One sample per expected category.
    samples = ["s_done", "s_missing", "s_corrupt", "s_readqc",
               "s_qcfail", "s_nolayers"]
    st = tmp / "samples.tsv"
    with open(st, "w") as fh:
        fh.write("sample\tsrr\tchemistry\n")
        for s in samples:
            fh.write(f"{s}\tSRR000\tv3\n")

    # done: h5ad + layers + in pass list + read-QC pass
    make_h5ad(wd / "h5ad" / "s_done.h5ad", with_layers=True)
    # qc_fail: h5ad present, NOT in pass list
    make_h5ad(wd / "h5ad" / "s_qcfail.h5ad", with_layers=True)
    # no_layers: h5ad present, in pass list, but no layers
    make_h5ad(wd / "h5ad" / "s_nolayers.h5ad", with_layers=False)
    # read_qc_fail: fully good on the alignment side (h5ad + layers + in pass list)
    # so we prove the read-QC FAIL wins over "done".
    make_h5ad(wd / "h5ad" / "s_readqc.h5ad", with_layers=True)
    # corrupt: h5ad present + in pass list + read-QC pass, but unreadable bytes.
    (wd / "h5ad" / "s_corrupt.h5ad").write_bytes(b"not-an-hdf5-file")
    # missing: no h5ad at all

    # QC gate outputs: pass list + gate table (with metrics for the ledger).
    (wd / "qc" / "qc_pass.txt").write_text(
        "s_done\ns_nolayers\ns_readqc\ns_corrupt\n")
    with open(wd / "qc" / "qc_gate.tsv", "w") as fh:
        fh.write("sample\tvalid_barcodes\treads_mapped_gene\testimated_cells"
                 "\tsaturation\tstatus\treasons\n")
        fh.write("s_done\t0.92\t0.64\t6201\t0.55\tPASS\t\n")
        fh.write("s_nolayers\t0.90\t0.61\t5000\t0.50\tPASS\t\n")
        fh.write("s_readqc\t0.88\t0.60\t4000\t0.50\tPASS\t\n")
        fh.write("s_corrupt\t0.90\t0.60\t5000\t0.50\tPASS\t\n")
        fh.write("s_qcfail\t0.10\t0.03\t50\t0.30\tFAIL\treads_mapped_gene=0.03<0.5\n")
        # s_missing has no Summary.csv -> not in gate table

    # Read-QC output: s_readqc FAILs; everyone else passes.
    with open(wd / "qc" / "read_qc.tsv", "w") as fh:
        fh.write("sample\tsrr\tn_reads\tcdna_len\tcdna_mean_q\tcdna_pct_adapter"
                 "\tcdna_pct_dup\tbc_len\tbc_mean_q\tflags\tstatus\treasons\n")
        for s in samples:
            if s == "s_readqc":
                fh.write(f"{s}\tSRR000\t1200000\t91\t18.4\t0.35\t95\t28\t30.1\t\t"
                         "FAIL\tlow_cdna_quality=18.4<28\n")
            else:
                fh.write(f"{s}\tSRR000\t30000000\t91\t35.5\t0.02\t65\t28\t36.2\t\t"
                         "PASS\t\n")

    report = tmp / "report.tsv"
    js = tmp / "out.json"
    ledger = tmp / "successful_samples.tsv"
    cmd = [sys.executable, str(RECON), "--samples", str(st), "--workdir", str(wd),
           "--accession", "GSE_TEST", "--feature", "Gene",
           "--report", str(report), "--json", str(js), "--ledger", str(ledger)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout)
    print(proc.stderr, file=sys.stderr)

    result = json.loads(js.read_text())
    got = {r["sample"]: r["category"] for r in result["rows"]}
    got_action = {r["sample"]: r["action"] for r in result["rows"]}
    expected = {
        "s_done": ("done", "none"),
        "s_missing": ("missing", "rerun"),
        "s_corrupt": ("corrupt", "rerun"),
        "s_readqc": ("read_qc_fail", "flag"),
        "s_qcfail": ("qc_fail", "flag"),
        "s_nolayers": ("no_layers", "flag"),
    }
    ok = True
    for s, (exp_cat, exp_act) in expected.items():
        cat_ok = got.get(s) == exp_cat
        act_ok = got_action.get(s) == exp_act
        status = "OK" if (cat_ok and act_ok) else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"  [{status}] {s}: expected {exp_cat}/{exp_act}, "
              f"got {got.get(s)}/{got_action.get(s)}")

    # Exit code must be 1 (work remains).
    ec_ok = proc.returncode == 1
    print(f"  [{'OK' if ec_ok else 'FAIL'}] exit code: expected 1, got {proc.returncode}")

    # qc_fail reason surfaced.
    reason = next(r["reason"] for r in result["rows"] if r["sample"] == "s_qcfail")
    r_ok = "0.03" in reason
    print(f"  [{'OK' if r_ok else 'FAIL'}] qc_fail reason surfaced: {reason!r}")

    # Ledger: only s_done is DONE, so exactly one data row, with real metrics.
    led = pd.read_csv(ledger, sep="\t", dtype=str)
    led_ok = (list(led["sample"]) == ["s_done"]
              and led.iloc[0]["accession"] == "GSE_TEST"
              and led.iloc[0]["feature"] == "Gene"
              and led.iloc[0]["chemistry"] == "v3"
              and led.iloc[0]["reads_mapped_gene"] == "0.64"
              and led.iloc[0]["read_qc_status"] == "PASS")
    print(f"  [{'OK' if led_ok else 'FAIL'}] ledger records only s_done with metrics: "
          f"{list(led['sample'])}")

    # Ledger idempotency: a second identical run must append 0 rows.
    subprocess.run(cmd, capture_output=True, text=True)
    led2 = pd.read_csv(ledger, sep="\t", dtype=str)
    idem_ok = len(led2) == len(led) == 1
    print(f"  [{'OK' if idem_ok else 'FAIL'}] ledger idempotent on re-run: "
          f"{len(led)} -> {len(led2)} row(s)")

    print(f"\nreport.tsv:\n{report.read_text()}")
    if ok and ec_ok and r_ok and led_ok and idem_ok:
        print("\nALL FIXTURE ASSERTIONS PASSED")
        sys.exit(0)
    print("\nFIXTURE TEST FAILED")
    sys.exit(1)


if __name__ == "__main__":
    main()
