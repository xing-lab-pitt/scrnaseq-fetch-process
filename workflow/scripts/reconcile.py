#!/usr/bin/env python3
"""Completeness reconciler for the scRNA-seq STARsolo pipeline.

Snakemake treats a sample as "done" when its per-sample .h5ad merely EXISTS.
That is not enough: the h5ad can exist while the sample FAILED the QC gate, or
while the Velocyto layers are silently absent (starsolo_to_h5ad.py writes zero
layers, without error, when Solo.out/Velocyto/ is missing). This script enforces
the fuller contract and produces the to-do list of what still needs (re)running.

A sample is DONE iff ALL of:
  1. <workdir>/h5ad/<sample>.h5ad exists, AND
  2. <sample> is listed in <workdir>/qc/qc_pass.txt (passed the QC gate), AND
  3. that .h5ad has layers spliced, unspliced, ambiguous.

Otherwise it is categorized, and each category carries an `action` that tells a
driver whether to RERUN it (transient — hand back to Snakemake, worth retrying) or
FLAG it (genuine — a rerun of the same inputs can't fix it, so a human decides):
  missing      - no h5ad (not run yet, or an upstream step crashed).      -> RERUN
  corrupt      - h5ad exists but is unreadable (truncated / killed write). -> RERUN
                 (the batch loop quarantines it first so Snakemake rebuilds it.)
  read_qc_fail - the raw reads failed read-QC (bad library / wrong
                 chemistry). The same reads re-aligned give the same bad
                 result, so this is a human call, not a retry.             -> FLAG
  qc_fail      - h5ad but not in qc_pass.txt (failed the alignment gate);
                 reason from qc_gate.tsv.                                  -> FLAG
  no_layers    - h5ad + QC pass but a velocity layer absent (mis-wire).    -> FLAG

Read-QC axis: if qc/read_qc.tsv exists (written by read_qc.py), a sample whose
reads genuinely FAILED read-QC is categorized read_qc_fail up front — rerunning
bad input reads is pointless. A "no_fastqc" read-QC status does NOT block (FastQC
is a separate branch); it just annotates the row. Absent read_qc.tsv -> axis
skipped entirely, so pre-read-QC workdirs reconcile exactly as before.

Two modes:
  Single workdir:  reconcile.py --samples config/samples.tsv --workdir <dir>
                   (both default from config/config.yaml if --config given/found)
  Manifest:        reconcile.py --manifest config/manifest.tsv

Success logbook (--ledger PATH): append the samples that are DONE this pass to an
append-only TSV (results/successful_samples.tsv), one row per (accession, sample),
skipping any pair already recorded. It is the durable record of what mapped
successfully, with the mapping metrics pulled from qc_gate.tsv.

Exit code: 0 when the to-do list is empty (everything DONE), 1 when work remains,
2 on a usage/IO error. The layer check reads only the /layers HDF5 group keys via
h5py -- O(1), no matrix load.
"""
import argparse
import csv
import json
import sys
from datetime import datetime
from pathlib import Path

REQUIRED_LAYERS = ("spliced", "unspliced", "ambiguous")

# Category ordering for summaries / exit decisions.
CATEGORIES = ("done", "missing", "corrupt", "read_qc_fail", "qc_fail", "no_layers")
# Categories that mean "not done yet".
TODO_CATEGORIES = ("missing", "corrupt", "read_qc_fail", "qc_fail", "no_layers")
# Transient -> hand back to Snakemake (worth an automatic retry).
RERUN_CATEGORIES = ("missing", "corrupt")
# Genuine -> a rerun of the same inputs can't fix it; a human decides.
FLAG_CATEGORIES = ("read_qc_fail", "qc_fail", "no_layers")


def category_action(category):
    """rerun (transient), flag (genuine), or none (done)."""
    if category == "done":
        return "none"
    return "rerun" if category in RERUN_CATEGORIES else "flag"


# --------------------------------------------------------------------------- #
# Config / sample-sheet loading
# --------------------------------------------------------------------------- #
def load_config(config_path):
    """Read config/config.yaml -> dict (best-effort; only keys we need)."""
    import yaml
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def expected_samples(samples_tsv):
    """Deduped, sorted sample names from the sheet's `sample` column.

    Mirrors Snakefile: SAMPLES = sorted(samples.groupby('sample')...).
    """
    path = Path(samples_tsv)
    if not path.exists():
        raise FileNotFoundError(f"samples sheet not found: {samples_tsv}")
    samples = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError(f"{samples_tsv}: no 'sample' column")
        for row in reader:
            s = (row.get("sample") or "").strip()
            if s:
                samples.append(s)
    return sorted(set(samples))


# --------------------------------------------------------------------------- #
# On-disk state readers
# --------------------------------------------------------------------------- #
def read_pass_set(workdir):
    """Sample names that passed the QC gate (qc/qc_pass.txt), or None if the
    gate has not run yet (file absent)."""
    p = Path(workdir) / "qc" / "qc_pass.txt"
    if not p.exists():
        return None
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def read_fail_reasons(workdir):
    """{sample: reason} for FAIL rows in qc/qc_gate.tsv (empty if absent)."""
    p = Path(workdir) / "qc" / "qc_gate.tsv"
    reasons = {}
    if not p.exists():
        return reasons
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("status") or "").strip().upper() == "FAIL":
                reasons[row["sample"].strip()] = (row.get("reasons") or "").strip()
    return reasons


def read_qc_metrics(workdir):
    """{sample: {valid_barcodes, reads_mapped_gene, estimated_cells, saturation}}
    for every row of qc/qc_gate.tsv (empty if absent). Feeds the success ledger."""
    p = Path(workdir) / "qc" / "qc_gate.tsv"
    out = {}
    if not p.exists():
        return out
    keep = ("valid_barcodes", "reads_mapped_gene", "estimated_cells", "saturation")
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if s:
                out[s] = {k: (row.get(k) or "").strip() for k in keep}
    return out


def read_readqc_status(workdir):
    """{sample: {status, reason, n_reads}} rolled up from qc/read_qc.tsv, or None if
    absent (axis not evaluated). Per-sample roll-up = WORST SRR: FAIL > no_fastqc >
    PASS; n_reads = summed cDNA reads across the sample's SRRs (for the ledger)."""
    p = Path(workdir) / "qc" / "read_qc.tsv"
    if not p.exists():
        return None
    rank = {"PASS": 0, "no_fastqc": 1, "FAIL": 2}
    agg = {}
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if not s:
                continue
            st = (row.get("status") or "").strip()
            reason = (row.get("reasons") or "").strip()
            rec = agg.setdefault(s, {"status": "", "reason": "", "n_reads": 0})
            if rank.get(st, 0) >= rank.get(rec["status"], -1):
                rec["status"], rec["reason"] = st, reason
            try:
                rec["n_reads"] += int(float(row.get("n_reads") or 0))
            except ValueError:
                pass
    return agg


def sample_meta(samples_tsv):
    """{sample: {srrs: [...], chemistry: str}} from the sheet (for the ledger)."""
    out = {}
    path = Path(samples_tsv)
    if not path.exists():
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if not s:
                continue
            rec = out.setdefault(s, {"srrs": [], "chemistry": ""})
            srr = (row.get("srr") or "").strip()
            if srr:
                rec["srrs"].append(srr)
            chem = (row.get("chemistry") or "").strip()
            if chem and not rec["chemistry"]:
                rec["chemistry"] = chem
    return out


def h5ad_layers(h5ad_path):
    """Layer names present in an .h5ad, read cheaply from the /layers HDF5 group.

    No matrix is loaded. Returns a set (empty if the group is absent). Raises on a
    corrupt/unreadable file so the caller can distinguish it from 'no layers'.
    """
    import h5py
    with h5py.File(h5ad_path, "r") as f:
        grp = f.get("layers")
        if grp is None:
            return set()
        return set(grp.keys())


# --------------------------------------------------------------------------- #
# Core per-sample evaluation
# --------------------------------------------------------------------------- #
def evaluate_sample(sample, workdir, pass_set, fail_reasons, readqc=None):
    """Return a dict: category + evidence for one sample.

    `readqc` is the read-QC roll-up dict (from read_readqc_status) or None when
    read_qc.tsv is absent (axis skipped for back-compat).
    """
    h5ad = Path(workdir) / "h5ad" / f"{sample}.h5ad"
    rec = {
        "sample": sample,
        "h5ad": str(h5ad),
        "has_h5ad": h5ad.exists(),
        "qc_pass": None,
        "read_qc": "",
        "layers_present": "",
        "category": None,
        "reason": "",
    }

    # Read-QC axis (input-read quality; annotates every row when evaluated).
    if readqc is not None:
        rq = readqc.get(sample, {"status": "no_fastqc", "reason": "no read_qc row"})
        rec["read_qc"] = rq.get("status", "")
        # A genuine read-QC FAIL is the root cause: rerunning the SAME reads gives
        # the same bad result, so flag for a human up front rather than looping it
        # through the rerun path. (no_fastqc does NOT block — FastQC is a separate
        # branch and may simply not have run for an otherwise-good sample.)
        if rec["read_qc"] == "FAIL":
            rec["category"] = "read_qc_fail"
            rec["reason"] = rq.get("reason") or "reads failed read-QC"
            return rec

    if not rec["has_h5ad"]:
        rec["category"] = "missing"
        rec["reason"] = "h5ad not produced yet (or upstream step crashed)"
        return rec

    # h5ad exists -> is it a QC pass? (pass_set is None until the gate runs)
    if pass_set is None:
        rec["category"] = "missing"
        rec["reason"] = "QC gate has not run (qc/qc_pass.txt absent)"
        return rec
    rec["qc_pass"] = sample in pass_set
    if not rec["qc_pass"]:
        rec["category"] = "qc_fail"
        rec["reason"] = fail_reasons.get(sample, "not in qc_pass.txt")
        return rec

    # QC passed -> verify the velocity layers actually landed.
    try:
        layers = h5ad_layers(str(h5ad))
    except Exception as e:  # unreadable file -> transient write failure, worth a retry
        rec["category"] = "corrupt"
        rec["reason"] = f"h5ad unreadable: {e}"
        return rec
    rec["layers_present"] = ",".join(sorted(layers))
    missing_layers = [l for l in REQUIRED_LAYERS if l not in layers]
    if missing_layers:
        rec["category"] = "no_layers"
        rec["reason"] = f"missing layers: {','.join(missing_layers)}"
        return rec

    rec["category"] = "done"
    return rec


def reconcile_workdir(samples_tsv, workdir, accession="", feature=""):
    """Evaluate every expected sample in one workdir. Returns a study result dict."""
    samples = expected_samples(samples_tsv)
    pass_set = read_pass_set(workdir)
    fail_reasons = read_fail_reasons(workdir)
    readqc = read_readqc_status(workdir)   # None if read_qc.tsv absent (axis skipped)
    rows = [evaluate_sample(s, workdir, pass_set, fail_reasons, readqc)
            for s in samples]
    for r in rows:
        r["accession"] = accession
        r["action"] = category_action(r["category"])

    counts = {c: 0 for c in CATEGORIES}
    for r in rows:
        counts[r["category"]] += 1

    combined = Path(workdir) / "combined.h5ad"
    return {
        "accession": accession,
        "feature": feature,
        "workdir": str(workdir),
        "samples_tsv": str(samples_tsv),
        "n_expected": len(samples),
        "counts": counts,
        "combined_h5ad": combined.exists(),
        "complete": all(r["category"] == "done" for r in rows) and len(rows) > 0,
        "rows": rows,
    }


# --------------------------------------------------------------------------- #
# Manifest mode
# --------------------------------------------------------------------------- #
def read_manifest(manifest_path, base_dir):
    """Rows of config/manifest.tsv. Resolves per-study samples_tsv + workdir.

    Columns: accession, feature, workdir[, samples_tsv, notes]. workdir and
    samples_tsv are resolved relative to base_dir (the pipeline root) if not
    absolute. samples_tsv defaults to <workdir>/samples.used.tsv then
    config/<accession>.samples.tsv then config/samples.tsv.
    """
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        raise FileNotFoundError(f"manifest not found: {manifest_path}")
    base = Path(base_dir)
    out = []
    with open(manifest_path, newline="") as fh:
        # Drop whole-line comments (incl. those ABOVE the header) so the first
        # real line becomes csv's header — otherwise a leading '#' line is taken
        # as the header and every data row is silently skipped.
        lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        acc = (row.get("accession") or "").strip()
        if not acc or acc.startswith("#"):
            continue
        workdir = (row.get("workdir") or "").strip()
        if not workdir:
            raise ValueError(f"manifest row {acc}: missing workdir")
        wd = Path(workdir)
        if not wd.is_absolute():
            wd = base / wd
        samples_tsv = (row.get("samples_tsv") or "").strip()
        if samples_tsv:
            st = Path(samples_tsv)
            if not st.is_absolute():
                st = base / st
        else:
            candidates = [
                wd / "samples.used.tsv",
                base / "config" / f"{acc}.samples.tsv",
                base / "config" / "samples.tsv",
            ]
            st = next((c for c in candidates if c.exists()), candidates[-1])
        out.append({
            "accession": acc,
            "feature": (row.get("feature") or "").strip(),
            "workdir": str(wd),
            "samples_tsv": str(st),
            "notes": (row.get("notes") or "").strip(),
        })
    return out


def reconcile_manifest(manifest_path, base_dir):
    entries = read_manifest(manifest_path, base_dir)
    studies = []
    for e in entries:
        try:
            res = reconcile_workdir(e["samples_tsv"], e["workdir"],
                                    accession=e["accession"], feature=e["feature"])
        except FileNotFoundError as exc:
            res = {
                "accession": e["accession"], "feature": e["feature"],
                "workdir": e["workdir"],
                "samples_tsv": e["samples_tsv"], "n_expected": 0,
                "counts": {c: 0 for c in CATEGORIES}, "combined_h5ad": False,
                "complete": False, "rows": [],
                "error": str(exc),
            }
        studies.append(res)
    complete = bool(studies) and all(s["complete"] for s in studies)
    return {"mode": "manifest", "complete": complete, "studies": studies}


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
REPORT_FIELDS = ["accession", "sample", "category", "action", "has_h5ad",
                 "qc_pass", "read_qc", "layers_present", "reason", "h5ad"]

# Success logbook: append-only, one row per (accession, sample) that is DONE.
LEDGER_FIELDS = ["recorded_at", "accession", "sample", "srrs", "chemistry",
                 "feature", "n_reads", "valid_barcodes", "reads_mapped_gene",
                 "estimated_cells", "saturation", "read_qc_status", "layers",
                 "workdir", "h5ad"]


def all_rows(result):
    if result.get("mode") == "manifest":
        for s in result["studies"]:
            yield from s["rows"]
    else:
        yield from result["rows"]


def write_report(result, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in all_rows(result):
            w.writerow(r)


def study_studies(result):
    """The per-study result list, for both manifest and single-workdir mode."""
    return result["studies"] if result.get("mode") == "manifest" else [result]


def print_summary(result):
    def study_line(s):
        c = s["counts"]
        tag = "OK " if s["complete"] else "TODO"
        err = f"  ERROR: {s['error']}" if s.get("error") else ""
        combined = "combined.h5ad" if s["combined_h5ad"] else "no-combined"
        return (f"[{tag}] {s['accession'] or s['workdir']}: "
                f"{c['done']}/{s['n_expected']} done "
                f"(missing={c['missing']} corrupt={c['corrupt']} "
                f"read_qc_fail={c['read_qc_fail']} qc_fail={c['qc_fail']} "
                f"no_layers={c['no_layers']}) [{combined}]{err}")

    print("=== reconcile summary ===")
    for s in study_studies(result):
        print(study_line(s))

    # Detail the non-done samples so a human sees WHY — and, crucially, whether the
    # driver should RERUN it or a human should look at it (the `action`).
    todo = [r for r in all_rows(result) if r["category"] != "done"]
    if todo:
        print("\n--- to-do (sample: category [action] — reason) ---")
        for r in todo:
            acc = f"{r['accession']}/" if r.get("accession") else ""
            print(f"  {acc}{r['sample']}: {r['category']} "
                  f"[{r.get('action', '?')}] — {r['reason']}")
    print(f"\nOverall: {'COMPLETE' if result['complete'] else 'INCOMPLETE'}")


# --------------------------------------------------------------------------- #
# Success logbook (append-only)
# --------------------------------------------------------------------------- #
def _ledger_keys(path):
    """Set of (accession, sample) already recorded, for idempotent appends."""
    keys = set()
    p = Path(path)
    if not p.exists():
        return keys
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            keys.add(((row.get("accession") or "").strip(),
                      (row.get("sample") or "").strip()))
    return keys


def write_ledger(result, path, recorded_at):
    """Append every DONE (accession, sample) not already in the ledger. Returns the
    number of new rows appended. Append-only + keyed, so re-running is a no-op for
    samples already recorded (idempotent)."""
    existing = _ledger_keys(path)
    new_rows = []
    for s in study_studies(result):
        meta = sample_meta(s["samples_tsv"])
        metrics = read_qc_metrics(s["workdir"])
        readqc = read_readqc_status(s["workdir"]) or {}
        feature = s.get("feature", "")
        acc = s.get("accession", "")
        for r in s["rows"]:
            if r["category"] != "done":
                continue
            if (acc, r["sample"]) in existing:
                continue
            existing.add((acc, r["sample"]))  # guard against dup rows within a run
            m = meta.get(r["sample"], {})
            qm = metrics.get(r["sample"], {})
            rq = readqc.get(r["sample"], {})
            new_rows.append({
                "recorded_at": recorded_at,
                "accession": acc,
                "sample": r["sample"],
                "srrs": ",".join(m.get("srrs", [])),
                "chemistry": m.get("chemistry", ""),
                "feature": feature,
                "n_reads": rq.get("n_reads", ""),
                "valid_barcodes": qm.get("valid_barcodes", ""),
                "reads_mapped_gene": qm.get("reads_mapped_gene", ""),
                "estimated_cells": qm.get("estimated_cells", ""),
                "saturation": qm.get("saturation", ""),
                "read_qc_status": rq.get("status", ""),
                "layers": r.get("layers_present", ""),
                "workdir": s["workdir"],
                "h5ad": r.get("h5ad", ""),
            })

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    write_header = not p.exists()
    with open(p, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        if write_header:
            w.writeheader()
        for row in new_rows:
            w.writerow(row)
    return len(new_rows)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="Reconcile every study in this manifest TSV.")
    ap.add_argument("--samples", help="samples.tsv (single-workdir mode).")
    ap.add_argument("--workdir", help="pipeline workdir (single-workdir mode).")
    ap.add_argument("--config", help="config/config.yaml to fill missing "
                    "--samples/--workdir (single-workdir mode).")
    ap.add_argument("--base-dir", default=".",
                    help="pipeline root for resolving manifest relative paths "
                         "(default: cwd).")
    ap.add_argument("--accession", default="",
                    help="accession label for the ledger (single-workdir mode); "
                         "manifest mode reads it per study.")
    ap.add_argument("--feature", default="",
                    help="feature recorded in the ledger (single-workdir mode); "
                         "manifest mode reads it per study.")
    ap.add_argument("--report", help="write per-sample TSV report here.")
    ap.add_argument("--json", dest="json_out", help="write machine-readable JSON here.")
    ap.add_argument("--ledger", help="append DONE samples to this success logbook "
                    "TSV (append-only, idempotent per accession+sample).")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="suppress the human summary (still writes report/json).")
    args = ap.parse_args()

    try:
        if args.manifest:
            result = reconcile_manifest(args.manifest, args.base_dir)
        else:
            samples, workdir = args.samples, args.workdir
            if (not samples or not workdir) and args.config:
                cfg = load_config(args.config)
                samples = samples or cfg.get("samples_tsv")
                workdir = workdir or cfg.get("workdir")
            if not samples or not workdir:
                ap.error("single-workdir mode needs --samples and --workdir "
                         "(or --config to supply them)")
            result = reconcile_workdir(samples, workdir,
                                       accession=args.accession, feature=args.feature)
    except (FileNotFoundError, ValueError) as e:
        print(f"reconcile: {e}", file=sys.stderr)
        sys.exit(2)

    if args.report:
        write_report(result, args.report)
    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(json.dumps(result, indent=2))
    if args.ledger:
        stamp = datetime.now().isoformat(timespec="seconds")
        n = write_ledger(result, args.ledger, stamp)
        if not args.quiet:
            print(f"ledger: +{n} newly-recorded sample(s) -> {args.ledger}")
    if not args.quiet:
        print_summary(result)

    sys.exit(0 if result["complete"] else 1)


if __name__ == "__main__":
    main()
