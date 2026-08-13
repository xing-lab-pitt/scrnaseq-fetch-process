#!/usr/bin/env python3
"""Deterministic read-QC flagger over per-SRR FastQC output.

FastQC runs on every sequencing run (SRR) but nothing reads its output — the
read-quality verdict is otherwise a human eyeballing the MultiQC HTML. This
script parses the FastQC zips and flags each SRR PASS / FAIL against FIXED,
config-supplied thresholds, so the read-quality decision is reproducible and
machine-readable (feeds the reconciler + the success ledger).

ROLE-AWARE (the important part). A 10x record has two reads:
  R1 = barcode+UMI read (16 bp CB + UMI, then poly-T). It is short and highly
       repetitive, so FastQC's per-base-content / GC / overrepresented / kmer
       modules FAIL *by design*. Judging R1 like a biological read floods false
       failures, so R1 gets only catastrophe checks (enough reads, sane length,
       quality not collapsed).
  R2 = cDNA read, the biologically meaningful one. Judged strictly: read count,
       mean base quality, adapter contamination.

Always exits 0 (this is the report layer, like qc_gate.py); enforcement / rerun
decisions live in the reconciler and the batch loop.

Report (qc/read_qc.tsv), one row per SRR:
  sample srr n_reads cdna_len cdna_mean_q cdna_pct_adapter cdna_pct_dup
  bc_len bc_mean_q flags status reasons
A sample rolls up as the WORST of its SRRs; qc/read_qc_pass.txt lists the samples
whose every SRR passed.
"""
import argparse
import csv
import io
import sys
import zipfile
from pathlib import Path

# FastQC modules that count as a real cDNA-read failure when FastQC marks them
# FAIL. Everything else (content/GC/overrepresented/kmer/duplication) is advisory
# for a cDNA read and ignored for a barcode read.
CDNA_FAIL_MODULES = {
    "per base sequence quality",
    "per sequence quality scores",
    "adapter content",
}


# --------------------------------------------------------------------------- #
# FastQC parsing
# --------------------------------------------------------------------------- #
def _read_member(zf, names, suffix):
    """Return the text of the single zip member ending in `suffix`, or None."""
    for n in names:
        if n.endswith(suffix):
            with zf.open(n) as fh:
                return io.TextIOWrapper(fh, encoding="utf-8", errors="replace").read()
    return None


def _parse_length(val):
    """'28' -> (28, 28); '35-76' -> (35, 76); junk -> (None, None)."""
    val = (val or "").strip()
    try:
        if "-" in val:
            lo, hi = val.split("-", 1)
            return int(lo), int(hi)
        n = int(val)
        return n, n
    except ValueError:
        return None, None


def parse_fastqc_zip(zip_path):
    """Parse one *_fastqc.zip into a flat metrics dict.

    Keys: modules {name.lower(): verdict}, n_reads, len_min, len_max, pct_gc,
    poor_quality, mean_q (avg of the per-base Mean column), max_adapter_frac
    (worst cumulative adapter % / 100), pct_dup (100 - total deduplicated %).
    Missing pieces come back as None. Raises on an unreadable zip so the caller
    can distinguish 'corrupt' from 'clean'.
    """
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        summary = _read_member(zf, names, "summary.txt")
        data = _read_member(zf, names, "fastqc_data.txt")

    modules = {}
    if summary:
        for line in summary.splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                modules[parts[1].strip().lower()] = parts[0].strip().upper()

    m = {
        "modules": modules, "n_reads": None, "len_min": None, "len_max": None,
        "pct_gc": None, "poor_quality": None, "mean_q": None,
        "max_adapter_frac": None, "pct_dup": None,
    }
    if not data:
        return m

    # Walk the module blocks (>>Name<TAB>verdict ... >>END_MODULE).
    block = None
    q_means, adapter_max = [], 0.0
    saw_adapter = False
    for line in data.splitlines():
        if line.startswith(">>END_MODULE"):
            block = None
            continue
        if line.startswith(">>"):
            block = line[2:].split("\t", 1)[0].strip().lower()
            continue
        if block == "basic statistics":
            if "\t" in line and not line.startswith("#"):
                k, _, v = line.partition("\t")
                k = k.strip()
                if k == "Total Sequences":
                    m["n_reads"] = int(float(v.strip()))
                elif k == "Sequence length":
                    m["len_min"], m["len_max"] = _parse_length(v)
                elif k == "%GC":
                    try:
                        m["pct_gc"] = float(v.strip())
                    except ValueError:
                        pass
                elif k == "Sequences flagged as poor quality":
                    try:
                        m["poor_quality"] = int(float(v.strip()))
                    except ValueError:
                        pass
        elif block == "per base sequence quality":
            if line and not line.startswith("#"):
                cols = line.split("\t")
                if len(cols) >= 2:
                    try:
                        q_means.append(float(cols[1]))
                    except ValueError:
                        pass
        elif block == "adapter content":
            if line and not line.startswith("#"):
                cols = line.split("\t")[1:]  # drop position column
                for c in cols:
                    try:
                        adapter_max = max(adapter_max, float(c))
                        saw_adapter = True
                    except ValueError:
                        pass
        elif block == "sequence duplication levels":
            if line.startswith("#Total Deduplicated Percentage"):
                try:
                    dedup = float(line.split("\t")[1])
                    m["pct_dup"] = round(100.0 - dedup, 2)
                except (ValueError, IndexError):
                    pass

    if q_means:
        m["mean_q"] = round(sum(q_means) / len(q_means), 2)
    if saw_adapter:
        m["max_adapter_frac"] = round(adapter_max / 100.0, 4)
    return m


# --------------------------------------------------------------------------- #
# Per-SRR evaluation
# --------------------------------------------------------------------------- #
def evaluate_srr(bc, cdna, thresholds):
    """Judge one SRR from parsed barcode (bc) and cDNA metrics dicts.

    bc/cdna may be None (zip missing) or a dict. `thresholds` keys:
    min_reads, min_mean_quality, max_adapter_fraction, min_barcode_len.
    Returns (row_dict, status, reasons_list).
    """
    min_reads = thresholds.get("min_reads")
    min_q = thresholds.get("min_mean_quality")
    max_adapter = thresholds.get("max_adapter_fraction")
    min_bc_len = thresholds.get("min_barcode_len", 24)

    reasons, flags = [], []

    if cdna is None:
        return ({}, "no_fastqc", ["cDNA-read FastQC output missing"])

    n_reads = cdna.get("n_reads")
    row = {
        "n_reads": n_reads,
        "cdna_len": cdna.get("len_max"),
        "cdna_mean_q": cdna.get("mean_q"),
        "cdna_pct_adapter": (None if cdna.get("max_adapter_frac") is None
                             else round(cdna["max_adapter_frac"] * 100, 2)),
        "cdna_pct_dup": cdna.get("pct_dup"),
        "bc_len": bc.get("len_max") if bc else None,
        "bc_mean_q": bc.get("mean_q") if bc else None,
    }

    # --- cDNA read (R2): strict ------------------------------------------- #
    if min_reads is not None and n_reads is not None and n_reads < min_reads:
        reasons.append(f"low_reads={n_reads}<{int(min_reads)}")
    if min_q is not None and cdna.get("mean_q") is not None and cdna["mean_q"] < min_q:
        reasons.append(f"low_cdna_quality={cdna['mean_q']}<{min_q}")
    if (max_adapter is not None and cdna.get("max_adapter_frac") is not None
            and cdna["max_adapter_frac"] > max_adapter):
        reasons.append(f"adapter={cdna['max_adapter_frac']:.3g}>{max_adapter}")
    for mod, verdict in (cdna.get("modules") or {}).items():
        if verdict == "FAIL" and mod in CDNA_FAIL_MODULES:
            flags.append(f"cdna:{mod}=FAIL")

    # --- barcode read (R1): catastrophe checks only ----------------------- #
    if bc is not None:
        bc_reads = bc.get("n_reads")
        if (n_reads and bc_reads and
                abs(bc_reads - n_reads) > 0.001 * max(bc_reads, n_reads)):
            reasons.append(f"read_count_mismatch bc={bc_reads} cdna={n_reads}")
        if bc.get("len_max") is not None and bc["len_max"] < min_bc_len:
            reasons.append(f"barcode_too_short={bc['len_max']}<{min_bc_len}")
        # Barcode quality must not be catastrophically low (5 below the cDNA
        # floor); its content/GC/overrepresented FAILs are expected, so ignored.
        if (min_q is not None and bc.get("mean_q") is not None
                and bc["mean_q"] < min_q - 5):
            reasons.append(f"barcode_quality_collapsed={bc['mean_q']}")
    else:
        flags.append("bc:no_fastqc")

    # A meaningful cDNA module FAIL escalates to FAIL only if it is not already
    # explained by a numeric reason (avoid double-counting adapter).
    if not reasons and flags:
        hard = [f for f in flags if f.startswith("cdna:")]
        if hard:
            reasons.append("fastqc_fail:" + ",".join(
                f.split(":", 1)[1] for f in hard))

    row["flags"] = ";".join(flags)
    status = "PASS" if not reasons else "FAIL"
    return row, status, reasons


# --------------------------------------------------------------------------- #
# Sheet + discovery
# --------------------------------------------------------------------------- #
def srr_to_sample(samples_tsv):
    """{srr: sample} from the sheet (both columns required)."""
    out = {}
    with open(samples_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "srr" not in reader.fieldnames \
                or "sample" not in reader.fieldnames:
            raise ValueError(f"{samples_tsv}: needs 'sample' and 'srr' columns")
        for row in reader:
            srr = (row.get("srr") or "").strip()
            sample = (row.get("sample") or "").strip()
            if srr:
                out[srr] = sample
    return out


def find_zip(fastqc_dir, srr, read):
    p = Path(fastqc_dir) / f"{srr}_R{read}_fastqc.zip"
    return p if p.exists() else None


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
REPORT_FIELDS = ["sample", "srr", "n_reads", "cdna_len", "cdna_mean_q",
                 "cdna_pct_adapter", "cdna_pct_dup", "bc_len", "bc_mean_q",
                 "flags", "status", "reasons"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastqc-dir", required=True, help="dir of *_fastqc.zip files")
    ap.add_argument("--samples", required=True, help="samples.tsv (srr -> sample)")
    ap.add_argument("--report", required=True)
    ap.add_argument("--passlist", required=True)
    ap.add_argument("--min-reads", type=float)
    ap.add_argument("--min-mean-quality", type=float)
    ap.add_argument("--max-adapter-fraction", type=float)
    ap.add_argument("--min-barcode-len", type=int, default=24)
    args = ap.parse_args()

    thresholds = {
        "min_reads": args.min_reads,
        "min_mean_quality": args.min_mean_quality,
        "max_adapter_fraction": args.max_adapter_fraction,
        "min_barcode_len": args.min_barcode_len,
    }

    s2s = srr_to_sample(args.samples)

    rows = []
    for srr, sample in s2s.items():
        z_bc = find_zip(args.fastqc_dir, srr, 1)   # R1 = barcode read
        z_cd = find_zip(args.fastqc_dir, srr, 2)   # R2 = cDNA read
        try:
            bc = parse_fastqc_zip(z_bc) if z_bc else None
            cdna = parse_fastqc_zip(z_cd) if z_cd else None
        except (zipfile.BadZipFile, OSError) as e:
            rows.append({"sample": sample, "srr": srr, "status": "no_fastqc",
                         "reasons": f"unreadable FastQC zip: {e}", "flags": ""})
            continue
        row, status, reasons = evaluate_srr(bc, cdna, thresholds)
        row.update({"sample": sample, "srr": srr, "status": status,
                    "reasons": ";".join(reasons)})
        rows.append(row)

    rows.sort(key=lambda r: (r.get("sample", ""), r.get("srr", "")))

    # Sample rolls up as the WORST of its SRRs: a sample passes iff EVERY SRR
    # has status PASS (a no_fastqc SRR keeps the sample out of the pass list).
    by_sample = {}
    for r in rows:
        by_sample.setdefault(r["sample"], []).append(r.get("status"))
    passing = sorted(s for s, sts in by_sample.items()
                     if sts and all(x == "PASS" for x in sts))

    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    Path(args.passlist).write_text("\n".join(passing) + ("\n" if passing else ""))

    n_fail = sum(1 for r in rows if r.get("status") != "PASS")
    print(f"read QC: {len(passing)}/{len(by_sample)} samples passed "
          f"({len(rows)} sequencing runs; {n_fail} run(s) flagged). "
          f"Report: {args.report}")
    for r in rows:
        if r.get("status") != "PASS":
            print(f"  {r.get('status')} {r['sample']}/{r['srr']}: "
                  f"{r.get('reasons') or r.get('flags')}", file=sys.stderr)
    sys.exit(0)  # report always succeeds; enforcement is separate


if __name__ == "__main__":
    main()
