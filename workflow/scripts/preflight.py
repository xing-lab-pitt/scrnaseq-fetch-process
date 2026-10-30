#!/usr/bin/env python3
"""Pre-alignment read diagnostics: predict mapping loss before spending an alignment.

Motivation. On GSE310923, 16.8% of cDNA reads began with the 10x template-switch oligo
(TSO), leaving 50 bp of real sequence against STAR's default 53 bp minimum. Those reads
aligned fine but were discarded by a length filter, costing ~12 points of gene mapping.
It took a 6 h alignment to notice. Every input needed to predict it was in the FASTQ.

FastQC does not catch this: its Adapter Content module knows only Illumina/Nextera/SOLiD
adapters, and its Overrepresented Sequences module counts exact strings, so one oligo
followed by many different cDNAs fragments into entries below the reporting threshold
(it reported ~0.9% against a true 16.8%). Anchored prefix matching is what finds it.

The core check is a consistency test between the reads and the aligner settings, not a
hardcoded rule: STAR drops a read when matched length < outFilterMatchNminOverLread x
read_len, so the same arithmetic run over the reads predicts the "too short" bucket.

Read-only and report-only; never blocks a DAG.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_runs import (  # noqa: E402  (path set above)
    DEFAULT_WHITELIST_DIR,
    classify_chemistry,
    match_whitelists,
)

# 5' oligos that leave a read unalignable at its start. The TSO is the 10x/SMARTer
# template-switch oligo; template switching adds a variable G run past its 3' end, so
# the real prefix is 30-36 bp, not a constant (see --report-boundary).
DEFAULT_OLIGOS = {
    "tso": "AAGCAGTGGTATCAACGCAGAGTACATGGG",
    "illumina_universal": "AGATCGGAAGAGC",
}
HOMOPOLYMER_MIN = 20            # run length that marks a read as poly-A/T/G
LOWCX_FRAC = 0.70               # one base at/above this share = low complexity
STAR_DEFAULT_MATCH_NMIN = 0.66  # STAR's own default when nothing overrides it

REPORT_FIELDS = [
    "sample", "srr", "n_sampled", "bc_len", "cdna_len", "chem_by_len",
    "wl_chem", "wl_frac", "tso_pct", "tso_ragged_pct", "illumina_universal_pct",
    "polyg_pct", "polya_pct", "polyt_pct", "lowcx_pct",
    "star_min_match", "pred_too_short_pct", "irrecoverable_pct", "sweep",
    "status", "flags", "notes",
]


def read_seqs(path, n_reads):
    """Sequence lines from the head of a gzipped FASTQ. Head, not random: enough to
    measure a percent-level effect, and it is what the aligner reads first anyway."""
    seqs = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.rstrip("\n"))
                if len(seqs) >= n_reads:
                    break
    return seqs


def modal_len(seqs):
    return Counter(len(s) for s in seqs).most_common(1)[0][0] if seqs else 0


def _trailing_run(seq, start, base):
    """Length of the run of `base` beginning at `start`."""
    n = 0
    while start + n < len(seq) and seq[start + n] == base:
        n += 1
    return n


def prefix_len(seq, oligos, anchor=25):
    """Leading unalignable length of one read.

    Matches each oligo anchored at position 0 using its first `anchor` bases (tolerating
    sequencing error in the tail), then absorbs any homopolymer run past the oligo — the
    variable G tail template switching leaves behind. Returns (length, oligo_name, ragged)
    where ragged means extra bases were absorbed beyond the nominal oligo.
    """
    for name, oligo in oligos.items():
        probe = oligo[:anchor]
        if seq.startswith(probe):
            end = len(oligo)
            extra = _trailing_run(seq, end, oligo[-1])
            return end + extra, name, extra > 0
    return 0, None, False


SWEEP_GRID = (0.66, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30)


def sweep_filter(usable_hist, read_len, grid=SWEEP_GRID):
    """Predicted 'too short' loss at each candidate outFilterMatchNminOverLread.

    Deliberately reports the whole curve rather than naming one value. A single
    recommendation hides the trade-off and is easy to over-trust: on GSE310923 the
    quantile-based pick was 0.50, while the actual A/B showed 0.30 recovered ~1.2 points
    more with no sign of admitting noise. The reads can say how much is lost at each
    setting; only an alignment can say whether what is gained is real. So: show the
    curve, let a human choose, confirm with strand_probe/an A/B if it matters.

    Returns [(T, predicted_too_short_pct), ...].
    """
    n = sum(usable_hist.values())
    if not n or not read_len:
        return []
    out = []
    for t in grid:
        thresh = math.ceil(t * read_len)
        lost = sum(c for u, c in usable_hist.items() if u < thresh)
        out.append((t, round(100.0 * lost / n, 2)))
    return out


def scan(seqs, oligos, min_match_frac):
    """Classify a cDNA read sample and predict STAR's too-short bucket."""
    n = len(seqs)
    if not n:
        return {}
    read_len = modal_len(seqs)
    # The threshold STAR will actually apply, derived from the configured filter.
    star_thresh = math.ceil(min_match_frac * read_len)

    counts = Counter()
    ragged = 0
    pred_short = 0
    homo = {"A": 0, "T": 0, "G": 0}
    usable_hist = Counter()   # alignable length per read -> recommended filter value

    for s in seqs:
        plen, name, is_ragged = prefix_len(s, oligos)
        if name:
            counts[name] += 1
            ragged += is_ragged

        usable = len(s) - plen

        # Whole-read artifacts leave nothing alignable regardless of prefix.
        for base in ("A", "T", "G"):
            if base * HOMOPOLYMER_MIN in s:
                homo[base] += 1
                usable = 0
        bases = Counter(s)
        if bases and max(bases.values()) / len(s) >= LOWCX_FRAC:
            counts["lowcx"] += 1
            usable = 0

        usable_hist[usable] += 1
        if usable < star_thresh:
            pred_short += 1

    pct = lambda k: round(100.0 * k / n, 3)
    sweep = sweep_filter(usable_hist, read_len)
    out = {
        "sweep": ";".join(f"{t:.2f}:{v}" for t, v in sweep),
        "_sweep": sweep,
        # Floor of the curve: reads with no alignable sequence at all (poly-A, low
        # complexity). No threshold recovers these, so the sweep cannot go below it.
        "irrecoverable_pct": round(100.0 * usable_hist.get(0, 0) / n, 3),
        "n_sampled": n,
        "cdna_len": read_len,
        "star_min_match": star_thresh,
        "tso_pct": pct(counts["tso"]),
        "tso_ragged_pct": pct(ragged),
        "illumina_universal_pct": pct(counts["illumina_universal"]),
        "polyg_pct": pct(homo["G"]),
        "polya_pct": pct(homo["A"]),
        "polyt_pct": pct(homo["T"]),
        "lowcx_pct": pct(counts["lowcx"]),
        "pred_too_short_pct": pct(pred_short),
    }
    return out


def judge(row, thr):
    """Report-only verdict. Returns (status, flags, notes).

    Only things that cost data AT THE CURRENT SETTINGS become flags. Contamination the
    configured filter already tolerates is a note: a gate that keeps shouting about a TSO
    after the filter has been fixed trains people to ignore it, and the one time it means
    something they will. `pred_too_short` is the load-bearing number because it is the
    one that folds the reads and the aligner settings together.
    """
    flags, notes = [], []
    if row.get("wl_frac") is not None and row["wl_frac"] < thr["min_whitelist_frac"]:
        flags.append(f"barcode_whitelist={row['wl_frac']:.3f}<{thr['min_whitelist_frac']}")
    if row.get("bc_len") and row["bc_len"] not in (24, 26, 27, 28):
        flags.append(f"barcode_len={row['bc_len']}_unexpected")

    costly = row.get("pred_too_short_pct", 0) > thr["max_pred_too_short_pct"]
    if costly:
        flags.append(
            f"pred_too_short={row['pred_too_short_pct']}%>{thr['max_pred_too_short_pct']}%")
        flags.append("see filter sweep below: pick --outFilterMatchNminOverLread")

    # Prefix contamination: actionable only while it is actually costing reads.
    for key in ("tso_pct", "illumina_universal_pct", "polyg_pct", "lowcx_pct"):
        if row.get(key, 0) > thr["max_prefix_pct"]:
            (flags if costly else notes).append(f"{key}={row[key]}%")
    if row.get("tso_ragged_pct", 0) > 1.0:
        notes.append("ragged_prefix_boundary:do_not_hard_clip")
    return ("FLAG" if flags else "OK"), ";".join(flags), ";".join(notes)


def star_min_match_from_config(cfg):
    """The outFilterMatchNminOverLread STARsolo will really use.

    Reading it from starsolo_extra (rather than assuming STAR's default) is what makes
    the prediction track the pipeline: relax the filter and the predicted bucket falls
    with it, instead of the check going stale.
    """
    extra = (cfg or {}).get("starsolo_extra", "") or ""
    m = re.search(r"--outFilterMatchNminOverLread\s+([0-9.]+)", extra)
    return float(m.group(1)) if m else STAR_DEFAULT_MATCH_NMIN


# Flags the starsolo rule sets itself from other config keys; setting one of these
# again inside starsolo_extra makes STAR receive it twice, which is a config-authoring
# bug, not a data problem, so this check exits instead of just flagging like judge() does.
FORBIDDEN_IN_STARSOLO_EXTRA = ("--soloFeatures", "--soloCellFilter")


def check_starsolo_extra_double_set(cfg):
    """Refuse a starsolo_extra that duplicates a flag the starsolo rule sets itself.

    --soloFeatures comes from the `feature` config key and --soloCellFilter from
    `starsolo_cell_filter`; STARsolo would otherwise receive either flag twice.
    """
    extra = (cfg or {}).get("starsolo_extra", "") or ""
    hits = [f for f in FORBIDDEN_IN_STARSOLO_EXTRA if f in extra]
    if hits:
        print(f"ERROR: starsolo_extra sets {', '.join(hits)}, which the starsolo rule "
              "already sets from a dedicated config key.", file=sys.stderr)
        print("Remove it from starsolo_extra to avoid STAR receiving it twice.",
              file=sys.stderr)
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir", help="pipeline workdir (expects <workdir>/fastq)")
    ap.add_argument("--samples", required=True, help="samples.tsv (for sample<-srr map)")
    ap.add_argument("--config", help="config.yaml (reads starsolo_extra + preflight block)")
    ap.add_argument("--srr", action="append", help="limit to these SRRs (repeatable)")
    ap.add_argument("--n-reads", type=int, default=1_000_000)
    ap.add_argument("--whitelist-dir", default=DEFAULT_WHITELIST_DIR)
    ap.add_argument("-o", "--output", help="TSV out (default: <workdir>/qc/preflight.tsv)")
    args = ap.parse_args()

    cfg = {}
    if args.config:
        import yaml
        cfg = yaml.safe_load(open(args.config)) or {}
    check_starsolo_extra_double_set(cfg)
    pf = cfg.get("preflight", {}) or {}
    oligos = pf.get("oligos") or DEFAULT_OLIGOS
    thr = {
        "min_whitelist_frac": pf.get("min_whitelist_frac", 0.70),
        "max_prefix_pct": pf.get("max_prefix_pct", 5.0),
        "max_pred_too_short_pct": pf.get("max_pred_too_short_pct", 10.0),
    }
    min_match = star_min_match_from_config(cfg)

    wd = Path(args.workdir)
    fq = wd / "fastq"
    with open(args.samples) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    srr_sample = {r["srr"]: r["sample"] for r in rows}
    targets = args.srr or [r["srr"] for r in rows]

    out_rows = []
    for srr in targets:
        r1, r2 = fq / f"{srr}_R1.fastq.gz", fq / f"{srr}_R2.fastq.gz"
        if not r2.exists():
            print(f"{srr}: no FASTQ yet, skipping", file=sys.stderr)
            continue
        cdna = read_seqs(r2, args.n_reads)
        bc = read_seqs(r1, min(args.n_reads, 200_000)) if r1.exists() else []

        row = {"sample": srr_sample.get(srr, "?"), "srr": srr}
        row.update(scan(cdna, oligos, min_match))
        row["bc_len"] = modal_len(bc)
        chem, _ = classify_chemistry(row["bc_len"], row.get("cdna_len"))
        row["chem_by_len"] = chem or "?"

        # Barcode identity beats barcode length: this is what proves the whitelist.
        wl_chem, wl_frac = (None, None)
        if bc:
            bcs = Counter(s[:16] for s in bc[:20000] if len(s) >= 16)
            wl_chem, wl_frac = match_whitelists(bcs, args.whitelist_dir, min_frac=0.0)
        row["wl_chem"] = wl_chem or "?"
        row["wl_frac"] = None if wl_frac is None else round(wl_frac, 4)

        row["status"], row["flags"], row["notes"] = judge(row, thr)
        out_rows.append(row)

    out = Path(args.output) if args.output else wd / "qc" / "preflight.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    hdr = f"{'srr':<13}{'bc':>4}{'cDNA':>6}{'wl':>6}{'wl_frac':>9}{'TSO%':>8}{'pred_short%':>13}  status"
    print(f"preflight: {len(out_rows)} run(s)   "
          f"outFilterMatchNminOverLread={min_match}  -> {out}\n")
    print(hdr)
    print("-" * len(hdr))
    for r in out_rows:
        print(f"{r['srr']:<13}{r['bc_len']:>4}{r['cdna_len']:>6}{r['wl_chem']:>6}"
              f"{r['wl_frac']:>9}{r['tso_pct']:>8}{r['pred_too_short_pct']:>13}  {r['status']}")
        if r["flags"]:
            print(f"{'':<13}  └─ {r['flags']}")
        if r["notes"]:
            print(f"{'':<13}  ·  {r['notes']}")
    n_flag = sum(1 for r in out_rows if r["status"] == "FLAG")
    print(f"\n{len(out_rows) - n_flag}/{len(out_rows)} clean")

    # The trade-off curve, for the flagged runs only. Predicted from read lengths alone:
    # how much of the library each candidate filter would discard. Choosing is a human
    # call -- the reads bound the loss, but only an alignment shows whether the reads a
    # looser setting admits are real signal (compare average mapped length) or noise.
    flagged = [r for r in out_rows if r["status"] == "FLAG" and r.get("_sweep")]
    for r in flagged:
        print(f"\nfilter sweep — {r['srr']} (read length {r['cdna_len']} bp, "
              f"currently {min_match:.2f})")
        print(f"  {'--outFilterMatchNminOverLread':<32}{'min match':>11}{'lost %':>9}")
        for t, lost in r["_sweep"]:
            mark = "  <- current" if abs(t - min_match) < 1e-9 else ""
            print(f"  {t:<32.2f}{math.ceil(t * r['cdna_len']):>11}{lost:>9.2f}{mark}")
        print(f"  floor: {r['irrecoverable_pct']}% of reads have no alignable sequence "
              "at all (poly-A / low complexity) — no setting recovers those.")
        print("  NOTE these are LOWER bounds. 'lost %' assumes a read matches its whole "
              "alignable span;\n       real alignments match less (mismatches, indels, "
              "splicing), so measured loss runs\n       higher — on GSE310923, 0.50 "
              "predicted 1.93% and measured 3.87%. Trust the SHAPE\n       (where the "
              "cliff falls), not the absolute values at permissive settings.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
