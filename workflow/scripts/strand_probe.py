#!/usr/bin/env python3
"""Decide --soloStrand empirically, on a read subsample, before a full alignment.

10x 3' GEX needs Forward; 5' GEX needs Reverse. Nothing in the reads distinguishes them:
barcode geometry and read lengths are identical, so preflight.py cannot help and the
sample sheet's `strand` column is only ever as good as the GEO record it was copied from.
Get it wrong and gene counts collapse (reads-mapped-to-gene well under 0.1) — today that
only surfaces after a multi-hour alignment.

So run the experiment: align the SAME subsample twice, changing ONLY --soloStrand, and
compare reads-mapped-to-gene. The correct setting wins by a wide margin; a narrow margin
means something else is wrong and neither answer should be trusted.

Read-only with respect to the pipeline: writes solely into its own tmp/output dirs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

STRANDS = ("Forward", "Reverse")
# Below this ratio the two strands are too close to call: the library is not behaving
# like a stranded 3'/5' protocol and the real problem is elsewhere (reference, chemistry).
DECISIVE_RATIO = 3.0


def subsample(src, dst, n_reads):
    with gzip.open(src, "rt") as fh, open(dst, "w") as out:
        for i, line in enumerate(fh):
            if i >= n_reads * 4:
                break
            out.write(line)


def gene_metric(solo_dir, feature):
    """'Reads Mapped to Gene: Unique+Multiple Gene' from a STARsolo Summary.csv."""
    summary = Path(solo_dir) / feature / "Summary.csv"
    if not summary.exists():
        return None
    for row in csv.reader(open(summary)):
        if row and row[0].startswith("Reads Mapped to Gene: Unique+Multiple"):
            return float(row[1])
    return None


def run_star(cdna, bc, strand, cfg, outdir):
    chem = cfg["chemistry"]
    cmd = [
        "STAR", "--runMode", "alignReads",
        "--genomeDir", cfg["index"],
        "--readFilesIn", str(cdna), str(bc),
        "--soloType", "CB_UMI_Simple",
        "--soloCBwhitelist", chem["whitelist"],
        "--soloCBstart", str(chem.get("cb_start", 1)),
        "--soloCBlen", str(chem.get("cb_len", 16)),
        "--soloUMIstart", str(chem.get("umi_start", 17)),
        "--soloUMIlen", str(chem.get("umi_len", 12)),
        "--soloStrand", strand,
        "--soloFeatures", cfg["feature"],
        "--soloCBmatchWLtype", "1MM_multi",
        "--outSAMtype", "None",
        "--runThreadN", str(cfg["threads"]),
        "--outFileNamePrefix", f"{outdir}/{strand}_",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return gene_metric(f"{outdir}/{strand}_Solo.out", cfg["feature"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir")
    ap.add_argument("--srr", required=True, help="one run to probe")
    ap.add_argument("--config", required=True)
    ap.add_argument("--n-reads", type=int, default=1_000_000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--tmp", default=None, help="scratch dir (default: system tmp)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    import yaml
    raw = yaml.safe_load(open(args.config)) or {}
    wd = Path(args.workdir)
    cfg = {
        "index": str(wd / "star_index"),
        "chemistry": raw.get("chemistry", {}),
        "feature": raw.get("feature", "Gene"),
        "threads": args.threads,
    }

    r1 = wd / "fastq" / f"{args.srr}_R1.fastq.gz"
    r2 = wd / "fastq" / f"{args.srr}_R2.fastq.gz"
    for f in (r1, r2):
        if not f.exists():
            sys.exit(f"missing {f}")

    tmp = tempfile.mkdtemp(dir=args.tmp, prefix="strand_probe_")
    try:
        sub_cdna, sub_bc = f"{tmp}/cdna.fastq", f"{tmp}/bc.fastq"
        subsample(r2, sub_cdna, args.n_reads)
        subsample(r1, sub_bc, args.n_reads)

        res = {}
        for strand in STRANDS:
            res[strand] = run_star(sub_cdna, sub_bc, strand, cfg, tmp)
            print(f"  {strand:8s} reads_mapped_gene = {res[strand]}", file=sys.stderr)

        fwd, rev = res["Forward"], res["Reverse"]
        if fwd is None or rev is None:
            verdict, ratio, best = "ERROR", None, None
        else:
            best = "Forward" if fwd >= rev else "Reverse"
            lo = min(fwd, rev) or 1e-9
            ratio = max(fwd, rev) / lo
            verdict = best if ratio >= DECISIVE_RATIO else "INCONCLUSIVE"

        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["srr", "forward", "reverse", "ratio", "recommended", "verdict"])
            w.writerow([args.srr, fwd, rev,
                        None if ratio is None else round(ratio, 1), best, verdict])

        print(f"\nstrand probe {args.srr}: Forward={fwd}  Reverse={rev}  "
              f"ratio={'?' if ratio is None else round(ratio, 1)}x  -> {verdict}")
        if verdict == "INCONCLUSIVE":
            print("  neither strand dominates: suspect the reference or chemistry, "
                  "not the strand setting.")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
