#!/usr/bin/env python3
"""Offline test for read_qc.py — synthetic FastQC zips, no network, no FastQC.

Builds tiny *_fastqc.zip files by hand (the exact summary.txt / fastqc_data.txt
layout read_qc.py parses) and asserts the ROLE-AWARE verdict:

  * a good cDNA read (R2) + a barcode read (R1) whose content/GC modules FAIL
    *by design* -> sample PASSes (R1 content FAILs must be ignored);
  * a cDNA read with collapsed base quality -> sample FAILs (low_cdna_quality);
  * the pass list contains only the fully-passing sample.
"""
import csv
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

# Repo root = parent of this tests/ dir, so the suite runs wherever the repo lives.
PIPE = Path(__file__).resolve().parents[1]
READ_QC = PIPE / "workflow" / "scripts" / "read_qc.py"


def fastqc_data(n_reads, length, mean_q, adapter_pct, dedup_pct, content_fail):
    """Render a minimal FastQC summary.txt + fastqc_data.txt pair."""
    verdict = "FAIL" if content_fail else "PASS"
    summary = (
        "PASS\tBasic Statistics\tx.fastq.gz\n"
        "PASS\tPer base sequence quality\tx.fastq.gz\n"
        "PASS\tPer sequence quality scores\tx.fastq.gz\n"
        f"{verdict}\tPer base sequence content\tx.fastq.gz\n"
        f"{verdict}\tPer sequence GC content\tx.fastq.gz\n"
        f"{verdict}\tOverrepresented sequences\tx.fastq.gz\n"
        "PASS\tAdapter Content\tx.fastq.gz\n"
    )
    data = [
        "##FastQC\t0.11.9",
        ">>Basic Statistics\tpass",
        "#Measure\tValue",
        "Filename\tx.fastq.gz",
        f"Total Sequences\t{n_reads}",
        f"Sequence length\t{length}",
        "%GC\t48",
        "Sequences flagged as poor quality\t0",
        ">>END_MODULE",
        ">>Per base sequence quality\tpass",
        "#Base\tMean\tMedian",
        f"1\t{mean_q}\t{mean_q}",
        f"2\t{mean_q}\t{mean_q}",
        ">>END_MODULE",
        ">>Adapter Content\tpass",
        "#Position\tIllumina Universal Adapter",
        f"1\t{adapter_pct}",
        f"2\t{adapter_pct}",
        ">>END_MODULE",
        ">>Sequence Duplication Levels\tpass",
        f"#Total Deduplicated Percentage\t{dedup_pct}",
        ">>END_MODULE",
    ]
    return summary, "\n".join(data) + "\n"


def write_zip(path, base, summary, data):
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{base}/summary.txt", summary)
        zf.writestr(f"{base}/fastqc_data.txt", data)


def main():
    tmp = Path(tempfile.mkdtemp(prefix="readqc_fx_"))
    fq = tmp / "fastqc"
    fq.mkdir()

    # samples.tsv: SRR_good -> s_good, SRR_bad -> s_bad
    st = tmp / "samples.tsv"
    st.write_text("sample\tsrr\ns_good\tSRR_good\ns_bad\tSRR_bad\n")

    # --- s_good: barcode R1 (28bp, q30, content modules FAIL by design),
    #             cDNA R2 (91bp, q35, 2% adapter). Expect PASS. ------------- #
    s, d = fastqc_data(30_000_000, 28, 30.0, 2.0, 65, content_fail=True)
    write_zip(fq / "SRR_good_R1_fastqc.zip", "SRR_good_R1_fastqc", s, d)
    s, d = fastqc_data(30_000_000, 91, 35.0, 2.0, 40, content_fail=False)
    write_zip(fq / "SRR_good_R2_fastqc.zip", "SRR_good_R2_fastqc", s, d)

    # --- s_bad: barcode R1 fine; cDNA R2 base quality collapsed (q18 < 28).
    #            Expect FAIL (low_cdna_quality). ------------------------------ #
    s, d = fastqc_data(30_000_000, 28, 30.0, 2.0, 65, content_fail=True)
    write_zip(fq / "SRR_bad_R1_fastqc.zip", "SRR_bad_R1_fastqc", s, d)
    s, d = fastqc_data(30_000_000, 91, 18.0, 2.0, 40, content_fail=False)
    write_zip(fq / "SRR_bad_R2_fastqc.zip", "SRR_bad_R2_fastqc", s, d)

    report = tmp / "read_qc.tsv"
    passlist = tmp / "read_qc_pass.txt"
    cmd = [sys.executable, str(READ_QC),
           "--fastqc-dir", str(fq), "--samples", str(st),
           "--report", str(report), "--passlist", str(passlist),
           "--min-reads", "1000", "--min-mean-quality", "28",
           "--max-adapter-fraction", "0.10", "--min-barcode-len", "24"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    print(proc.stdout)
    print(proc.stderr, file=sys.stderr)

    rows = {r["sample"]: r for r in csv.DictReader(open(report), delimiter="\t")}
    passing = [x for x in passlist.read_text().split("\n") if x]

    ok = True

    def check(cond, msg):
        nonlocal ok
        ok = ok and cond
        print(f"  [{'OK' if cond else 'FAIL'}] {msg}")

    check(proc.returncode == 0, f"read_qc.py exits 0 (report layer); got {proc.returncode}")
    check(rows.get("s_good", {}).get("status") == "PASS",
          f"s_good PASS despite R1 content FAILs; got {rows.get('s_good', {}).get('status')}")
    check("low_cdna_quality" in rows.get("s_bad", {}).get("reasons", ""),
          f"s_bad FAILs on cDNA quality; reasons={rows.get('s_bad', {}).get('reasons')!r}")
    check(rows.get("s_bad", {}).get("status") == "FAIL",
          f"s_bad status FAIL; got {rows.get('s_bad', {}).get('status')}")
    check(passing == ["s_good"], f"pass list is only s_good; got {passing}")

    print(f"\nread_qc.tsv:\n{report.read_text()}")
    if ok:
        print("\nALL READ-QC ASSERTIONS PASSED")
        sys.exit(0)
    print("\nREAD-QC TEST FAILED")
    sys.exit(1)


if __name__ == "__main__":
    main()
