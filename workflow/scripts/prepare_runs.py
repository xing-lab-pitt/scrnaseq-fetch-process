#!/usr/bin/env python3
"""Resolve a GEO (GSE) or SRA study (SRP) accession into a Snakemake sample sheet.

Reuses the project's ncbi_utils for metadata + ENA URL resolution, then writes a
TSV with one row per sequencing run:

    sample  srr  source  fastq_1_url  fastq_2_url  fastq_1_md5  fastq_2_md5  chemistry  strand

`source` is "ena" when ENA mirrors pre-split FASTQs (downloaded via curl) or
"sra" when the run is only in SRA (downloaded via prefetch + fasterq-dump).
Every run known to NCBI appears, so runs absent from ENA are NOT dropped.

Convention: fastq_1 = the barcode/UMI read (10x R1), fastq_2 = the cDNA read (R2).
STARsolo aligns them in the standard two-file mode (cDNA read first, barcode read
second). Runs are grouped into samples by GSM so multiple runs/lanes of one sample
are merged at alignment time.

IMPORTANT: SRA/ENA sometimes exports three files (I1 index + R1 + R2) and the
numeric suffix does NOT always map to R1/R2. Always verify read lengths of the
chosen files (the barcode read should be ~26-28 bp). See select_reads().
"""
import argparse
import csv
import sys
from pathlib import Path


def select_reads(urls, expected_bc_len=28):
    """Pick (barcode_url, cdna_url) from an ENA fastq_ftp URL list.

    2 files  -> (_1, _2): _1 is R1/barcode, _2 is R2/cDNA.
    3 files  -> (_2, _3): _1 is the I1 index; _2 is R1/barcode, _3 is R2/cDNA.
    Returns (barcode_url, cdna_url). Raises on anything unexpected.
    """
    ordered = sorted(urls, key=lambda u: u.rsplit("_", 1)[-1])
    if len(ordered) == 2:
        return ordered[0], ordered[1]
    if len(ordered) == 3:
        return ordered[1], ordered[2]
    # 1 file (ENA dropped the barcode read as "technical") or >3: ENA can't give
    # us a usable barcode+cDNA pair. Signal the caller to use the SRA fallback
    # instead of crashing the whole sheet.
    return None, None


def peek_read_length(url, n_reads=2000, n_bytes=1_000_000):
    """Return the most common read length in the first ~n_bytes of a gzipped
    FASTQ, via an HTTP range request (no full download). None on failure."""
    import gzip, io
    from collections import Counter
    from urllib.request import Request, urlopen
    try:
        req = Request(url, headers={"Range": f"bytes=0-{n_bytes}"})
        with urlopen(req, timeout=60) as r:
            blob = r.read()
        # gzip may complain about the truncated tail; read what we can.
        seqs, gz = [], gzip.GzipFile(fileobj=io.BytesIO(blob))
        for i, line in enumerate(gz):
            if i % 4 == 1:
                seqs.append(len(line.rstrip(b"\n")))
            if len(seqs) >= n_reads:
                break
    except (EOFError, OSError):
        pass
    except Exception:
        return None
    if not seqs:
        return None
    return Counter(seqs).most_common(1)[0][0]


# v3 and v4 (GEM-X) share identical barcode geometry: 16 bp CB + 12 bp UMI = 28 bp
# barcode read. ONLY the whitelist differs (v3: 3M-february-2018; v4/GEM-X:
# 3M-3pgex-may-2023), so read length alone CANNOT distinguish them — the classifier
# reports "v3" for 28 bp and the guard treats the two as geometry-equivalent.
SAME_GEOMETRY = {"v3", "v4"}


def classify_chemistry(bc_len, cdna_len=None):
    """Map a barcode read length to a 10x chemistry, or flag non-droplet data.

    Returns (chem, message) where chem is 'v3' | 'v2' | None.
      * 28 bp -> v3 (also v4/GEM-X — same 16 CB + 12 UMI geometry; whitelist differs)
      * 26 bp -> v2 (16 CB + 10 UMI)
      * a barcode read as long as the cDNA read (symmetric mates) -> not 10x
        droplet (bulk / full-length / Smart-seq) -> chem None, refusal message.
    """
    if bc_len in (28, 27):        # some v3 exports trim by 1
        return "v3", (f"barcode read {bc_len} bp -> 10x 3' v3 or v4/GEM-X "
                      "(16 CB + 12 UMI; whitelist distinguishes them)")
    if bc_len in (26, 24):        # v2: 16 CB + 10 UMI (some trimmed)
        return "v2", f"barcode read {bc_len} bp -> 10x 3' v2 (16 CB + 10 UMI)"
    if cdna_len is not None and bc_len >= cdna_len - 5:
        return None, (f"both reads ~{bc_len} bp (symmetric mates): this is NOT "
                      "10x droplet data (looks like bulk / full-length / "
                      "Smart-seq). STARsolo has no barcode read to parse here.")
    return None, (f"barcode read {bc_len} bp: not a recognized 10x length "
                  "(expected 26 for v2 or 28 for v3).")


def fetch_ena_table(study_accession):
    """ENA filereport incl. md5. Returns {srr: {'urls': [...], 'md5': [...]}}."""
    from urllib.request import urlopen
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={study_accession}"
           "&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv")
    with urlopen(url, timeout=60) as r:
        lines = r.read().decode().strip().split("\n")
    out = {}
    header = lines[0].split("\t")
    ri, fi, mi = (header.index(c) for c in ("run_accession", "fastq_ftp", "fastq_md5"))
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) <= max(ri, fi, mi) or not f[fi]:
            continue
        out[f[ri]] = {
            "urls": [f"https://{u}" for u in f[fi].split(";") if u],
            "md5": f[mi].split(";") if len(f) > mi else [],
        }
    return out


def resolve_ena_study(accession, runs_meta):
    """Find the ENA-queryable study accession for any input accession.

    For GSE we resolve GSE->SRP; for SRX/SRR input the study isn't the input
    itself, so read it from the run metadata (sra_study). Falls back to the
    accession as-is (works for SRP/PRJNA input).
    """
    if accession.upper().startswith("GSE"):
        from ncbi_utils import fetch_sra_study_accession
        return fetch_sra_study_accession(accession) or accession
    for run in runs_meta or []:
        if run.get("sra_study"):
            return run["sra_study"]
    return accession


def build_rows(accession, srr_to_gsm, runs_meta=None):
    """Build sample-sheet rows for EVERY run, preferring ENA URLs.

    Each row gets a `source` column:
      * "ena" — ENA mirrors a usable barcode+cDNA FASTQ pair; download via curl.
      * "sra" — run is not in ENA, or ENA only has a single (cDNA) file because
                the barcode read was submitted as "technical"; download from SRA
                via prefetch + fasterq-dump --include-technical (Snakefile branch).

    Each row also gets a `chemistry` column (v2/v3) based on barcode read length,
    or empty string if undetectable (SRA-only; verified after download).

    Each row gets a `strand` column, defaulting to "Forward" (10x 3' GEX).
    This is NOT auto-detected — barcode geometry cannot distinguish 3' from 5' —
    so set it to "Reverse" by hand for 10x 5' GEX datasets (see the Snakefile /
    config note on --soloStrand).

    The full SRR list comes from NCBI runinfo (srr_to_gsm), so runs missing
    from ENA still appear in the sheet instead of being silently dropped.
    """
    study = resolve_ena_study(accession, runs_meta)
    ena = fetch_ena_table(study)

    # Scope the run set to the input accession.
    #   Study-level input (GSE / SRP / PRJNA): take every run, unioning NCBI's
    #     list with ENA's keys so a run missing from either source isn't dropped.
    #   Sub-study input (SRX experiment / SRR run / GSM sample): the ENA
    #     filereport is study-wide, so union would pull in unrelated sibling
    #     runs. Restrict to exactly the runs NCBI resolved for this accession.
    a = accession.upper()
    is_substudy = a.startswith(("SRX", "SRR", "ERX", "ERR", "GSM"))
    if a.startswith(("SRR", "ERR")):
        # An explicit run accession: NCBI resolves it via its parent experiment
        # UID and returns all sibling runs, so keep only the one asked for.
        all_srr = [s for s in srr_to_gsm if s.upper() == a]
    elif is_substudy:
        all_srr = sorted(srr_to_gsm)
    else:
        all_srr = sorted(set(srr_to_gsm) | set(ena))

    rows = []
    for srr in all_srr:
        info = ena.get(srr)
        bc_url = cdna_url = None
        detected_chem = ""  # empty for SRA (detected after download)

        if info and info.get("urls"):
            bc_url, cdna_url = select_reads(info["urls"])

        if bc_url and cdna_url:
            # Peek BOTH files to find which is barcode (26-28bp) vs cDNA (50-100bp)
            len1 = peek_read_length(bc_url)
            len2 = peek_read_length(cdna_url)

            # Find the barcode read (26-28bp for 10x)
            bc_len = None
            if len1 and len1 in (24, 26, 27, 28):
                bc_len = len1
            elif len2 and len2 in (24, 26, 27, 28):
                bc_len = len2
                # Swap URLs if _2 is actually the barcode
                bc_url, cdna_url = cdna_url, bc_url

            if bc_len:
                chem, _ = classify_chemistry(bc_len, cdna_len=len2 if bc_len == len1 else len1)
                detected_chem = chem if chem else ""

            by_suffix = dict(zip(info["urls"], info["md5"])) if info["md5"] else {}
            rows.append({
                "sample": srr_to_gsm.get(srr, srr),
                "srr": srr,
                "source": "ena",
                "fastq_1_url": bc_url,
                "fastq_2_url": cdna_url,
                "fastq_1_md5": by_suffix.get(bc_url, ""),
                "fastq_2_md5": by_suffix.get(cdna_url, ""),
                "chemistry": detected_chem,
                "strand": "Forward",   # 10x 3' default; set Reverse for 10x 5' GEX
            })
        else:
            # Not in ENA, or ENA has no usable barcode+cDNA pair (single file)
            # -> SRA fallback. No URLs/md5; chemistry detected after download.
            rows.append({
                "sample": srr_to_gsm.get(srr, srr),
                "srr": srr,
                "source": "sra",
                "fastq_1_url": "",
                "fastq_2_url": "",
                "fastq_1_md5": "",
                "fastq_2_md5": "",
                "chemistry": "",  # detected after download
                "strand": "Forward",   # 10x 3' default; set Reverse for 10x 5' GEX
            })
    return rows


def check_chemistry(rows, expected="v3"):
    """Peek the barcode read of the first ENA run and warn/refuse on a mismatch.

    Only ENA rows can be peeked cheaply (range request on a public URL); SRA-only
    runs are skipped with a note (their barcode read is verified after download).
    Returns the detected chemistry ('v2'/'v3') or None. Raises SystemExit if the
    data is clearly not 10x droplet (symmetric mates).
    """
    ena_rows = [r for r in rows if r["source"] == "ena" and r["fastq_1_url"]]
    if not ena_rows:
        print("  chemistry: no ENA rows to peek (SRA-only); verify barcode read "
              "length (~28 bp v3 / ~26 bp v2) after the first download.")
        return None
    r = ena_rows[0]
    bc = peek_read_length(r["fastq_1_url"])
    cdna = peek_read_length(r["fastq_2_url"]) if r["fastq_2_url"] else None
    if bc is None:
        print("  chemistry: could not peek read length (network?); skipping guard.")
        return None
    chem, msg = classify_chemistry(bc, cdna)
    print(f"  chemistry check ({r['srr']}): {msg}")
    if chem is None:
        sys.exit(
            "REFUSING to write sheet: this pipeline is 10x droplet STARsolo only.\n"
            f"  {msg}\n"
            "  For bulk / full-length / Smart-seq data use a different pipeline "
            "(plain STAR alignReads + featureCounts), not this one.")
    # v3 and v4 look identical here (both 28 bp) — the length classifier always
    # says "v3". Don't warn when the user declared v4; just remind them the
    # whitelist is what makes it v4, since we can't verify that from read length.
    if {chem, expected} <= SAME_GEOMETRY and chem != expected:
        print(f"  NOTE: read length can't tell v3 from v4 (both 28 bp). You declared "
              f"{expected}; make sure config.yaml points at the {expected} whitelist "
              "(v4/GEM-X: 3M-3pgex-may-2023.txt; v3: 3M-february-2018.txt).")
    elif chem != expected:
        print(f"  WARNING: detected {chem} but config.yaml is set for {expected}. "
              f"Update config chemistry (v2: umi_len 10, 737K-august-2016 "
              "whitelist) before running, or matrices will be near-empty.")
    return chem


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("accession", help="GSE or SRP accession")
    ap.add_argument("-o", "--output", default="config/samples.tsv")
    ap.add_argument("--tools-dir", default=None,
                    help="Optional dir to prepend to sys.path for ncbi_utils. "
                         "Not needed when ncbi_utils.py sits next to this script.")
    ap.add_argument("--chem", default="v3", choices=["v2", "v3", "v4"],
                    help="Chemistry the config is set for; guard warns on mismatch. "
                         "v3 and v4 (GEM-X) share geometry — v4 only changes the "
                         "whitelist (3M-3pgex-may-2023.txt).")
    ap.add_argument("--no-chem-check", action="store_true",
                    help="Skip the pre-write chemistry guard (peek + refuse).")
    args = ap.parse_args()

    # ncbi_utils.py is vendored alongside this script, so its directory is
    # already importable. --tools-dir is an optional override for an external copy.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    if args.tools_dir:
        sys.path.insert(0, args.tools_dir)

    # Enumerate every run for the accession (GSE / GSM / SRP / SRX / SRR / PRJNA).
    # The detailed (runinfo CSV) resolver reliably fills GSM even on the
    # BioProject fallback path, and fetch_sra_run_info_detailed now searches SRA
    # accessions directly (plain term) instead of only GSE ([GEO]).
    from ncbi_utils import fetch_sra_run_info_detailed
    runs_meta = fetch_sra_run_info_detailed(args.accession)
    srr_to_gsm = {}
    for run in runs_meta:
        if run.get("srr"):
            srr_to_gsm[run["srr"]] = run.get("gsm") or run["srr"]

    rows = build_rows(args.accession, srr_to_gsm, runs_meta)
    if not rows:
        sys.exit(f"No runs resolved for {args.accession}")

    # Chemistry guard: peek the barcode read and refuse non-10x data before we
    # write a sheet the STARsolo pipeline can't process. Skip with --no-chem-check.
    if not args.no_chem_check:
        check_chemistry(rows, expected=args.chem)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_ena = sum(1 for r in rows if r["source"] == "ena")
    n_sra = sum(1 for r in rows if r["source"] == "sra")
    print(f"Wrote {out}: {len(rows)} runs across "
          f"{len({r['sample'] for r in rows})} samples "
          f"({n_ena} via ENA, {n_sra} via SRA fallback)")
    if n_sra:
        print(f"  NOTE: {n_sra} run(s) not mirrored on ENA -> fetched from SRA "
              "(prefetch + fasterq-dump). Requires sra_tools_bin on PATH.")


if __name__ == "__main__":
    main()
