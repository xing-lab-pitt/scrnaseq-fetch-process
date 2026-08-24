#!/usr/bin/env python3
"""Resolve a GEO (GSE) or SRA study (SRP) accession into a Snakemake sample sheet.

Reuses the project's ncbi_utils for metadata + ENA URL resolution, then writes a
TSV with one row per sequencing run:

    sample  srr  source  fastq_1_url  fastq_2_url  fastq_1_md5  fastq_2_md5  chemistry  strand  barcode_read_length  software

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
import os
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


# chem -> 10x cell-barcode whitelist filename (same names the Snakefile/config use).
WHITELIST_SUFFIX = {
    "v2":  "737K-august-2016.txt",
    "v3":  "3M-february-2018.txt",
    "v4":  "3M-3pgex-may-2023.txt",
    "arc": "737K-arc-v1.txt",
}
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.yaml"


def whitelist_dir_from_config(config_path=DEFAULT_CONFIG):
    """The directory holding the 10x whitelists, taken from `chemistry.whitelist`.

    Config already carries that path for STARsolo, and the Snakefile reads the
    directory off it the same way, so this stays the single place a site records
    where its whitelists live. Returns "" when the config or the key is absent —
    the barcode probe then finds no whitelists and goes quiet, which the caller
    reports.
    """
    try:
        import yaml
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh) or {}
    except Exception:
        return ""
    wl = (cfg.get("chemistry") or {}).get("whitelist") or ""
    return str(Path(wl).parent) if wl else ""


def resolve_whitelist_dir(cli_value="", config_path=DEFAULT_CONFIG):
    """--whitelist-dir wins, then $SCRNASEQ_WHITELIST_DIR, then config."""
    return (cli_value
            or os.environ.get("SCRNASEQ_WHITELIST_DIR", "")
            or whitelist_dir_from_config(config_path))


def _fetch_read_seqs(url, n_reads=2000, n_bytes=2_000_000):
    """Return up to n_reads sequence lines from the head of a gzipped FASTQ, via an
    HTTP range request (no full download). Empty list on failure."""
    import gzip, io
    from urllib.request import Request, urlopen
    seqs = []
    try:
        req = Request(url, headers={"Range": f"bytes=0-{n_bytes}"})
        with urlopen(req, timeout=60) as r:
            blob = r.read()
        gz = gzip.GzipFile(fileobj=io.BytesIO(blob))
        for i, line in enumerate(gz):
            if i % 4 == 1:
                seqs.append(line.rstrip(b"\n").decode("ascii", "replace"))
            if len(seqs) >= n_reads:
                break
    except (EOFError, OSError):
        pass
    except Exception:
        return []
    return seqs


def probe_barcode_chem(read_url, whitelist_dir, cb_len=16, n_reads=2000, min_frac=0.5):
    """Identify 10x chemistry by matching a read's leading barcode against whitelists.

    For reads too long to classify by length (an over-sequenced R1, e.g. 150 bp on a
    10x run where the sequencer read through the 28 bp barcode into the poly-T tail),
    read length alone cannot tell 10x-droplet from Smart-seq/bulk. The decisive test
    is whether the read's first `cb_len` bp ARE cell barcodes: sample reads, take the
    leading cb_len bp of each, and match against every candidate whitelist. A real 10x
    read matches its whitelist at a high fraction (~0.8); random 16-mers match a 737K–
    3M entry list at ~1e-4. A hit therefore both proves droplet data AND names the
    chemistry (whichever whitelist the barcodes came from).

    Streams each whitelist file once (O(sample) memory), so the multi-million-line
    v3/v4 lists are never fully loaded. Returns (chem, frac) for the best match, or
    (None, best_frac) if no whitelist matches at >= min_frac.
    """
    from collections import Counter
    seqs = _fetch_read_seqs(read_url, n_reads=n_reads)
    bcs = Counter(s[:cb_len] for s in seqs if len(s) >= cb_len)
    total = sum(bcs.values())
    if not total:
        return None, 0.0
    best_chem, best_frac = None, 0.0
    for chem, suffix in WHITELIST_SUFFIX.items():
        wl = Path(whitelist_dir) / suffix
        if not wl.exists():
            continue
        hits = 0
        with open(wl) as fh:
            for line in fh:
                bc = line.strip()
                if bc in bcs:
                    hits += bcs[bc]
        frac = hits / total
        if frac > best_frac:
            best_chem, best_frac = chem, frac
    if best_frac >= min_frac:
        return best_chem, best_frac
    return None, best_frac


# GEO serves every sample in machine-readable SOFT text; `!Sample_data_processing`
# is the submitter's own free-text statement of how the data was made. targ=self
# keeps the payload to just this GSM (no platform/series dump).
GEO_SOFT_URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
                "acc={acc}&targ=self&form=text&view=quick")

# Ordered (regex, software, chem_hint, is_droplet); FIRST match wins, so the most
# specific / disambiguating patterns come first: cellranger-arc before plain
# cellranger, alevin before bulk salmon, bustools before bulk kallisto. `chem_hint`
# is set only when the software pins a chemistry read length can't (arc is 28 bp,
# identical to v3/v4 by length — only cellranger-ARC in the record tells them apart).
# `is_droplet` False marks plate-based / bulk software so the guard can refuse with
# a precise reason instead of the generic symmetric-mates message.
_SOFTWARE_RULES = [
    (r"cell\s*-?\s*ranger[\s-]*arc",                    "cellranger-arc", "arc", True),
    (r"cell\s*-?\s*ranger",                             "cellranger",     "",    True),
    (r"star\s*-?\s*solo",                               "starsolo",       "",    True),
    (r"kallisto[\s|]*bustools|\bkb[-\s]?python\b|\bbustools\b", "kb-python", "", True),
    (r"salmon[\s|]*alevin|alevin[-\s]?fry|\balevin\b",  "alevin",         "",    True),
    (r"drop-?seq|dropseqtools",                         "dropseq",        "",    True),
    (r"smart-?seq",                                     "smartseq",       "",    False),
    (r"\brsem\b",                                       "rsem",           "",    False),
    (r"feature\s*counts|\bsubread\b",                   "featurecounts",  "",    False),
    (r"\bhisat2?\b",                                    "hisat2",         "",    False),
    (r"\btophat\b",                                     "tophat",         "",    False),
    (r"\bkallisto\b",                                   "kallisto",       "",    False),
    (r"\bsalmon\b",                                     "salmon",         "",    False),
]
# software token -> is_droplet, so the guard can re-derive the verdict from the
# `software` column alone (build_rows already resolved it).
_SOFTWARE_DROPLET = {sw: drop for _, sw, _, drop in _SOFTWARE_RULES}


def software_is_droplet(software):
    """True / False / None (unknown) for a software token from `_SOFTWARE_RULES`."""
    return _SOFTWARE_DROPLET.get(software, None)


def fetch_geo_software(gsm):
    """Read a GSM's GEO record and classify the processing software the submitter
    named in the free-text `!Sample_data_processing` field.

    Returns {"software", "chem_hint", "is_droplet", "text"}: a recognized aligner
    (is_droplet True/False), an empty software with is_droplet None if the field is
    present but names nothing we know, or None if the record can't be fetched or has
    no data_processing field at all.

    This is a metadata ASSERT (see the detect-vs-assert principle): the submitter's
    own statement of how the data was made, reachable even for SRA-only runs the
    FASTQ range-probe can't peek. cellranger-arc pins `arc` (read length can't tell
    it from v3/v4); a Smart-seq/RSEM/HISAT2 mention confirms NON-droplet.
    """
    import re
    from urllib.request import urlopen
    if not gsm or not gsm.upper().startswith("GSM"):
        return None
    try:
        with urlopen(GEO_SOFT_URL.format(acc=gsm), timeout=60) as r:
            txt = r.read().decode("utf-8", "replace")
    except Exception:
        return None
    dp = " ".join(re.findall(r"^!Sample_data_processing\s*=\s*(.*)$", txt, re.M))
    if not dp.strip():
        return None
    low = dp.lower()
    for rx, sw, chem, drop in _SOFTWARE_RULES:
        if re.search(rx, low):
            return {"software": sw, "chem_hint": chem, "is_droplet": drop, "text": dp}
    return {"software": "", "chem_hint": "", "is_droplet": None, "text": dp}


# v3, v4 (GEM-X) and arc (Multiome) all share identical barcode geometry: 16 bp CB
# + 12 bp UMI = 28 bp barcode read. ONLY the whitelist differs (v3: 3M-february-2018;
# v4/GEM-X: 3M-3pgex-may-2023; arc: 737K-arc-v1), so read length alone CANNOT
# distinguish them — the classifier reports "v3" for 28 bp and the guard treats all
# three as geometry-equivalent. arc is never auto-detected from length: pass --multiome.
SAME_GEOMETRY = {"v3", "v4", "arc"}


def classify_chemistry(bc_len, cdna_len=None):
    """Map a barcode read length to a 10x chemistry, or flag non-droplet data.

    Returns (chem, message) where chem is 'v3' | 'v2' | None.
      * 28 bp -> v3 (also v4/GEM-X — same 16 CB + 12 UMI geometry; whitelist differs)
      * 26 bp -> v2 (16 CB + 10 UMI)
      * a barcode read as long as the cDNA read (symmetric mates) -> not 10x
        droplet (bulk / full-length / Smart-seq) -> chem None, refusal message.
    """
    if bc_len in (28, 27):        # some v3 exports trim by 1
        return "v3", (f"barcode read {bc_len} bp -> 10x 3' v3, v4/GEM-X, or arc/Multiome "
                      "(16 CB + 12 UMI; whitelist distinguishes them)")
    if bc_len in (26, 24):        # v2: 16 CB + 10 UMI (some trimmed)
        return "v2", f"barcode read {bc_len} bp -> 10x 3' v2 (16 CB + 10 UMI)"
    if cdna_len is not None and bc_len >= cdna_len - 5:
        return None, (f"both reads ~{bc_len} bp (symmetric mates): NOT 10x droplet "
                      "data — no short (~26-28 bp) barcode read for STARsolo to parse. "
                      "(bulk, or full-length/plate-based single cell like Smart-seq)")
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


def build_rows(accession, srr_to_gsm, runs_meta=None, chem_override="", whitelist_dir=None,
               use_geo_software=True):
    """Build sample-sheet rows for EVERY run, preferring ENA URLs.

    Each row gets a `source` column:
      * "ena" — ENA mirrors a usable barcode+cDNA FASTQ pair; download via curl.
      * "sra" — run is not in ENA, or ENA only has a single (cDNA) file because
                the barcode read was submitted as "technical"; download from SRA
                via prefetch + fasterq-dump --include-technical (Snakefile branch).

    Each row also gets a `chemistry` column (v2/v3) based on barcode read length,
    or empty string if undetectable (SRA-only; verified after download).

    Each row gets a `barcode_read_length` column: "0" when the barcode read is
    over-sequenced (longer than the 28 bp CB+UMI, e.g. a 150 bp 10x R1 read through
    into poly-T), which the Snakefile turns into `--soloBarcodeReadLength 0` so
    STARsolo skips its "barcode read == 28 bp" assertion (CB/UMI still parsed from
    the fixed first-28 positions). Empty for standard-length runs. When a long-mate
    pair can't be classified by length and chemistry isn't overridden, `whitelist_dir`
    enables probe_barcode_chem to decide 10x (+ which chemistry) vs non-droplet.

    `chem_override` (e.g. "arc" from --multiome) forces the chemistry on EVERY row,
    bypassing the length-based call. Use it for chemistries that are geometrically
    indistinguishable from v3/v4 (arc/Multiome is 28 bp) and so cannot be detected
    from read length — the caller asserts the whole dataset is that chemistry.

    Each row gets a `strand` column, defaulting to "Forward" (10x 3' GEX).
    This is NOT auto-detected — barcode geometry cannot distinguish 3' from 5' —
    so set it to "Reverse" by hand for 10x 5' GEX datasets (see the Snakefile /
    config note on --soloStrand).

    Each row gets a `software` column: the processing software named in the run's
    GEO `!Sample_data_processing` field (cellranger-arc, starsolo, smartseq, ...),
    fetched once per study when `use_geo_software` is set. It's provenance in the
    sheet, but also drives detection: cellranger-arc pins `arc` chemistry (which read
    length can't tell from v3/v4) even for SRA-only runs, and non-droplet software
    (Smart-seq/RSEM/HISAT2) lets the guard refuse with a precise reason.

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

    # The barcode-probe verdict is uniform across a study's runs (same chemistry,
    # same which-mate-is-R1 orientation), so compute it once on the first long-mate
    # row and reuse. long_verdict: None (not yet probed) | {"chem":.., "swap":bool}
    # | {"chem": None} (probed, no whitelist match -> not droplet).
    long_verdict = None

    def _decide_long_reads(url_a, url_b):
        """Probe a long-mate pair to identify 10x chemistry and which mate is R1.
        Returns {"chem", "swap"} on a hit (swap=True means url_b is the barcode read),
        or {"chem": None} on a miss."""
        chem_a, frac_a = probe_barcode_chem(url_a, whitelist_dir)
        if chem_a:
            return {"chem": chem_a, "swap": False, "frac": frac_a}
        chem_b, frac_b = probe_barcode_chem(url_b, whitelist_dir)
        if chem_b:
            return {"chem": chem_b, "swap": True, "frac": frac_b}
        return {"chem": None, "frac": max(frac_a, frac_b)}

    # GEO `!Sample_data_processing` names the processing software (cellranger-arc,
    # STARsolo, Smart-seq, ...). It's study-uniform, so fetch it once from the first
    # resolvable GSM and reuse: it (a) pins arc — which read length can't tell from
    # v3/v4 — so it AUTO-detects Multiome without --multiome, (b) fills chemistry for
    # SRA-only runs the FASTQ probe can't reach, and (c) confirms non-droplet software
    # so the guard refuses precisely. See fetch_geo_software / detect-vs-assert.
    sw_verdict = None
    if use_geo_software:
        for s in all_srr:
            if srr_to_gsm.get(s, "").upper().startswith("GSM"):
                sw_verdict = fetch_geo_software(srr_to_gsm[s])
                if sw_verdict is not None:
                    break
    sw = sw_verdict or {}
    sw_str = sw.get("software", "")
    sw_chem = sw.get("chem_hint", "")

    rows = []
    for srr in all_srr:
        info = ena.get(srr)
        bc_url = cdna_url = None
        # --multiome forces every run to `arc`; otherwise start empty (ENA rows get
        # filled from the peek below; SRA rows stay empty, detected after download).
        detected_chem = chem_override
        bc_read_len = ""   # "0" only when the barcode read is over-sequenced (>28 bp)
        probe_rejected = False   # True if the barcode probe actively ruled out 10x

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
                # Standard-length barcode read: classify by length as before.
                if not chem_override:
                    chem, _ = classify_chemistry(bc_len, cdna_len=len2 if bc_len == len1 else len1)
                    detected_chem = chem if chem else ""
            else:
                # No short barcode read among the mates. Could be an over-sequenced-R1
                # 10x run (CB/UMI still in the leading 28 bp) or genuine non-droplet
                # data. The barcode read, whichever it is, is longer than 28 bp -> tell
                # STARsolo to skip its length assertion.
                bc_read_len = "0"
                if chem_override:
                    # User asserted the chemistry (e.g. --multiome); trust it and keep
                    # _1 as the barcode read (10x R1 convention).
                    pass
                elif whitelist_dir:
                    if long_verdict is None:
                        long_verdict = _decide_long_reads(bc_url, cdna_url)
                    if long_verdict["chem"]:
                        detected_chem = long_verdict["chem"]
                        if long_verdict.get("swap"):
                            bc_url, cdna_url = cdna_url, bc_url
                    else:
                        # Probe found no whitelist match -> not 10x droplet. Leave
                        # chemistry empty so the guard refuses; drop the flag (moot).
                        # Mark it rejected so a GEO "cellranger-arc" label can't
                        # resurrect data the sequence probe actively ruled out.
                        detected_chem = ""
                        bc_read_len = ""
                        probe_rejected = True

            # GEO software can name a chemistry read length can't resolve (arc, v3
            # and v4 are all 28 bp). If the record pins one and our length/probe call
            # was the ambiguous default (or empty), upgrade to it. A confident v2
            # (26 bp) is never overridden — v2 isn't in SAME_GEOMETRY — and a probe
            # that actively rejected the run wins over the metadata label.
            if sw_chem and not probe_rejected and (not detected_chem or detected_chem in SAME_GEOMETRY):
                detected_chem = sw_chem

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
                "barcode_read_length": bc_read_len,
                "software": sw_str,
            })
        else:
            # Not in ENA, or ENA has no usable barcode+cDNA pair (single file)
            # -> SRA fallback. No URLs/md5, and the FASTQ probe can't reach it. Fill
            # chemistry from the GEO data_processing software (chem_override wins if
            # --multiome was passed); else leave empty (verified after download).
            rows.append({
                "sample": srr_to_gsm.get(srr, srr),
                "srr": srr,
                "source": "sra",
                "fastq_1_url": "",
                "fastq_2_url": "",
                "fastq_1_md5": "",
                "fastq_2_md5": "",
                "chemistry": chem_override or sw_chem,
                "strand": "Forward",   # 10x 3' default; set Reverse for 10x 5' GEX
                "barcode_read_length": "",   # unknown for SRA-only; verified after download
                "software": sw_str,
            })

    # Study-uniformity backfill: if any ENA run in this study proved over-sequenced
    # (barcode_read_length="0"), its SRA-only siblings almost certainly share that
    # library prep. Propagate "0" to SRA rows whose chemistry we managed to fill, so
    # STARsolo skips its 28 bp assertion for them too. --soloBarcodeReadLength 0 is a
    # safe superset (harmless on a genuine 28 bp R1: CB/UMI still parse from the
    # leading positions), so this can't break a standard-length SRA run.
    if any(r["source"] == "ena" and r["barcode_read_length"] == "0" for r in rows):
        for r in rows:
            if r["source"] == "sra" and r["chemistry"] and not r["barcode_read_length"]:
                r["barcode_read_length"] = "0"
    return rows


def _refuse_not_droplet(msg, single_cell):
    """Hard-stop: the data isn't 10x droplet, so this STARsolo pipeline can't run it.
    The message splits bulk from full-length/plate-based single cell (Smart-seq)."""
    if single_cell:
        sys.exit(
            "REFUSING to write sheet: this pipeline is 10x droplet STARsolo only.\n"
            f"  {msg}\n"
            "  SRA marks this as SINGLE CELL, and that's not a contradiction: the "
            "symmetric mates say it is full-length / plate-based single cell "
            "(Smart-seq family), NOT 10x droplet. There is no in-read cell barcode "
            "to demultiplex here — each cell is its own library/file.\n"
            "  Process it with a full-length single-cell workflow (per-cell STAR "
            "alignReads + featureCounts), not this droplet pipeline.")
    sys.exit(
        "REFUSING to write sheet: this pipeline is 10x droplet STARsolo only.\n"
        f"  {msg}\n"
        "  For bulk / full-length / Smart-seq data use a different pipeline "
        "(plain STAR alignReads + featureCounts), not this one.")


def check_chemistry(rows, expected="v3", single_cell=False):
    """Warn/refuse on the first ENA run based on the chemistry build_rows detected.

    Only ENA rows are resolvable pre-download; SRA-only runs are skipped with a note
    (their barcode read is verified after download). Returns the detected chemistry
    ('v2'/'v3'/'v4'/'arc') or None. Raises SystemExit if the data is clearly not 10x
    droplet — a symmetric/long-mate pair whose leading bases don't match any 10x
    whitelist (build_rows already ran the barcode probe and left chemistry empty).

    `single_cell` (from SRA LibrarySource == "TRANSCRIPTOMIC SINGLE CELL") sharpens
    the refusal: symmetric mates that are still single cell are full-length /
    plate-based (Smart-seq family), NOT bulk — a distinct, but still non-droplet,
    pipeline. Both SRA's "single cell" label and our "not droplet" verdict are then
    correct; the message says so rather than implying the data is bulk.
    """
    ena_rows = [r for r in rows if r["source"] == "ena" and r["fastq_1_url"]]
    if not ena_rows:
        # No ENA row to peek. Fall back to the GEO data_processing software signal,
        # which reaches SRA-only runs the FASTQ range-probe can't.
        sw_str = rows[0].get("software", "") if rows else ""
        chem0 = rows[0].get("chemistry", "") if rows else ""
        if software_is_droplet(sw_str) is False:
            _refuse_not_droplet(
                f"GEO data_processing names non-droplet software ({sw_str}); this "
                "is not 10x droplet data.", single_cell)
        if chem0:
            print(f"  chemistry ({rows[0]['srr']}): SRA-only, set to 10x {chem0} from "
                  f"GEO data_processing software ({sw_str or 'stated'}).")
            return chem0
        print("  chemistry: no ENA rows to peek (SRA-only) and GEO data_processing "
              "named no known aligner; verify barcode read length (~28 bp v3 / "
              "~26 bp v2) after the first download.")
        return None
    r = ena_rows[0]
    chem = r.get("chemistry") or ""
    over_sequenced = r.get("barcode_read_length") == "0"

    if chem:
        # build_rows already resolved the chemistry (by length, override, or the
        # barcode-whitelist probe for over-sequenced R1). Report and validate it.
        if over_sequenced:
            print(f"  chemistry check ({r['srr']}): 10x {chem} via barcode-whitelist "
                  "probe — R1 is over-sequenced (>28 bp); will set "
                  "--soloBarcodeReadLength 0 so STARsolo reads CB/UMI from the "
                  "leading 28 bp and ignores the tail.")
        else:
            print(f"  chemistry check ({r['srr']}): 10x {chem} (barcode read length).")
    else:
        # No chemistry resolved. Peek the reads to build an explanatory message, then
        # refuse: a short unrecognized length, or a long/symmetric pair the barcode
        # probe already rejected (leading bases aren't 10x whitelist barcodes).
        bc = peek_read_length(r["fastq_1_url"])
        cdna = peek_read_length(r["fastq_2_url"]) if r["fastq_2_url"] else None
        if bc is None:
            print("  chemistry: could not peek read length (network?); skipping guard.")
            return None
        _, msg = classify_chemistry(bc, cdna)
        print(f"  chemistry check ({r['srr']}): {msg}")
        _refuse_not_droplet(msg, single_cell)
    # v3 and v4 look identical here (both 28 bp) — the length classifier always
    # says "v3". Don't warn when the user declared v4; just remind them the
    # whitelist is what makes it v4, since we can't verify that from read length.
    if {chem, expected} <= SAME_GEOMETRY and chem != expected:
        print(f"  NOTE: read length can't tell v3/v4/arc apart (all 28 bp). You declared "
              f"{expected}; make sure config.yaml points at the {expected} whitelist "
              "(v3: 3M-february-2018.txt; v4/GEM-X: 3M-3pgex-may-2023.txt; "
              "arc/Multiome: 737K-arc-v1.txt).")
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
    ap.add_argument("--chem", default="v3", choices=["v2", "v3", "v4", "arc"],
                    help="Chemistry the config is set for; guard warns on mismatch. "
                         "v3, v4 (GEM-X) and arc (Multiome) share geometry — they "
                         "differ only by whitelist (v4: 3M-3pgex-may-2023.txt; "
                         "arc: 737K-arc-v1.txt).")
    ap.add_argument("--multiome", action="store_true",
                    help="Mark EVERY run as 10x Multiome (Chromium ARC) GEX: writes "
                         "chemistry=arc into the sheet (STARsolo then uses the "
                         "737K-arc-v1 whitelist). Multiome GEX is 28 bp, IDENTICAL to "
                         "v3/v4 by read length, so it cannot be auto-detected — this "
                         "flag is how you assert it. Pair with feature: GeneFull in "
                         "config/config.yaml for single-nucleus (snRNA) counting.")
    ap.add_argument("--no-chem-check", action="store_true",
                    help="Skip the pre-write chemistry guard (peek + refuse).")
    ap.add_argument("--no-geo-software", action="store_true",
                    help="Skip the GEO !Sample_data_processing lookup that names the "
                         "processing software (cellranger-arc, STARsolo, Smart-seq, "
                         "...). That lookup auto-detects arc/Multiome (which read "
                         "length can't tell from v3/v4), fills chemistry for SRA-only "
                         "runs, and sharpens non-droplet refusals. Use for offline "
                         "runs or to fall back to length/probe detection only.")
    ap.add_argument("--whitelist-dir", default="",
                    help="Dir of 10x cell-barcode whitelists (737K-arc-v1.txt, "
                         "3M-february-2018.txt, …). Used to identify chemistry AND "
                         "confirm droplet data when a run's R1 is over-sequenced "
                         "(>28 bp) so read length alone can't classify it. Defaults "
                         "to the directory of chemistry.whitelist in --config.")
    ap.add_argument("--config", default=str(DEFAULT_CONFIG),
                    help="Pipeline config the whitelist directory is read from "
                         f"(default: {DEFAULT_CONFIG}).")
    args = ap.parse_args()

    # --multiome forces chemistry=arc on all rows and aligns the guard's expected
    # chemistry, so the geometry NOTE (not a spurious mismatch WARNING) fires.
    chem_override = "arc" if args.multiome else ""
    guard_expected = "arc" if args.multiome else args.chem

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

    whitelist_dir = resolve_whitelist_dir(args.whitelist_dir, args.config)
    if whitelist_dir:
        print(f"  whitelist dir: {whitelist_dir}")
    else:
        print("  WARNING: no whitelist dir (chemistry.whitelist missing from "
              f"{args.config}) — the barcode-sequence probe is disabled, so an "
              "over-sequenced R1 can't be confirmed as 10x.")

    rows = build_rows(args.accession, srr_to_gsm, runs_meta, chem_override=chem_override,
                      whitelist_dir=whitelist_dir,
                      use_geo_software=not args.no_geo_software)
    if not rows:
        sys.exit(f"No runs resolved for {args.accession}")

    # Chemistry guard: peek the barcode read and refuse non-10x data before we
    # write a sheet the STARsolo pipeline can't process. Skip with --no-chem-check.
    # SRA LibrarySource tells us whether a refusal is bulk vs full-length single cell
    # (Smart-seq) so the guard message doesn't wrongly imply single-cell data is bulk.
    if not args.no_chem_check:
        single_cell = any(
            "SINGLE CELL" in (run.get("library_source") or "").upper()
            for run in runs_meta)
        check_chemistry(rows, expected=guard_expected, single_cell=single_cell)

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
    if args.multiome:
        print("  chemistry: forced to `arc` (10x Multiome / Chromium ARC GEX) on all "
              "rows -> STARsolo uses the 737K-arc-v1 whitelist. REMEMBER to set "
              "`feature: GeneFull` in config/config.yaml for single-nucleus (snRNA) "
              "counting (nuclear RNA is mostly unspliced pre-mRNA).")


if __name__ == "__main__":
    main()
