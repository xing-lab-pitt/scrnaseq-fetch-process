---
name: scrnaseq-fetch-process
description: |
  Fetch and process raw 10x Genomics scRNA-seq data from a public accession into
  a per-sample .h5ad, using the STARsolo Snakemake pipeline this skill ships with
  (the repository root, referred to below as $PIPE). Use when:
  (1) Turning a GEO/SRA accession (GSE / SRP / PRJNA / SRX / SRR / GSM) into count
      matrices + spliced/unspliced/ambiguous layers (STARsolo Velocyto feature)
  (2) Building or extending config/samples.tsv from an accession
  (3) Launching the pipeline on SLURM and diagnosing controller hangs
  (4) Deciding ENA vs SRA download, or checking that a dataset is actually 10x
      droplet (not bulk / Smart-seq) before spending a STARsolo run
  This is the UPSTREAM half (accession -> .h5ad); downstream analysis of the
  resulting .h5ad (QC, clustering, UMAP, DE) is a separate scanpy step.
---

# scRNA-seq Fetch & Process (accession → .h5ad)

This skill drives the existing Snakemake pipeline — it does **not** reimplement
it. Snakemake owns the DAG, resume-on-failure, and SLURM right-sizing; this skill
is the runbook for the judgment steps around it (accession resolution, chemistry
guard, launch, hang diagnosis, post-run verification).

**`$PIPE`** = the repository root (the directory containing `workflow/`,
`config/`, `run_slurm.sh` — the parent of this `skill/` folder). All commands
below run from there.

> **THIS MACHINE (xing lab cluster) — local activation note.**
> Actual runs use the established working directory, which already has a
> filled-in `config/config.yaml`, the prebuilt human index, and prior `results/`:
> ```
> $PIPE = /net/capricorn/home/xing/lul176/skills_agent/fetch_process_snakemake
> ```
> The shared clean repo (`.../skills_agent/scrnaseq-fetch-process`) ships
> `/path/to` placeholders for distribution — do NOT run there. Environment for
> `run_slurm.sh` on this cluster:
> ```
> export SCRNASEQ_VENV=/net/capricorn/home/xing/lul176/mskcc/blood_combined/.venv
> export SCRNASEQ_EXTRA_PATH=/opt/FastQC:/net/capricorn/home/xing/soh29/libraries/sratoolkit.3.0.2-ubuntu64/bin
> ```
> STAR 2.7.10a + samtools are already in `/usr/bin`.
>
> Reference + whitelists on this cluster:
> ```
> STAR index (human): /net/capricorn/home/xing/lul176/reference/GRCh38/STAR
> 10x whitelists:     /net/capricorn/home/xing/lul176/reference/10x_whitelists/
>   737K-august-2016.txt   (v2, current config — 26bp R1)
>   3M-february-2018.txt   (v3 — 28bp R1)
>   3M-3pgex-may-2023.txt  (v4/GEM-X — same geometry, swap whitelist only)
>   737K-arc-v1.txt        (10x Multiome/ARC GEX — snRNA; set chemistry: arc)
> refdata bundles:    /net/capricorn/home/xing/lul176/reference/refdata-gex-{GRCh38,GRCm39}-2024-A.tar.gz
> To keep a rebuilt index: set reference.star_index to a writable dir (e.g. a
> sibling STAR_rebuild/); resolve_star_index builds into it once, then reuses it.
> ```
> This note is local-only (this copy of the skill is NOT the repo file); keep
> machine paths out of the repo.

**Environment (site-specific — read the repo's README + `config/config.yaml`):**
- Python env with `snakemake` + the SLURM executor plugin and the pins in
  `requirements-pipeline.txt`. `run_slurm.sh` activates `$SCRNASEQ_VENV` if set;
  otherwise `snakemake` must already be on PATH.
- Binary tools (STAR, samtools, FastQC, MultiQC; sra-tools only for `source=sra`)
  on PATH — via `$SCRNASEQ_EXTRA_PATH`, conda, or modules.
- A STAR index + 10x whitelist configured in `config/config.yaml` (copied from
  `config/config.example.yaml`). `reference.star_index` is only a HINT/path: rule
  `resolve_star_index` checks its `genomeVersion` against the pipeline's STAR and
  either symlinks it into `<workdir>/star_index` (compatible), or rebuilds from
  fasta+gtf (mismatch/absent). On a rebuild, if `star_index` names a writable path
  it builds the index THERE (durable — reused, no rebuild, next run) and symlinks
  `<workdir>/star_index` to it; if `star_index` is `""` it builds into the workdir.
  The aligner only ever reads that pipeline-owned `<workdir>/star_index`, so the
  index and aligner always come from the same STAR.

Never edit outside `$PIPE`. All fixes here stay in `$PIPE` (mainly `profiles/slurm/config.yaml` and `config/`).

## Phase 1 — Build the sample sheet from an accession

`prepare_runs.py` resolves any accession type and writes `config/samples.tsv`
(one row per run: `sample srr source fastq_1_url fastq_2_url fastq_1_md5 fastq_2_md5 chemistry strand`).

```bash
cd "$PIPE"
# activate your Python env first (e.g. source "$SCRNASEQ_VENV/bin/activate")
python workflow/scripts/prepare_runs.py <ACCESSION> -o config/samples.tsv
```

Accession scoping (handled automatically):
- **GSE / SRP / PRJNA** → every run in the study.
- **SRX** (experiment) → only that experiment's runs.
- **SRR / ERR** (run) → exactly that one run.
- **GSM** → that sample's runs. GSM labels are preserved as the `sample` column.

`source` per run:
- **`ena`** — ENA mirrors a usable barcode+cDNA FASTQ pair → curl download + md5 (fast, preferred).
- **`sra`** — not in ENA, or ENA only mirrors the single cDNA file (barcode read
  submitted as "technical") → prefetch + `fasterq-dump --include-technical`.

### Auto-chemistry detection (runs automatically, per-sample)
`prepare_runs.py` peeks BOTH reads via HTTP range request and finds which is the
barcode (26-28bp) vs cDNA (50-100bp). For each ENA sample:
- 28 bp barcode → **v3** (or v4/GEM-X — same geometry) — writes `chemistry: v3`
- 26 bp barcode → **v2** — writes `chemistry: v2`
- No 26-28bp read found → likely not 10x droplet (bulk/Smart-seq) — **refuses** to write sheet
- SRA-only runs → `chemistry: ""` (empty; detected after download, not yet implemented)

The Snakefile reads the `chemistry` column and dynamically selects whitelist + umi_len
per sample (`CHEM_PARAMS` maps `v2`/`v3`/`v4`/`arc` → whitelist + UMI length). Mixed-chemistry
datasets (v2 + v3 samples in one study) work automatically. Falls back to `config["chemistry"]`
for empty or unrecognized chemistry values.

**10x Multiome (ARC) / single-nucleus — pass `--multiome`.** Multiome GEX has the SAME 28bp
geometry as v3/v4 (16 CB + 12 UMI), so it is **not** auto-detectable — the length detector would
call it `v3`. For a Multiome dataset, run `prepare_runs.py <ACC> --multiome`: this stamps
`chemistry: arc` on **every** row (STARsolo then uses whitelist `737K-arc-v1.txt`, from
cellranger-arc) instead of hand-editing the column. Then set **`feature: GeneFull`** in
`config/config.yaml` (intronic reads dominate nuclei — `prepare_runs.py` prints this reminder).
If a Multiome run were left at `v3`, its valid-barcode fraction would collapse (wrong whitelist)
and the QC gate would fail it.

> **Why chemistry, not "snRNA", is the detectable thing.** `--multiome`/`feature: GeneFull` are
> two orthogonal decisions and NEITHER is detectable from raw-read *length*:
> - **feature (cells vs nuclei → Gene vs GeneFull)** is a library-prep fact, invisible in the
>   FASTQ — it only shows up *after* alignment as a high intronic fraction (GeneFull ≫ Gene
>   counts). So it stays a manual call from the GEO/SRA record ("single-nucleus"/"snRNA"/"nuclei").
> - **chemistry (arc vs v3 vs v4)** *is* recoverable from reads, but by matching barcode
>   **sequences** to each candidate whitelist (the right list matches ~all observed CBs, the
>   wrong ones ~none) — NOT from the 28bp length, which is identical across all three. This is
>   how Arc's scRecounter picks a whitelist (empirical parameter scan). The current detector only
>   peeks length, so arc needs the `--multiome` assertion; a sequence-matching auto-detector is a
>   possible future add.

### Strand (3′ vs 5′) — set explicitly, NOT auto-detected
`--soloStrand` must be `Forward` for 10x **3′** GEX and `Reverse` for 10x **5′** GEX.
This is **not** auto-detected: barcode/UMI geometry is identical for 3′ and 5′, so read
length cannot tell them apart. `prepare_runs.py` writes `strand: Forward` in every row;
**flip the rows that are 5′ to `Reverse` by hand** (or set `chemistry.strand: Reverse`
in config for an all-5′ study). The Snakefile reads the per-sample `strand` column and
falls back to `config["chemistry"]["strand"]` (default `Forward`) when it is blank.

How to decide the value:
- **Metadata**: check the GEO/SRA record — "Chromium … 5′", "5′ GEX", "5 prime", "VDJ 5′"
  → `Reverse`; "3′" or unstated → `Forward`.
- **Empirical (definitive)**: after a run, read `star/<sample>/Solo.out/Gene/Summary.csv`
  → "Reads Mapped to Gene: Unique+Multiple Gene". Correct strand ≈ 0.5–0.7; a wrong
  strand collapses it well below 0.1. The QC gate (`qc.min_reads_mapped_gene`) already
  fails a wrong-strand run, so a mistake is caught rather than silently shipped.

**Chemistry guard (pre-write check, can be skipped with `--no-chem-check`):**
- Refuses symmetric-mate data (both reads ~same length → not 10x droplet)
- Warns if detected chemistry mismatches config's expected chemistry (e.g. v2 data but config set for v3)
- symmetric equal-length mates (e.g. 101/101, 151/151) → **refuses, exit 1, no sheet
  written**: this is bulk / full-length / Smart-seq, not droplet. Use a different
  pipeline (plain STAR alignReads + featureCounts), not this one.
- SRA-only rows can't be peeked cheaply → note printed; verify after first download (Phase 4).

Flags: `--chem v2` (config is set for v2, adjust the warn threshold), `--no-chem-check` (skip the guard).

### Feature: Gene (scRNA) vs GeneFull (snRNA) — config toggle
`config/config.yaml` key **`feature`** selects the STARsolo count feature for one run:
- **`Gene`** (default) — exonic counts, correct for standard **scRNA-seq** (single-cell).
- **`GeneFull`** — counts over the **full gene body incl. introns**, required for
  **single-nucleus (snRNA-seq)**: nuclear RNA is largely unspliced pre-mRNA, so exon-only
  counting throws away most single-nucleus signal.

Set `feature: GeneFull` for single-nucleus datasets; leave `Gene` otherwise. The `starsolo`
rule emits `--soloFeatures {feature} Velocyto`, so the spliced/unspliced/ambiguous layers are
produced in **both** modes, and STARsolo writes matrices under `Solo.out/<feature>/`. The QC
gate reads the matching "Reads Mapped to Gene**Full**" column automatically (no extra config).
This is a per-run toggle; it is **not** auto-detectable — decide it from the assay (single-cell
vs single-nucleus / nuclei prep) in the GEO/SRA record. Design precedent: Arc Institute's
scRecounter uses `GeneFull` throughout for exactly this reason.

**Multiple / custom sample lists:** either run `prepare_runs.py` per accession and
concatenate the TSVs (one header), or hand-write rows. Back up any existing sheet
first (`cp config/samples.tsv config/samples.<tag>.bak.tsv`).

## Phase 2 — Launch on SLURM

`run_slurm.sh` runs the Snakemake controller as a small job; it activates the venv,
puts fastqc + sra-tools on PATH, and calls `snakemake --profile profiles/slurm`.
The SLURM executor submits each rule as its own right-sized job (STARsolo → big_memory/200G).

```bash
cd "$PIPE"
sbatch run_slurm.sh                 # whole sheet, all samples
./run_slurm.sh -n                   # dry-run on the login node first (recommended)
```

**Smoke-test ONE sample end-to-end** (5 jobs: check_versions → download_fastq →
fastqc → starsolo → to_h5ad) by targeting that sample's outputs directly:
```bash
sbatch run_slurm.sh \
  results/h5ad/<SAMPLE>.h5ad \
  results/qc/fastqc/<SRR>_R1_fastqc.zip \
  results/qc/fastqc/<SRR>_R2_fastqc.zip
```
Do **not** add `qc_gate` / `merge` / `multiqc` to a single-sample target — those
rules `expand` over ALL samples and won't resolve.

## Phase 3 — Watch progress (don't confuse "running" with "hung")

```bash
squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"      # ctl + per-rule jobs
LOG=$(ls -t "$PIPE"/.snakemake/log/*.snakemake.log | head -1); tail -20 "$LOG"
```
Per-rule slurm logs: `$PIPE/.snakemake/slurm_logs/rule_<name>/<wildcards>/<jobid>.log`.

The account-guessing lines (`No SLURM account given, trying to guess` /
`sacct: invalid option`) are **cosmetic** — the executor proceeds without an account.
They are NOT the hang. A download job legitimately takes 20–40 min (prefetch +
fasterq-dump + gzip of a multi-GB run); growing `results/fastq/*.fastq.gz` = healthy.

## Phase 4 — Verify a completed run

```bash
# 1. Barcode read length (proves chemistry parsed): expect 28 bp for v3
zcat "$PIPE"/results/fastq/<SRR>_R1.fastq.gz | head -2 | awk 'NR==2{print length($0)" bp"}'

# 2. STARsolo valid-barcode fraction should be healthy (not ~0)
grep -i "Reads With Valid Barcodes" "$PIPE"/results/star/<SAMPLE>/Solo.out/Gene/Summary.csv

# 3. h5ad has the spliced/unspliced/ambiguous layers
python - <<'PY'
import anndata as ad
a = ad.read_h5ad("<PIPE>/results/h5ad/<SAMPLE>.h5ad")
print(a.shape, "layers:", list(a.layers))   # expect spliced/unspliced/ambiguous
PY
```
A near-zero valid-barcode fraction means the chemistry/whitelist is wrong (e.g. v2
data run with a v3 config) — go back to Phase 1 and fix `config/config.yaml`.

The resulting `.h5ad` (Gene matrix in `X`, velocyto layers) is ready for
downstream scanpy analysis (QC, clustering, UMAP, DE).

## Troubleshooting (all fixes stay inside `$PIPE`)

| Symptom | Cause | Fix |
|---|---|---|
| Controller submits first job then never advances, though jobs succeed | SLURM executor polls `sacct`; slurmdbd down (`Connection refused …:6819`) | `profiles/slurm/config.yaml`: `slurm-status-command: squeue` + `slurm-no-account: true` (already set). `squeue` is DB-independent; cluster `MinJobAge=300s` keeps finished jobs pollable. |
| `LockException` on start | stale `.snakemake` locks from a killed controller | Verify no live controller (`squeue`), then `snakemake --profile profiles/slurm --unlock` and resubmit. |
| `prepare_runs.py` HTTP 400 from ENA | transient API error | Re-run; it usually succeeds on retry. |
| SRX/SRR gives whole-study runs | scoping regression | `build_rows` restricts sub-study accessions; confirm `is_substudy` branch in `prepare_runs.py`. |
| ENA has only 1 file for a 10x run | barcode read marked "technical", not mirrored | Expected — that run is routed to `source=sra` automatically. |
| STARsolo matrices near-empty | chemistry mismatch (v2 data, v3 config) | Set config to detected chemistry (v2: `umi_len 10`, 737K whitelist) and rerun. |
| STARsolo aborts `EXITING because of FATAL ERROR: Genome version … is INCOMPATIBLE` | a hand-supplied `star_index` was built by a different STAR (e.g. a 10x refdata `star/` = 2.7.1a vs your 2.7.10a) | Should not happen: `resolve_star_index` validates `genomeVersion` and rebuilds on mismatch. If you see it, the aligner read an index other than `<workdir>/star_index` — confirm the `starsolo` rule's `index` input is `STAR_INDEX`, then delete `<workdir>/star_index` and rerun so the resolver regenerates it. |
| Rebuilt STAR index disappears / rebuilds every run | `star_index: ""`, so the resolver builds into the ephemeral `<workdir>/star_index` | Point `reference.star_index` at a writable dir to keep the built index in (e.g. a sibling `STAR_rebuild/`). The resolver builds into that path once, then just symlinks it on later runs. |
| High genome mapping but ~3% reads mapped to **gene** | someone re-added `--soloBarcodeMate`/`--clip5pNbases` — that clips the separate barcode read down to 0bp of cDNA | Use the standard two-file `starsolo` rule (cDNA read first, barcode read second, no clip). See the comment in `workflow/Snakefile`. |
| fasterq-dump "cannot connect to external services" on login node | run on login node | Harmless for the rule — it prefetches the `.sra` first, then fasterq-dumps the local file on the compute node (no network). |

## Decision: skill vs. rewriting the pipeline

Keep Snakemake as the execution engine. It provides the DAG, resume, SLURM
right-sizing, and per-sample parallelism that an agent driving raw steps would
have to re-implement and babysit across multi-hour 200 GB jobs. This skill wraps
that engine with the judgment layer (resolution, guard, launch, diagnosis). The
failures seen in practice were infra (slurmdbd down), input (non-10x data), or
operational (stale locks) — none argue for abandoning Snakemake.

## Key files in `$PIPE`
- `workflow/scripts/prepare_runs.py` — accession → sample sheet + chemistry guard.
- `workflow/scripts/ncbi_utils.py` — GEO/SRA metadata + ENA URL resolution (vendored).
- `workflow/scripts/starsolo_to_h5ad.py` — STARsolo Gene/GeneFull + Velocyto feature → single .h5ad (`--feature Gene|GeneFull`).

> **"Velocyto" here = STARsolo's `--soloFeatures Velocyto` mode, NOT the separate
> velocyto.py tool.** STARsolo itself classifies each read as spliced / unspliced /
> ambiguous (exon vs intron overlap) and writes three matrices under
> `Solo.out/Velocyto/`. No separate velocyto tool runs. `starsolo_to_h5ad.py` loads
> those into `adata.layers['spliced'/'unspliced'/'ambiguous']`, alongside the
> standard `Gene/` matrix as `adata.X`.
- `workflow/Snakefile` — rules; `download_fastq` branches on `source`.
- `profiles/slurm/config.yaml` — SLURM profile + per-rule resources + sacct/squeue fix.
- `config/config.yaml` — chemistry (umi_len, whitelist), `feature` (Gene/GeneFull), reference paths.
- `run_slurm.sh` — controller launcher (venv + PATH + `snakemake --profile`).
