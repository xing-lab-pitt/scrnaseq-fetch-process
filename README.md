# scrnaseq-fetch-process

A Snakemake pipeline that turns a public **GEO/SRA accession** into per-sample
count matrices with RNA-velocity layers, using **STARsolo**:

```
accession ──prepare_runs.py──> samples.tsv
                                   │
 check_versions (preflight: env must satisfy requirements)
                                   │
 download_fastq ─ fastqc          resolve_star_index (symlink if compatible, else build)
        │                              │
        └────────────► starsolo (per sample) ──► qc_gate ──► qc_check (strict)
                             │                       │
                        to_h5ad (per sample) ────► merge ──► [optional downstream]
                             │
                          multiqc report
```

Output per sample: an `.h5ad` with the standard Gene count matrix as `adata.X`
and `spliced` / `unspliced` / `ambiguous` layers (from STARsolo's `Velocyto`
feature — this is STARsolo itself, **not** the separate velocyto.py tool).
`adata.X` is STARsolo's **`Gene`** feature (exonic UMI counts — or **`GeneFull`**,
intron-inclusive, when `feature: GeneFull`), the matrix you normalize and cluster on.
`Gene` and `Velocyto` are independent counting modes, so `X` is **not** spliced,
unspliced, or their sum — the layers sit alongside it for RNA velocity.

## Features

- **Any accession type** — `prepare_runs.py` resolves GSE / SRP / PRJNA / SRX /
  SRR / GSM into a sample sheet, scoping correctly (a study gives all its runs;
  an experiment/run gives just its own).
- **ENA-first download with SRA fallback** — per run, prefers ENA's pre-split,
  md5-checked FASTQs over HTTPS (`curl`); falls back to NCBI `prefetch` +
  `fasterq-dump --include-technical` for runs ENA hasn't mirrored. Chosen
  automatically per run (the `source` column); no global switch.
- **Automatic chemistry detection (v2 / v3 / v4 / arc)** — `prepare_runs.py` learns
  each run's 10x chemistry from remote metadata, no download, via three channels:
  barcode-read **length** (ENA peek), a barcode-**sequence** whitelist probe (rescues
  over-sequenced R1), and the **processing software** named in the run's GEO
  `data_processing` record (reaches SRA-only runs). It writes the result per run and
  refuses genuinely non-10x-droplet data (bulk / Smart-seq) with a precise reason.
- **10x Multiome (ARC) / single-nucleus** — the `arc` chemistry is supported and
  auto-detected (or asserted with `--multiome`); pair with `feature: GeneFull` for
  nuclei. Over-sequenced barcode reads (e.g. a 150 bp R1) are handled automatically
  via `--soloBarcodeReadLength 0`.
- **Per-run toggles** — chemistry, strand (3′/5′), feature (Gene/GeneFull), and the
  over-sequenced-R1 flag are all per-sample sheet columns, so mixed studies (e.g. v2 +
  v3 samples) process correctly in a single run.
- **Deterministic QC gate** — pass/fail thresholds live in `config.yaml`
  (version-controlled), applied to STARsolo's `Summary.csv`.
- **STAR index resolver** — the index is a pipeline output: `resolve_star_index`
  probes what your STAR speaks (`genomeVersion`) and symlinks a compatible supplied
  index or rebuilds from FASTA+GTF, so an index/aligner version mismatch is impossible
  by construction.
- **SLURM-ready** — each rule is submitted as a right-sized job; the included
  profile is resilient to a down `slurmdbd` (polls `squeue`, not `sacct`).

## Requirements

**Python** (a venv or conda env): `snakemake` ≥ 9.24 with
`snakemake-executor-plugin-slurm`, plus the analysis pins in
`requirements-pipeline.txt` (`anndata`, `scanpy`, `numpy`, `pandas`, `scipy`,
`h5py`, `pyyaml`, `requests`). `check_versions.py` enforces these at preflight.

**Binary tools** (not pip; from bioconda, modules, or system): `STAR` (2.7.10a
tested), `samtools`, `FastQC`, `MultiQC`, and — only for `source=sra` runs —
`sra-tools` (`prefetch`, `fasterq-dump`). Conda users can let Snakemake manage
these via `--use-conda` (env specs in `workflow/envs/`).

**Reference data** (not shipped — obtain separately):
- A STAR index, **or** a genome FASTA + gene GTF for the pipeline to build one
  (`star_index` rule). Two easy sources:
  - **10x Cell Ranger reference bundles** (convenient — one download gives FASTA +
    GTF, already gene-filtered to match 10x conventions):
    - Human: <https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz>
    - Mouse: <https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz>

    Each bundle extracts to `refdata-gex-<genome>-2024-A/` containing
    `fasta/genome.fa` and `genes/genes.gtf.gz` (plus a prebuilt `star/` index).
    Point `reference.fasta` / `reference.gtf` in `config.yaml` at these (gunzip
    the GTF first — the pipeline expects an uncompressed `.gtf`). You *may* point
    `reference.star_index` at the bundle's `star/` index: `resolve_star_index`
    checks its `genomeVersion` against your STAR and only reuses it if compatible,
    otherwise it rebuilds from fasta+gtf automatically. (Cell Ranger 2024-A ships
    a STAR 2.7.1a index; a 2.7.10a STAR can't load it, so the resolver rebuilds —
    no manual intervention, no "Genome version INCOMPATIBLE" abort.) When a rebuild
    is needed and `reference.star_index` names a *writable* path (e.g. a sibling
    `STAR_rebuild/` dir), the resolver builds the fresh index **into that path** so
    it persists — the next run finds a compatible index there and just symlinks it
    (no repeated ~30-min rebuild). Leave `reference.star_index` as `""` to always
    build into the (ephemeral) workdir instead.
  - **Ensembl** (raw): e.g. GRCh38 primary-assembly FASTA + release GTF, if you
    prefer unfiltered annotation.
- A 10x barcode whitelist matching your chemistry (see [Chemistry](#chemistry)).
  These ship with Cell Ranger (`lib/python/cellranger/barcodes/`); the
  [Teichmann lab scg_lib_structs](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)
  project also mirrors them:
  - v3 / v3.1: `3M-february-2018.txt`
  - **v4 / GEM-X**: `3M-3pgex-may-2023.txt`
    ([download](https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-3pgex-may-2023.txt.gz))
  - v2: `737K-august-2016.txt`
  - **arc (10x Multiome GEX)**: `737K-arc-v1.txt` — ships with Cell Ranger ARC /
    Cell Ranger (`lib/python/cellranger/barcodes/`); not in the raw refdata bundles.

  Put all the whitelists you might need in one directory; `prepare_runs.py` reads it
  (`--whitelist-dir`, default `…/reference/10x_whitelists`) for the sequence probe,
  and the Snakefile picks the one matching each sample's detected chemistry.

## Setup

```bash
git clone <your-repo-url> scrnaseq-fetch-process
cd scrnaseq-fetch-process

# 1. Get a reference (example: 10x human GRCh38-2024-A). Mouse: GRCm39-2024-A.
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
gunzip -k refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz   # pipeline needs an uncompressed .gtf

# 2. Create your config from the template and edit the paths in it.
cp config/config.example.yaml config/config.yaml
$EDITOR config/config.yaml     # set reference.* (fasta/gtf above; leave star_index
                               #   an empty dir to build), chemistry.whitelist, sra_tools_bin

# 3. (SLURM) point run_slurm.sh at your environment via env vars, and edit the
#    partition names in profiles/slurm/config.yaml to match your cluster.
export SCRNASEQ_VENV=/path/to/your/venv           # optional if snakemake is already on PATH
export SCRNASEQ_EXTRA_PATH=/opt/FastQC:/opt/sratoolkit/bin   # tool dirs not in the venv
```

`config/config.yaml` and any real `samples.tsv` are gitignored, so your local
paths never get committed.

## Usage

```bash
# 1. Resolve an accession to a sample sheet (network step, run once).
python workflow/scripts/prepare_runs.py GSE123456 -o config/samples.tsv
#    SRX/SRR/GSM/SRP/PRJNA also work. Chemistry is auto-detected; add --chem v2 if
#    your config is v2, or --multiome for 10x Multiome / single-nucleus data.

# 2. Dry-run, then submit on SLURM (controller runs as a small job; each rule
#    becomes its own job).
sbatch run_slurm.sh -n        # dry run
sbatch run_slurm.sh           # full run

# ...or run Snakemake directly (e.g. on a workstation, with tools on PATH):
snakemake --profile profiles/slurm -n
```

### Smoke-test one sample first (recommended)

Target one sample's `.h5ad` — this pulls the whole scientific path
(check_versions → download_fastq → fastqc → starsolo → to_h5ad) without the
aggregating rules (qc_gate/merge/multiqc `expand` over ALL samples, so they
can't run on a single sample):

```bash
sbatch run_slurm.sh \
  results/h5ad/<SAMPLE>.h5ad \
  results/qc/fastqc/<SRR>_R1_fastqc.zip \
  results/qc/fastqc/<SRR>_R2_fastqc.zip
```

Then verify: barcode read ~28 bp (v3) —
`zcat results/fastq/<SRR>_R1.fastq.gz | head -2` — and the `.h5ad` has
`spliced`/`unspliced`/`ambiguous` layers.

## Chemistry

Defaults are 10x Chromium 3′ **v3**. Supported chemistries:

| Chemistry | Barcode read | CB / UMI | `umi_len` | Whitelist |
|---|---|---|---|---|
| v2 | 26 bp | 16 / 10 | `10` | `737K-august-2016.txt` |
| **v3** / v3.1 (default) | 28 bp | 16 / 12 | `12` | `3M-february-2018.txt` |
| **v4** / GEM-X | 28 bp | 16 / 12 | `12` | `3M-3pgex-may-2023.txt` |
| **arc** (Multiome GEX) | 28 bp | 16 / 12 | `12` | `737K-arc-v1.txt` |

The barcode+UMI read and the cDNA read are separate FASTQs, so STARsolo runs in
the standard two-file mode — the barcode read is never clipped, and there is no
`clip5p`/`barcode_mate` setting to tune.

### Automatic detection (three channels)

`prepare_runs.py` fills the `chemistry` column from remote metadata (no download),
and the Snakefile maps it to the right whitelist + `umi_len` per sample — so
**mixed-chemistry studies process in one run**. It falls back to `config.yaml`'s
static `chemistry` for any row it leaves empty.

1. **Barcode-read length (ENA)** — peeks both reads via HTTP range request: 26 bp
   barcode → v2, 28 bp → v3/v4/arc.
2. **Barcode-sequence probe (ENA)** — when no mate is a barcode length (e.g. an
   **over-sequenced R1**: a 150 bp read that ran through the 28 bp barcode into the
   poly-T tail), length alone can't tell 10x from bulk/Smart-seq. The probe samples
   ~2000 reads and matches their leading 16 bp against each whitelist (streamed, so
   the multi-million-line lists are never fully loaded). A hit both proves droplet
   data and names the chemistry, and flags the run `barcode_read_length=0` →
   `--soloBarcodeReadLength 0` (STARsolo reads CB/UMI from the leading positions and
   ignores the tail). A miss → refuse.
3. **GEO `data_processing` software** — the processing software the submitter named
   in the run's GEO record (`cellranger-arc`, `starsolo`, … vs non-droplet
   `smartseq`/`rsem`/`hisat2`/…), read once per study. This is the only channel that
   reaches **SRA-only runs** the FASTQ probe can't peek: it auto-detects `arc`, fills
   SRA-only chemistry, and drives a precise non-droplet refusal that names the software.

Guardrails: a confident v2 (26 bp) is never overridden; a sequence probe that
*actively rejects* a run beats a `cellranger-arc` metadata label (data > metadata).
Wrong chemistry (or the wrong whitelist for the right geometry) silently yields
near-empty matrices — the QC gate's `min_valid_barcodes` catches it.

### v3 vs v4

Identical read geometry — the *only* difference is the whitelist. To run v4 data,
keep the v3 geometry and just point `chemistry.whitelist` at `3M-3pgex-may-2023.txt`.
Because read length can't distinguish them, detection reports both as "v3"; pass
`--chem v4` to `prepare_runs.py` and it prints a reminder to confirm the whitelist
rather than a false mismatch warning.

### 10x Multiome (ARC) / single-nucleus

Multiome GEX has the same 28 bp geometry as v3/v4, so it can't be told apart by
**length** — but it is resolved to `arc` automatically when the GEO record names
`cellranger-arc` (channel 3) or the barcode-sequence probe matches `737K-arc-v1.txt`
(channel 2). You can also assert it explicitly with `prepare_runs.py <ACC>
--multiome` (stamps `chemistry: arc` on every row). For nuclei, also set
**`feature: GeneFull`** in `config.yaml` — `prepare_runs.py` prints this reminder.
Leaving Multiome data at `v3` uses the wrong whitelist, collapses the valid-barcode
fraction, and the QC gate fails it.

### `prepare_runs.py` chemistry flags

| Flag | Effect |
|---|---|
| `--multiome` | Assert `chemistry: arc` on every row (10x Multiome / single-nucleus). |
| `--chem {v2,v3,v4,arc}` | Expected chemistry for the guard's warn threshold. |
| `--no-chem-check` | Skip the pre-write chemistry guard entirely. |
| `--no-geo-software` | Skip the GEO `data_processing` lookup (offline / fallback to length+probe). |
| `--whitelist-dir DIR` | Directory of 10x whitelists used by the sequence probe. |

## Feature & strand (per-run settings)

- **`feature`** (`config.yaml`) — `Gene` = exonic counts (standard **scRNA-seq**,
  default); `GeneFull` = full gene body incl. introns (**single-nucleus / snRNA-seq**,
  where unspliced nuclear RNA dominates). Both modes still emit the Velocyto
  spliced/unspliced/ambiguous layers. Not detectable from raw reads (it only shows as
  a high intronic fraction after alignment), so it's a manual call from the record.
- **`strand`** (per-sample sheet column, or `chemistry.strand` default) — `Forward`
  for 10x **3′** GEX, `Reverse` for 10x **5′** GEX. Not auto-detectable (3′ and 5′
  share the same barcode geometry); set it from the GEO/SRA record. A wrong call
  collapses the "reads mapped to gene" fraction, which the QC gate catches.

## Batch completeness & the agent loop

Snakemake calls a sample done as soon as its `.h5ad` merely *exists*. A thin
operational layer on top — `reconcile.py` + `run_batch.sh` — enforces the fuller
contract, records what actually mapped, and drives reruns across many studies. It
keeps **mechanics in scripts and judgment with the human**; nothing here downloads
or aligns — it only ever hands re-runnable work back to Snakemake.

**"Done" contract.** A sample is *done* iff **all three** hold: its
`<workdir>/h5ad/<sample>.h5ad` exists, `<sample>` is in `<workdir>/qc/qc_pass.txt`
(passed the alignment QC gate), and that h5ad has layers `spliced`/`unspliced`/
`ambiguous`. Everything else is categorized, and each category carries an
**action** — **rerun** (transient, worth an automatic retry) or **flag** (genuine;
re-running the same inputs gives the same result, so a human decides):

| category       | meaning                                                  | action |
|----------------|----------------------------------------------------------|--------|
| `missing`      | no h5ad (not run yet / upstream crash)                   | **rerun** — hand back to Snakemake |
| `corrupt`      | h5ad exists but unreadable (truncated / killed write)    | **rerun** — quarantined first, then regenerated |
| `read_qc_fail` | raw reads failed read-QC (bad library / wrong chemistry) | **flag** — same reads re-align to the same failure |
| `qc_fail`      | h5ad exists but not in `qc_pass.txt` (failed align gate) | **flag** — reason in `qc/qc_gate.tsv` |
| `no_layers`    | QC passed but a velocity layer is absent                 | **flag** — silent mis-wire in `Solo.out/Velocyto/` |

**Read-QC axis (raw-read quality).** The Snakemake `read_qc` rule (`read_qc.py`)
parses each run's FastQC output into `qc/read_qc.tsv` + `qc/read_qc_pass.txt`. It is
**role-aware**: R1 = barcode read (short/repetitive — FastQC content/GC/kmer FAILs
are *by design*, judged leniently), R2 = cDNA read (judged strictly on read count,
mean base quality, adapter). Always non-fatal (report layer). If a run genuinely
fails, the reconciler marks the sample `read_qc_fail` up front — rerunning bad reads
can't help. If `qc/read_qc.tsv` is absent (older workdir), the axis is simply skipped.

### Files to check

| file | what it tells you |
|------|-------------------|
| `results/successful_samples.tsv` | **the success logbook** — one append-only, idempotent row per *done* `(accession, sample)` with its mapping metrics, chemistry, feature, SRRs, read-QC status. The durable record of what mapped. |
| your `--report` TSV (e.g. `qc_reconcile.tsv`) | every sample's `category` + `action` + reason for the latest reconcile. |
| `<workdir>/qc/qc_gate.tsv` | alignment metrics + PASS/FAIL reason (valid barcodes, reads-mapped-gene, cells, saturation). |
| `<workdir>/qc/read_qc.tsv` | per-run raw-read verdict (R1/R2 length, cDNA quality, adapter). |
| `<workdir>/.quarantine/<stamp>/` | corrupt outputs moved aside (never deleted) so Snakemake regenerates them. |
| `config/manifest.tsv` | one row per study (`accession  feature  workdir  [samples_tsv]  [notes]`); comment lines allowed. Copy `config/manifest.example.tsv` to start. |

### Use it

```bash
source "$SCRNASEQ_VENV/bin/activate"    # any env with h5py + pyyaml (snakemake's is fine)
cp config/manifest.example.tsv config/manifest.tsv   # then edit for your studies

# Check ONE workdir (read-only; launches nothing). exit 0 = done, 1 = work remains, 2 = error.
python workflow/scripts/reconcile.py --samples config/samples.tsv \
       --workdir results/GSE123456 --accession GSE123456 --feature Gene \
       --report qc_reconcile.tsv --ledger results/successful_samples.tsv

# Check EVERY study in the manifest, recording successes to the logbook:
python workflow/scripts/reconcile.py --manifest config/manifest.tsv --base-dir "$PWD" \
       --json recon.json --ledger results/successful_samples.tsv

# ONE batch cycle: reconcile -> record ledger -> rerun missing + quarantine/rerun
# corrupt -> flag the rest. NOT a daemon; call it repeatedly.
SCRNASEQ_PY=$(command -v python) ./run_batch.sh          # add DRY_RUN=1 to preview (no sbatch)
PIPE=$PWD ./snakemake_status.sh                          # watch jobs (or: squeue -u $USER)
```

### The loop (thin, two human gates)

```
GATE 1 (human) ─ approve config/manifest.tsv; ensure each study has a sample sheet
                 (prepare_runs.py per accession — chemistry guard wants a human eye)
      │
      ▼
  ┌─ run_batch.sh → reconcile → record ledger → rerun missing + quarantine/rerun
  │     │            corrupt → flag read_qc_fail / qc_fail / no_layers
  │     ▼
  │   watch snakemake_status.sh (or squeue) until controllers are idle
  └─────┘  repeat while run_batch.sh exits 1 AND it launched something this cycle
      │
      ▼
  stop when: reconciler COMPLETE (exit 0), OR only flagged / retry-capped studies remain
      │
      ▼
GATE 2 (human) ─ review the reconcile report + qc_gate.tsv + read_qc.tsv; decide on
                 flagged studies (fix chemistry & re-prepare, accept, or drop).
                 Successful samples are already in the logbook.
```

- **Why `prepare_runs.py` is Gate 1, not in the loop:** building a sheet does
  network metadata calls + a chemistry guard that wants human confirmation. The loop
  only relaunches Snakemake for studies whose sheet already exists — it never guesses
  chemistry. A study with no sheet shows up flagged with "samples sheet not found".
- **Retry cap:** a study with `missing`/`corrupt` samples is relaunched at most
  `MAX_RETRIES` (default 2, counted in `.batch_state.json`), then escalated to the
  human instead of looping forever.
- **Tests:** `tests/test_reconcile.py` (every category + action + ledger idempotency)
  and `tests/test_read_qc.py` (offline role-aware read-QC verdict).

## Repository layout

```
config/
  config.example.yaml     # template — copy to config.yaml and edit
  samples.example.tsv     # sample-sheet schema example
  manifest.example.tsv    # batch manifest schema — copy to manifest.tsv (agent loop)
workflow/
  Snakefile               # the DAG
  scripts/                # prepare_runs, ncbi_utils, starsolo_to_h5ad, qc_gate,
                          #   read_qc, reconcile, ...
  envs/                   # conda env specs (scanpy.yaml, tools.yaml)
tests/                    # offline unit tests (test_reconcile.py, test_read_qc.py)
profiles/slurm/           # SLURM profile (edit partitions for your cluster)
run_slurm.sh              # SLURM launcher (paths from env vars / script location)
run_batch.sh              # batch agent loop: reconcile -> rerun/flag across studies
snakemake_status.sh       # watch running SLURM jobs
requirements-pipeline.txt # Python version pins (enforced at preflight)
skill/                    # optional Claude Code agent runbook (see below)
```

## Claude Code skill (optional)

`skill/SKILL.md` is an agent runbook for driving this pipeline with
[Claude Code](https://claude.com/claude-code) — accession resolution, launch,
hang diagnosis, and post-run verification. To use it, symlink it into your
Claude skills directory:

```bash
ln -s "$(pwd)/skill" ~/.claude/skills/scrnaseq-fetch-process
```

It's entirely optional — the pipeline runs fine without Claude Code.

## License

MIT — see [LICENSE](LICENSE). Update the copyright line with your name/lab.
