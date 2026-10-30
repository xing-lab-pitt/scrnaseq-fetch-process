# How this pipeline compares to bollito and scRecounter

All three pipelines run STARsolo on 10x scRNA-seq. They aim at different jobs:

| | goal |
|---|---|
| **scrnaseq-fetch-process** (this repo) | Public accession -> checked per-sample and merged `.h5ad` with velocity layers, on any SLURM cluster |
| **[bollito](https://github.com/cnio-bu/bollito)** | Local FASTQs or matrices -> full Seurat analysis (QC, clustering, DE, trajectories) |
| **[scRecounter](https://github.com/ArcInstitute/scRecounter)** | Bulk re-processing of SRA for Arc Institute's scBaseCamp atlas, driven from their cloud database |

Versions compared: bollito `4ea46b0` (last commit 2023-04-12), scRecounter `92891ee`
(2026-07-24; last code change about 2025-02). Each claim below was checked against the
source code, not just the README.

## At a glance

| | this pipeline | bollito | scRecounter |
|---|---|---|---|
| **Input** | GSE / GSM / SRP / PRJNA / SRX / SRR, resolved to runs automatically | Local FASTQs or matrices only | SRX/SRR CSV, or a query of Arc's database. No GEO lookup |
| **Download** | ENA first (md5-checked), SRA fallback, chosen per run | None | `fasterq-dump`, falling back to `fastq-dump` capped at 200M read pairs |
| **Download integrity** | Truncation, mid-file record dropout and corrupt-prefetch checks; atomic writes; disk preflight; cluster-wide download limit | n/a | None beyond retries |
| **Chemistry (v2/v3/v4/arc)** | Detected before download from three sources (read length, barcode whitelist match, GEO processing notes); per sample, so mixed studies work | Set by hand; no v3.1/v4/arc | Trial alignment of 1M reads against 5 whitelists |
| **Over-sequenced barcode read (e.g. 150 bp R1)** | Detected and handled | No | Always `--soloBarcodeReadLength 0` |
| **Non-10x data** | Refused up front with a reason (bulk, Smart-seq) | Not checked | Only filtered when using Arc's database; a user CSV is not checked |
| **Strand (3′ / 5′)** | Per-sample setting, plus on-demand empirical test (`strand_probe`) | Not set (STAR default, 3′ only) | Tested on every sample; ambiguous samples are dropped |
| **snRNA (intron counting)** | `feature: GeneFull` | No | GeneFull variants always produced |
| **Read diagnostics before alignment** | `preflight`: predicts reads STAR will drop (TSO, adapters, poly-A/T/G) | FastQC / RSeQC only | None |
| **STAR index** | Pipeline checks the index matches the installed STAR, else rebuilds | Built without the GTF; no version check | Prebuilt paths only |
| **Cell calling** | Explicit `--soloCellFilter` (EmptyDrops_CR in the example); re-call without re-aligning | Raw matrix, keeping barcodes with >100 and then ≥200 detected genes; no EmptyDrops | EmptyDrops_CR |
| **QC gate** | Thresholds in config; strict mode stops the run | None on alignment | Valid-barcode fraction ≥ 0.3 in trials |
| **"Did every dataset finish?"** | `reconcile.py` + `run_batch.sh` | No | Rows in Arc's database |
| **Output** | Per-sample + merged `.h5ad`; `X` = Gene or GeneFull; spliced / unspliced / ambiguous layers | Seurat 3 `.rds` per step; `.h5ad` only after clustering; velocity only in R | Gzipped MTX per feature; `.h5ad` via a separate cloud pipeline |
| **Runs on a generic SLURM cluster** | Yes; site settings via env vars | No SLURM profile; an interactive y/N prompt can stall cluster jobs | No: every step connects to Arc's Google Cloud SQL database and Secret Manager |
| **Analysis after counting** | Optional scanpy step; analysis lives in separate skills | Extensive (DE, gene-set scoring, Slingshot, VISION, velocity) | None |
| **Organisms** | One reference per config (any species with FASTA + GTF) | One per config | 26 prebuilt indexes |

## Where this pipeline is stronger

### 1. It starts from the accession, not from files
Give it a GEO series and `prepare_runs.py` finds every run, picks ENA or SRA for each,
and writes the sample sheet. Scoping is exact: a series gives all its runs, a single
sample gives only its own. `run_slurm.sh` prints which config and workdir it is about to
run, so a wrong study is visible in the log.
bollito expects you to have the FASTQs already. scRecounter takes SRR lists, or rows from
a database only Arc can reach.

### 2. Wrong data is stopped before it costs an alignment
- **Chemistry** is learned from metadata before any download: barcode-read length, a
  barcode-sequence match against the 10x whitelists, and the software named in the GEO
  record. scRecounter instead spends trial alignments to find the whitelist, and bollito
  trusts the user.
- **Non-droplet data** (bulk, Smart-seq) is refused with a stated reason. bollito never
  checks. scRecounter only filters when it pulls from its own database.
- **`preflight`** predicts, from a read subsample, how many reads STAR will throw away as
  too short (for example reads starting with the 10x template-switch oligo, which FastQC
  does not detect). On one study this was 21.8% of reads and would otherwise have shown up
  only after a 6-hour alignment. Neither other pipeline looks for it.

### 3. Downloads are verified, not trusted
SRA and ENA transfers can be cut short or drop records partway through a file without any
error. This pipeline checks file integrity and the FASTQ record structure. It also writes
outputs atomically, checks node disk space first and limits concurrent downloads.
scRecounter retries and falls back to another tool. When it falls back, it silently caps
the run at 200M read pairs.

### 4. Version mismatches can't happen
The STAR index is a pipeline output. `resolve_star_index` reads which index format the
installed STAR expects, then links a compatible index or builds a new one. bollito and
scRecounter both trust whatever index path they are given.

### 5. Every sample is gated, and every study is accounted for
- A deterministic QC gate (thresholds in version-controlled config) checks each sample's
  `Summary.csv`. In strict mode a failing sample stops the run rather than sliding into
  the merged object.
- `reconcile.py` and `run_batch.sh` answer "did every dataset in this batch finish, and if
  not, why?" across many studies.
- Read-only checks (`inspect_qc.py`, `check_layers.py`) explain failures. For example,
  `check_layers.py` tells "velocity was never counted" apart from "it was counted but not
  attached to the h5ad".

### 6. Output is ready for Python analysis
One `.h5ad` per sample plus a merged one. `X` is Gene (or GeneFull for nuclei), and the
spliced / unspliced / ambiguous layers sit alongside it, ready for scanpy, scVelo or
dynamo. Gene IDs are kept unique, with symbols in `var['gene_symbols']`. bollito's velocity
stays inside an R object. scRecounter stops at MTX files.

### 7. It runs anywhere you have SLURM
No cloud account, no private database, no interactive prompts. Site settings (venv,
tool paths, partitions, whitelist directory) come from environment variables and the
SLURM profile.

### 8. Cell calling is a decision, not an accident
STAR's default cell filter assumes about 3,000 cells. On one 12-sample study it kept
38,556 cells, against 65,469 with `EmptyDrops_CR`; the EmptyDrops_CR calls matched the
authors' own cell calls at 99.93%. This pipeline makes the filter an explicit config key.
`recall_cells.py` re-applies a different filter to finished runs in minutes, with no
re-alignment.

## Where the others are stronger

We state these plainly, because they mark real gaps or deliberate scope limits.

- **STARsolo settings that mimic Cell Ranger** (scRecounter). It runs
  `--clipAdapterType CellRanger4 --soloUMIdedup 1MM_CR --soloCBmatchWLtype
  1MM_multi_Nbase_pseudocounts --soloMultiMappers EM`. Its counts are therefore closer to
  Cell Ranger output, and template-switch-oligo and poly-A tails are clipped rather than
  just detected. We keep STAR's defaults so existing studies stay comparable. Adopting
  these flags would change counts and would need a re-run of every study to avoid a
  pipeline-made batch effect.
- **Strand tested on every sample** (scRecounter). Our `strand_probe` exists but is
  opt-in. The default strand comes from the sample sheet.
- **Many organisms out of the box** (scRecounter, 26 prebuilt indexes). We support any
  species but build or link one reference per config.
- **Downstream analysis** (bollito). It goes from counts to clusters, markers, gene-set
  scores and trajectories in one workflow. We stop at the `.h5ad` by design; analysis is
  handled by separate scanpy / scVelo / dynamo tooling.
- **Matrix input** (bollito). It can start from existing count matrices. We start from
  reads (author-processed h5ad files can be fetched with `fetch_geo_supp.py`, outside the
  workflow).

## Which to use

- **Public GEO/SRA 10x data -> analysis-ready `.h5ad` on your own cluster:** this pipeline.
- **You already have FASTQs or matrices and want a full Seurat analysis in R:** bollito.
  Expect to set chemistry and whitelist by hand. Use it only for 3′ data and only with v1-v3
  chemistry.
- **Re-processing SRA at atlas scale on Arc's infrastructure:** scRecounter. Outside Arc,
  it needs code changes to remove the database dependency first.
