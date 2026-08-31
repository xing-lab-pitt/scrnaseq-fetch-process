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

> **Site activation — fill these in for your machine, and keep the values out of
> the repo.** Record them wherever you keep local notes (a gitignored file, your
> shell profile), not here.
>
> Run from a working copy that has a filled-in `config/config.yaml`, a built STAR
> index, and its own `results/` — not from a fresh clone, which ships `/path/to`
> placeholders:
> ```
> $PIPE = <your working copy of this pipeline>
> ```
> Environment `run_slurm.sh` reads:
> ```
> export SCRNASEQ_VENV=<venv with snakemake + the SLURM executor plugin>
> export SCRNASEQ_EXTRA_PATH=<dirs holding fastqc, sra-tools, … if not on PATH>
> export SCRNASEQ_WHITELIST_DIR=<dir of 10x barcode whitelists>
> ```
> STAR and samtools must be on `PATH` (often already in `/usr/bin`); check with
> `STAR --version`.
>
> Whitelist files the barcode probe looks for in `$SCRNASEQ_WHITELIST_DIR`
> (obtain from your Cell Ranger installation — `lib/python/cellranger/barcodes/`):
> ```
>   737K-august-2016.txt   (v2 — 26bp R1)
>   3M-february-2018.txt   (v3 — 28bp R1)
>   3M-3pgex-may-2023.txt  (v4/GEM-X — same geometry, swap whitelist only)
>   737K-arc-v1.txt        (10x Multiome/ARC GEX — snRNA; set chemistry: arc)
> ```
> References: point `reference.fasta` / `reference.gtf` at a 10x `refdata-gex-*`
> bundle (or your own genome + GTF). To keep a rebuilt index, set
> `reference.star_index` to a writable dir (e.g. a sibling `STAR_rebuild/`);
> `resolve_star_index` builds into it once, then reuses it.
>
> **Non-10x genome: build the index yourself first.** `resolve_star_index` builds
> with STAR's defaults, and two of them are wrong outside the 10x bundles. A small
> genome needs `--genomeSAindexNbases 12` (the default 14 is sized for ~3 Gb). A GTF
> that names genes in an attribute other than `gene_name` needs
> `--sjdbGTFtagExonParentGeneName <attr>` (FlyBase uses `gene_symbol`) — without it
> every gene name comes out blank and `adata.var` carries bare accession IDs, with no
> error anywhere. Check yours with `grep -m1 gene_ your.gtf`. Build once, then point
> `reference.star_index` at that directory so the resolver symlinks instead of
> rebuilding. The README has a worked example (Drosophila, FlyBase r6.31).

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
`prepare_runs.py` learns each run's chemistry from remote metadata — no download — using
**three complementary channels**, and writes the result into the `chemistry` (and, when
relevant, `barcode_read_length` / `software`) columns of the sheet.

**Channel 1 — barcode-read LENGTH (ENA).** Peeks BOTH reads via HTTP range request and finds
which is the barcode (26-28bp) vs cDNA (50-100bp):
- 28 bp barcode → **v3** (or v4/GEM-X / arc — same geometry) — writes `chemistry: v3`
- 26 bp barcode → **v2** — writes `chemistry: v2`
- No 26-28bp read among the mates → hand off to Channel 2 (below) before deciding.

**Channel 2 — barcode SEQUENCE probe (ENA).** When neither mate is a barcode length (e.g. an
**over-sequenced R1**: a 150 bp read that ran through the 28 bp barcode into the poly-T tail),
length alone can't tell 10x from Smart-seq/bulk. `probe_barcode_chem` samples ~2000 reads, takes
each read's leading 16 bp, and **streams** each candidate whitelist once (O(sample) memory — never
loads the 3M-line lists) to find the best match fraction. A real 10x read matches its own whitelist
at ~0.8; random 16-mers match at ~1e-4. So a hit both **proves droplet data** and **names the
chemistry** (incl. `arc`). On a hit the run gets `barcode_read_length: 0`, which the Snakefile turns
into `--soloBarcodeReadLength 0` so STARsolo reads CB/UMI from the fixed leading positions and
ignores the tail. On a miss, `chemistry` stays empty and the guard **refuses** (not 10x droplet).
The whitelist directory comes from `--whitelist-dir`, defaulting to `$SCRNASEQ_WHITELIST_DIR`.
Unset and not passed, this channel finds no whitelists and stays silent — the other two
channels still run.

**Channel 3 — GEO `!Sample_data_processing` software (reaches SRA-only runs).** GEO serves every
sample as machine-readable SOFT text naming the submitter's own processing software.
`fetch_geo_software` (run once per study, via the run→GSM map) classifies it — `cellranger-arc`,
`cellranger`, `starsolo`, `kb-python`, `alevin`, … vs non-droplet `smartseq`/`rsem`/`hisat2`/… —
and writes it to the `software` column. This is the only channel that reaches **SRA-only runs the
FASTQ probe can't peek**. It: (a) **auto-detects `arc`** — `cellranger-arc` pins the arc chemistry
that read length can't tell from v3/v4, so `--multiome` is no longer required for arc studies whose
record names the software; (b) **fills chemistry for SRA-only rows**; (c) drives a **precise
non-droplet refusal** (names the Smart-seq/bulk software). Guardrails: it never overrides a
confident v2 (26 bp), and a Channel-2 probe that *actively rejected* a run wins over a
"cellranger-arc" label (data beats metadata). Disable with `--no-geo-software` (offline/fallback).

The Snakefile reads the `chemistry` column and dynamically selects whitelist + umi_len per sample
(`CHEM_PARAMS` maps `v2`/`v3`/`v4`/`arc` → whitelist + UMI length). Mixed-chemistry datasets
(v2 + v3 samples in one study) work automatically. Falls back to `config["chemistry"]` for empty or
unrecognized chemistry values. `software` is provenance only — the Snakefile ignores it.

**10x Multiome (ARC) / single-nucleus.** Multiome GEX has the SAME 28bp geometry as v3/v4 (16 CB +
12 UMI), so it is not distinguishable by **length** — the length detector would call it `v3`. It is
now resolved to `arc` automatically by Channel 3 (GEO says `cellranger-arc`) or Channel 2 (barcodes
match `737K-arc-v1.txt`, from cellranger-arc). **`--multiome` remains as an explicit assertion**:
pass it to stamp `chemistry: arc` on **every** row when the metadata is silent or you want to force
it, instead of hand-editing the column. Either way, also set **`feature: GeneFull`** in
`config/config.yaml` (intronic reads dominate nuclei — `prepare_runs.py` prints this reminder).
If a Multiome run were left at `v3`, its valid-barcode fraction would collapse (wrong whitelist) and
the QC gate would fail it.

> **Why chemistry, not "snRNA", is the detectable thing.** `--multiome`/`feature: GeneFull` are
> two orthogonal decisions:
> - **feature (cells vs nuclei → Gene vs GeneFull)** is a library-prep fact, invisible in the raw
>   FASTQ — it only shows up *after* alignment as a high intronic fraction (GeneFull ≫ Gene
>   counts). So it stays a manual call from the GEO/SRA record ("single-nucleus"/"snRNA"/"nuclei").
> - **chemistry (arc vs v3 vs v4)** is NOT distinguishable by the 28bp length (identical across all
>   three), but IS recoverable two ways, both now implemented: from the barcode **sequences**
>   (Channel 2, matching each candidate whitelist — the same idea as Arc's scRecounter parameter
>   scan) and from the **stated software** (Channel 3, `cellranger-arc` in the GEO record).

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
- Warns if detected chemistry mismatches config's expected chemistry (e.g. v2 data but config set for v3).
- Symmetric equal-length mates (e.g. 101/101, 151/151) are **NOT auto-refused** — the barcode
  SEQUENCE probe (Channel 2) runs first. If the leading 16 bp match a 10x whitelist, it's an
  over-sequenced-R1 10x run and is processed (`barcode_read_length: 0`); only if the probe finds
  no whitelist match does the guard **refuse, exit 1, no sheet written** (genuine bulk /
  full-length / Smart-seq → use plain STAR alignReads + featureCounts, not this pipeline). The
  refusal message is sharpened by SRA `LibrarySource` (bulk vs plate-based single cell).
- SRA-only rows can't be FASTQ-peeked, but their chemistry is filled from the GEO
  `data_processing` software (Channel 3): non-droplet software → refuse; `cellranger-arc` → arc.
  If the record names no known aligner, a note is printed and it's verified after first download.

Flags: `--chem v2` (config is set for v2, adjust the warn threshold), `--no-chem-check` (skip the
guard), `--multiome` (assert `chemistry: arc` on every row), `--no-geo-software` (skip the GEO
software lookup — offline/fallback), `--whitelist-dir DIR` (10x whitelist dir for the barcode probe).

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
sbatch run_slurm.sh                 # whole sheet, all samples (uses config/config.yaml)
./run_slurm.sh -n                   # dry-run on the login node first (recommended)
```

`run_slurm.sh` now always prints, loudly, which configfile it resolved and that config's
`workdir`/`samples_tsv` near the top of its log — check that against what you intended before
trusting the run. **The moment you keep more than one study's config side by side in `config/`**
(common once you're not starting from a fresh clone), pass `--configfile` explicitly on every
launch instead of relying on the default:
```bash
sbatch run_slurm.sh --configfile config/config.<STUDY>.yaml
```
A bare invocation with multiple configs present will still run — silently defaulting to
`config/config.yaml` — and that default may be a *different, fully valid* study. Nothing inside the
pipeline can catch this: the "genomeVersion" check in `resolve_star_index` is STAR's own
index-format compatibility (binary vs. index), not species, and configs carry no
`species`/`organism` field. Each config is fully self-contained (own `fasta`/`gtf`/`star_index`/
`workdir`/`samples_tsv`), so the wrong one runs cleanly to completion with zero errors — it just
processes the wrong study.

**Smoke-test ONE sample end-to-end** (5 jobs: check_versions → download_fastq →
fastqc → starsolo → to_h5ad) by targeting that sample's outputs directly. Put targets *before*
`--configfile` — anything after it is consumed as another config file:
```bash
sbatch run_slurm.sh \
  results/h5ad/<SAMPLE>.h5ad \
  results/qc/fastqc/<SRR>_R1_fastqc.zip \
  results/qc/fastqc/<SRR>_R2_fastqc.zip \
  --configfile config/config.<STUDY>.yaml
```
Do **not** add `qc_gate` / `merge` / `multiqc` to a single-sample target — those
rules `expand` over ALL samples and won't resolve.

## Phase 3 — Watch progress (don't confuse "running" with "hung")

```bash
squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"      # ctl + per-rule jobs
LOG=$(ls -t "$PIPE"/.snakemake/log/*.snakemake.log | head -1); tail -20 "$LOG"
```
Per-rule slurm logs: `$PIPE/.snakemake/slurm_logs/rule_<name>/<wildcards>/<jobid>.log`.

**Before trusting a run, don't stop at the "Config file(s):" / launch banner line — confirm it's
actually processing the study you intended:**
```bash
grep -aoE "SRR[0-9]+" "$LOG" | sort -u | head           # accessions actually in flight
grep -Fc "<one of those SRRs>" config/samples.<STUDY>.tsv   # do they belong to this study?
readlink -f <workdir>/star_index                          # resolves to the right reference?
ls <workdir>/fastq/ | tail                                 # outputs landing in the right workdir?
```
A config filename alone isn't proof it's the right *study* — only that some valid config loaded.

The account-guessing lines (`No SLURM account given, trying to guess` /
`sacct: invalid option`) are **cosmetic** — the executor proceeds without an account.
They are NOT the hang. A download job legitimately takes 20–40 min (prefetch +
fasterq-dump + gzip of a multi-GB run); growing `results/fastq/*.fastq.gz` = healthy.

**Check node-local `/scratch` usage on any node running SRA downloads, especially during a long
multi-sample run.** `download_fastq`'s SRA branch stages on node-local disk, and orphaned
`sra_<jobid>` directories from past killed/crashed jobs (`trap ... EXIT` doesn't fire on SIGKILL)
accumulate there over days until a node fills up — a download that starts writing into an
almost-full disk produces a **silently truncated FASTQ**, not a clean error. This has caused a
whole day of hard-to-diagnose corruption on a real run: nodes at 90%+ full from dozens of orphaned
staging dirs going back days, and every truncated file traced back to landing on one of them.
```bash
squeue -u $USER -h -o "%N" | sort -u                              # nodes your jobs are on
for n in <nodes>; do ssh "$n" 'df -h /scratch'; done               # check each
ssh <node> 'du -sh /scratch/sra_* 2>/dev/null | sort -rh | head'   # find the big ones
```
Cross-check `sra_<jobid>` directory names against `squeue -u $USER -h -o "%i"` — any not currently
running is orphaned and safe to delete. The rule has a preflight check that refuses to stage if
free space is under 100GB (fails loudly instead of truncating), but that doesn't clean up existing
orphaned space — do the sweep above periodically on any node this pipeline uses heavily.

## Phase 4 — Verify a completed run

Two read-only scripts cover this. Both take the **workdir**, launch nothing, work
part-way through a run, and exit non-zero when something is wrong — so they compose
into a shell gate. Prefer them over ad-hoc `cat`/`grep` of Solo.out: they emit one
compact table instead of dumping files.

```bash
# 1. QC metrics for every aligned sample (imports qc_gate.py, so these numbers are
#    the same code path the in-DAG gate enforces — they cannot drift apart).
python workflow/scripts/inspect_qc.py <WORKDIR> --config config/config.yaml

# 2. Velocity layers present? Names the CAUSE, not just the symptom.
python workflow/scripts/check_layers.py <WORKDIR>
```

`inspect_qc.py` reads thresholds and `feature` from `--config`, overridable with
`--min-valid-barcodes` / `--min-reads-mapped-gene` / `--min-estimated-cells`.
A near-zero valid-barcode fraction means the chemistry/whitelist is wrong (e.g. v2
data run with a v3 config) — go back to Phase 1 and fix `config/config.yaml`.

`check_layers.py` distinguishes the two ways layers go missing, because they need
fixes in different files:

| status | meaning | fix |
|---|---|---|
| `OK` | all layers present, non-empty | — |
| `NOT_ATTACHED` | Solo.out has counts, h5ad doesn't | `starsolo_to_h5ad.py` reads `Velocyto/<sub>` following the Gene matrix (`filtered`), but STARsolo only ever writes `Velocyto/raw/` |
| `UPSTREAM_EMPTY` | Velocyto `.mtx` exist but all `nnz=0` | STARsolo's Velocyto counter consumes the **Gene** pass. `--soloFeatures GeneFull Velocyto` yields zeros *silently*; it needs `--soloFeatures Gene GeneFull Velocyto` |
| `NO_VELOCYTO` | no Velocyto dir | Velocyto not requested in `--soloFeatures` |
| `MISSING` | no h5ad yet | sample hasn't reached `to_h5ad` |

`UPSTREAM_EMPTY` is the one to watch on **snRNA/Multiome runs**, since that is exactly
the `feature: GeneFull` path — an all-zero Velocyto is easy to miss because STAR exits 0
and the h5ad is written normally, just with no layers.

It reads only HDF5 group structure and MatrixMarket headers, so cost is flat in file
size — safe to run over a whole study.

The resulting `.h5ad` (count matrix in `X`, velocyto layers) is ready for
downstream scanpy analysis (QC, clustering, UMAP, DE).

### Checking many large FASTQs for corruption (mod4 line-count sweep)

A truncated FASTQ (record boundary cut mid-write, e.g. by a job hitting its time
limit or a disk filling up) usually still passes `gzip -t` — the gzip container
closes cleanly either way. The only cheap, reliable check is: decompressed line
count must be a multiple of 4, and R1/R2 line counts must match.

```bash
zcat "$f" | wc -l   # mod 4 == 0, and R1 count == R2 count for the pair
```

Running this serially (or in a single sandboxed process) does not scale: each
~800M-line file takes ~9 minutes, so a 70-file study takes 10+ hours serially.
Submit it as its own SLURM job and fan out with `xargs -P`:

```bash
printf '%s\n' "${files[@]}" | xargs -P 8 -I{} bash -c \
  'n=$(zcat "{}" | wc -l); echo "{} lines=$n mod4=$((n % 4))"'
```

8-way parallel brings a 70-file study sweep down to well under an hour.

### Reading the pipeline's own code cheaply

`workflow/scripts/symctx.py` prints just the slice of Python you need, so reviewing a
script doesn't pull the whole file into context. Stdlib `ast` only.

```bash
python workflow/scripts/symctx.py outline workflow/scripts/qc_gate.py   # signatures, bodies folded
python workflow/scripts/symctx.py show   workflow/scripts/qc_gate.py evaluate
python workflow/scripts/symctx.py find   workflow/scripts build_matrices
```

## Phase 5 — Batch completeness & the agent loop (many studies)

Snakemake calls a sample done as soon as its `.h5ad` merely *exists*. For a real
"is every dataset finished?" check across many studies, use `reconcile.py` +
`run_batch.sh`. This layer only ever hands re-runnable work back to Snakemake — it
never downloads or aligns itself.

**"Done" contract** — a sample is done iff **all three**: `<workdir>/h5ad/<sample>.h5ad`
exists, `<sample>` is in `<workdir>/qc/qc_pass.txt`, and that h5ad has layers
`spliced`/`unspliced`/`ambiguous`. Everything else gets a **category** + **action**:

| category | action | why |
|---|---|---|
| `missing` | **rerun** | not produced yet — hand back to Snakemake |
| `corrupt` | **rerun** | unreadable/truncated h5ad — quarantined first, then regenerated |
| `read_qc_fail` | **flag** | raw reads failed read-QC — re-aligning the same reads can't help |
| `qc_fail` | **flag** | failed the alignment gate — reason in `qc/qc_gate.tsv` |
| `no_layers` | **flag** | QC passed but a velocity layer is absent (silent mis-wire) |

**Auto-rerun ONLY the transient categories (`missing`/`corrupt`); always flag the
rest for the human** — never auto-retry a genuine failure.

```bash
# activate an env with h5py + pyyaml (snakemake's is fine)
# One workdir (read-only; launches nothing). exit 0 = done, 1 = work remains, 2 = error.
python workflow/scripts/reconcile.py --samples config/samples.tsv \
       --workdir "$PIPE"/results/<ACC> --accession <ACC> --feature Gene \
       --report qc_reconcile.tsv --ledger results/successful_samples.tsv

# Every study in the manifest (copy config/manifest.example.tsv -> manifest.tsv first):
python workflow/scripts/reconcile.py --manifest config/manifest.tsv --base-dir "$PIPE" \
       --json recon.json --ledger results/successful_samples.tsv

# ONE batch cycle (reconcile -> record logbook -> rerun missing + quarantine/rerun
# corrupt -> flag the rest). NOT a daemon; call it repeatedly. DRY_RUN=1 previews.
SCRNASEQ_PY=$(command -v python) ./run_batch.sh          # add DRY_RUN=1 (no sbatch)
```

**Two human gates** (the loop is thin between them):
- **Gate 1 (before the loop):** approve `config/manifest.tsv`, and make sure each
  study has a sample sheet — `prepare_runs.py` (Phase 1) is a human step because of
  the chemistry guard, so the loop never guesses chemistry. A study with no sheet is
  flagged "samples sheet not found", not run.
- **Loop:** call `run_batch.sh` repeatedly while it exits 1 **and** launched something
  this cycle. A study's `missing`/`corrupt` samples are relaunched at most
  `MAX_RETRIES` (default 2, tracked in `.batch_state.json`), then escalated.
- **Gate 2 (after):** review the `--report` TSV + `qc_gate.tsv` + `read_qc.tsv` and
  decide on flagged studies (fix chemistry & re-prepare, accept, or drop). Successful
  samples are already recorded in the logbook.

**The logbook** — `results/successful_samples.tsv` is the durable, append-only,
idempotent record of what mapped: one row per done `(accession, sample)` with its
metrics, chemistry, feature, SRRs, and read-QC status. Check it to answer "what
finished?" without re-scanning the trees.

## Troubleshooting (all fixes stay inside `$PIPE`)

**In plain terms, corrupted-looking FASTQs usually come from one of three causes,
none of which throw a clean error:** (1) the compute node's local disk was nearly
full and a write got cut off mid-way; (2) compression was too slow (single-threaded
`gzip` on a huge file) and the job's time limit killed it before it finished
writing; (3) a long-running controller loaded the pipeline's rule code once at
startup, so a bug fix you made on disk never reached jobs it was still dispatching
— restarting the controller is what makes a fix "live." All three leave a file that
looks fine (exists, valid gzip) but is truncated — the reliable check is that the
decompressed line count must be a multiple of 4, not `gzip -t`.

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
| Same corruption bug recurs on new SRRs after you already fixed and re-verified the Snakefile | a long-running controller process loads the Snakefile into memory once at startup; editing the file on disk does NOT change what jobs it's still dispatching | Kill and restart the controller after any Snakefile edit meant to change in-flight rule behavior — a fix on disk is not live until the controller that's actually submitting jobs is restarted. |
| A redownloaded/recompressed FASTQ looks corrupted even though the recovery step "succeeded" | compressing straight into the final output path is not atomic; a reader (e.g. Snakemake) can observe a partially-written file mid-write | Write to `<output>.tmp` then `mv` into place — `mv` on the same filesystem is atomic. Already done in `download_fastq`'s SRA branch and in any standalone recovery script that recompresses FASTQs. |
| Only some files from a batch downloaded during a disk-full/contention window failed loudly downstream | corruption during a resource-constrained window is not limited to the files that happen to error out immediately — a truncated file with valid record boundaries can pass silently into `starsolo` and just produce quietly-wrong counts | When you find ANY corrupt file from a given time window, sweep every file downloaded in that same window, not just the ones that already failed. In one incident, 2 of 6 SRRs from the same window failed loudly, but a full sweep found all 6 were corrupt. |
| Several large SRA downloads running at once take much longer than expected, occasionally past their time limit, even though each one's own extraction/compression is otherwise fine | multiple concurrent `download_fastq` jobs landed on the same node's local disk, saturating its I/O (`iostat -x` showing `%util` near 100%) — a distinct cause from CPU contention, and `ssh <node> nproc` can mislead here since it reflects the ssh session's own cgroup quota, not the node's real CPU count (use `scontrol show node <name>` instead) | `download_fastq` declares `resources: download_io=1`; pass `--resources download_io=4` (run_slurm.sh does this by default) to cap concurrent downloads cluster-wide, regardless of which node each lands on, so they can't stack up on one disk. |

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
- `workflow/scripts/fetch_geo_supp.py` — list/download a GEO record's **supplementary**
  files (the authors' own `.h5ad`, `.rds`, `.loom`, mtx bundles). Standalone, NOT in the DAG:
  ```bash
  python workflow/scripts/fetch_geo_supp.py GSE218661                       # list only
  python workflow/scripts/fetch_geo_supp.py GSE218661 -o DIR --download \
      --pattern '*.h5ad.gz' --gunzip
  ```
  Resumes (`curl -C -`), verifies Content-Length, skips files already complete.
  Exit 0 ok / 1 transfer failed / 2 accession or pattern matched nothing. An author h5ad is
  **not** equivalent to pipeline output — different reference and filtering, and normally no
  spliced/unspliced/ambiguous layers. Use it for a quick look; reprocess (Phases 1–2) when the
  matrix must be comparable across studies or carry velocity layers.
- `workflow/scripts/starsolo_to_h5ad.py` — STARsolo Gene/GeneFull + Velocyto feature → single .h5ad (`--feature Gene|GeneFull`).

> **"Velocyto" here = STARsolo's `--soloFeatures Velocyto` mode, NOT the separate
> velocyto.py tool.** STARsolo itself classifies each read as spliced / unspliced /
> ambiguous (exon vs intron overlap) and writes three matrices under
> `Solo.out/Velocyto/`. No separate velocyto tool runs. `starsolo_to_h5ad.py` loads
> those into `adata.layers['spliced'/'unspliced'/'ambiguous']`, alongside the
> standard `Gene/` matrix as `adata.X`.
- `workflow/scripts/read_qc.py` — role-aware raw-read QC from FastQC (rule `read_qc`); writes `qc/read_qc.tsv`.
- `workflow/scripts/reconcile.py` — completeness checker (done-contract, categories/actions); single-workdir or `--manifest`.
- `workflow/Snakefile` — rules; `download_fastq` branches on `source`; `read_qc` reporting rule.
- `profiles/slurm/config.yaml` — SLURM profile + per-rule resources + sacct/squeue fix + `keep-going`/`retries`.
- `config/config.yaml` — chemistry (umi_len, whitelist), `feature` (Gene/GeneFull), `read_qc` thresholds, reference paths.
- `config/manifest.tsv` — one row per study for the batch loop (copy from `manifest.example.tsv`).
- `run_slurm.sh` — controller launcher (venv + PATH + `snakemake --profile`).
- `run_batch.sh` — batch agent loop (reconcile → rerun/flag); `snakemake_status.sh` watches jobs.
- `results/successful_samples.tsv` — the success logbook (append-only record of what mapped).
