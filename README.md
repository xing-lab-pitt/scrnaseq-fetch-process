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

## Features

- **Any accession type** — `prepare_runs.py` resolves GSE / SRP / PRJNA / SRX /
  SRR / GSM into a sample sheet, scoping correctly (a study gives all its runs;
  an experiment/run gives just its own).
- **ENA-first download with SRA fallback** — per run, prefers ENA's pre-split,
  md5-checked FASTQs over HTTPS (`curl`); falls back to NCBI `prefetch` +
  `fasterq-dump --include-technical` for runs ENA hasn't mirrored. Chosen
  automatically per run (the `source` column); no global switch.
- **Chemistry guard** — before writing the sheet, peeks the barcode read length
  and refuses non-10x-droplet data (bulk / Smart-seq), warns on a v2/v3 mismatch.
- **Deterministic QC gate** — pass/fail thresholds live in `config.yaml`
  (version-controlled), applied to STARsolo's `Summary.csv`.
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
#    SRX/SRR/GSM/SRP/PRJNA also work. Add --chem v2 if your config is v2.

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

Defaults are 10x Chromium 3′ **v3**. Supported 3′ chemistries:

| Chemistry | Barcode read | CB / UMI | `umi_len` | Whitelist |
|---|---|---|---|---|
| v2 | 26 bp | 16 / 10 | `10` | `737K-august-2016.txt` |
| **v3** / v3.1 (default) | 28 bp | 16 / 12 | `12` | `3M-february-2018.txt` |
| **v4** / GEM-X | 28 bp | 16 / 12 | `12` | `3M-3pgex-may-2023.txt` |

The barcode+UMI read and the cDNA read are separate FASTQs, so STARsolo runs in
the standard two-file mode — the barcode read is never clipped, and there is no
`clip5p`/`barcode_mate` setting to tune.

**v3 vs v4:** identical read geometry — the *only* difference is the whitelist. To
run v4 data, keep the v3 geometry and just point `chemistry.whitelist` at
`3M-3pgex-may-2023.txt`. Because read length can't distinguish them, the chemistry
guard reports both as "v3"; pass `--chem v4` to `prepare_runs.py` and it prints a
reminder to confirm the whitelist rather than a false mismatch warning.

Wrong chemistry (or the wrong whitelist for the right geometry) silently yields
near-empty matrices — the QC gate's `min_valid_barcodes` catches this.

## Repository layout

```
config/
  config.example.yaml     # template — copy to config.yaml and edit
  samples.example.tsv     # sample-sheet schema example
workflow/
  Snakefile               # the DAG
  scripts/                # prepare_runs, ncbi_utils, starsolo_to_h5ad, qc_gate, ...
  envs/                   # conda env specs (scanpy.yaml, tools.yaml)
profiles/slurm/           # SLURM profile (edit partitions for your cluster)
run_slurm.sh              # SLURM launcher (paths from env vars / script location)
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
