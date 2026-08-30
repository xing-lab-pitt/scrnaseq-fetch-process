#!/bin/bash
# =============================================================================
# Launch the scRNA-seq pipeline on SLURM (Snakemake submits each rule as a job).
#
#   sbatch run_slurm.sh                 # run the controller itself as a job
#   ./run_slurm.sh -n                   # or dry-run directly on a login node
#
# Extra args are passed through to snakemake (e.g. -n for a dry run, or a
# specific target path to run one rule).
#
# WHERE THE STAR INDEX IS STORED (custom directory):
#   The index location is a CONFIG setting, not a flag here. In config/config.yaml
#   set reference.star_index to the directory you want the built index kept in:
#     reference:
#       star_index: "/path/to/refs/STAR_rebuild"   # writable dir to build+keep it in
#   rule resolve_star_index then builds a fresh index INTO that dir (once) and
#   symlinks <workdir>/star_index -> it; later runs reuse it (no rebuild). If it
#   already holds a genomeVersion-compatible index, it is symlinked as-is. Leave
#   star_index "" to build into the (ephemeral) workdir instead. To force one build
#   without a full run, target that path directly:
#     sbatch run_slurm.sh <workdir>/star_index      # <workdir> = config workdir
#
# CONFIGURE FOR YOUR SITE via environment variables (all optional):
#   SCRNASEQ_VENV       path to a Python venv to activate (snakemake, executor
#                       plugin, multiqc, scanpy). If unset, snakemake must already
#                       be on PATH (e.g. an activated conda env or module).
#   SCRNASEQ_EXTRA_PATH colon-separated dirs to prepend to PATH before running,
#                       for binary tools not in the venv — e.g. a FastQC dir or
#                       an sra-tools bin dir:
#                         export SCRNASEQ_EXTRA_PATH=/opt/FastQC:/opt/sratoolkit/bin
#   SCRNASEQ_PROFILE    Snakemake profile dir (default: profiles/slurm)
#
# Adjust the #SBATCH lines below (partition, time) to your cluster. Keep
# --cpus-per-task at 8 or more: the controller's SLURM executor plugin spawns one
# polling thread per concurrently-tracked job (see profiles/slurm's `jobs:`), and
# on cgroup-managed clusters a too-small CPU allocation can cap the controller's
# thread budget, crashing it with "RuntimeError: can't start new thread" under
# sustained concurrency/retries. 1 CPU was not enough at `jobs: 24`.
# =============================================================================
#SBATCH -J snakemake_ctl
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH -t 7-00:00:00
#SBATCH --output=snakemake_ctl.%j.log

set -euo pipefail

# Run from the repo root (this script's directory), wherever it was cloned.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

# Visibility guard: a bare invocation with no --configfile silently falls back to
# Snakemake's default config/config.yaml. That's fine for a single-study clone (see
# README quickstart), but the moment more than one study's config lives side by
# side in config/, that silent default can run the WRONG study with zero errors —
# it's a fully valid, self-consistent config, just not the intended one. Rather
# than block the documented single-config workflow, always resolve and print
# LOUDLY which config is about to run, so a wrong-or-defaulted choice is obvious
# at a glance in the log instead of invisible.
CONFIGFILE=""
prev=""
for arg in "$@"; do
    [ "$prev" = "--configfile" ] && CONFIGFILE="$arg"
    prev="$arg"
done
: "${CONFIGFILE:=config/config.yaml}"
if [ ! -f "$CONFIGFILE" ]; then
    echo "ERROR: configfile '$CONFIGFILE' not found." >&2
    echo "First run: cp config/config.example.yaml config/config.yaml (see README)," >&2
    echo "or pass --configfile config/config.<STUDY>.yaml explicitly." >&2
    exit 1
fi
echo "=============================================================="
echo "LAUNCHING WITH configfile: $CONFIGFILE"
grep -E '^workdir:|^samples_tsv:' "$CONFIGFILE" | sed 's/^/  /'
echo "If more than one study's config lives in config/, pass --configfile explicitly"
echo "next time — a bare invocation always defaults to config/config.yaml."
echo "=============================================================="

PROFILE="${SCRNASEQ_PROFILE:-profiles/slurm}"

# Activate a venv if one was provided.
if [[ -n "${SCRNASEQ_VENV:-}" ]]; then
    # shellcheck disable=SC1091
    source "${SCRNASEQ_VENV}/bin/activate"
fi

# Prepend any site-specific tool dirs (FastQC, sra-tools, ...) so they propagate
# to every SLURM job (SLURM exports the submit environment by default).
if [[ -n "${SCRNASEQ_EXTRA_PATH:-}" ]]; then
    export PATH="${SCRNASEQ_EXTRA_PATH}:$PATH"
fi

# Sanity: fail early if a required binary is missing. sra-tools (prefetch) is
# only needed for source=sra runs, so it is checked but non-fatal.
for t in snakemake STAR samtools fastqc multiqc; do
    command -v "$t" >/dev/null || { echo "MISSING on PATH: $t (see run_slurm.sh header / README)" >&2; exit 1; }
done
command -v prefetch >/dev/null || echo "NOTE: prefetch not on PATH (only needed for source=sra runs)" >&2

exec snakemake --profile "$PROFILE" "$@"
