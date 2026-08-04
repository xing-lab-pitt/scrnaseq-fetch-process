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
# Adjust the #SBATCH lines below (partition, time) to your cluster.
# =============================================================================
#SBATCH -J snakemake_ctl
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 7-00:00:00
#SBATCH --output=snakemake_ctl.%j.log

set -euo pipefail

# Run from the repo root (this script's directory), wherever it was cloned.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

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
