#!/usr/bin/env bash
# snakemake_status.sh — DB-free status of a Snakemake + SLURM run.
# Needs only squeue + the log files; does NOT use sacct or seff.
#
# Classifies the run as RUNNING / SUCCESS / FAILED / INCOMPLETE from the
# controller log's terminal markers, and on failure prints the failing rule(s)
# and tails their real stderr (the SLURM per-rule logs).
#
# Usage:  PIPE=/path/to/repo ./snakemake_status.sh      (PIPE defaults to $PWD)
set -uo pipefail

PIPE="${PIPE:-$(pwd)}"                    # Snakemake workdir (holds .snakemake/)
CTL_NAME="${CTL_NAME:-snakemake_ctl}"     # controller job name (from #SBATCH -J)
LOGDIR="$PIPE/.snakemake/log"
SLURMLOGS="$PIPE/.snakemake/slurm_logs"

LOG="$(ls -t "$LOGDIR"/*.snakemake.log 2>/dev/null | head -1)"
[ -z "${LOG:-}" ] && { echo "no controller log under $LOGDIR"; exit 2; }
echo "controller log: $LOG"

# 1) Controller still alive? (liveness — no accounting DB needed)
if squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -qx "$CTL_NAME"; then
  echo "STATE: RUNNING"
  grep -oE '[0-9]+ of [0-9]+ steps \([0-9]+%\) done' "$LOG" | tail -1
  tail -3 "$LOG"
  exit 0
fi

# 2) Controller gone — classify from terminal markers in the log
if grep -qE 'Exiting because a job execution failed|Error in rule ' "$LOG"; then
  echo "STATE: FAILED"
elif grep -qE '\(100%\) done|Nothing to be done' "$LOG"; then
  echo "STATE: SUCCESS"; exit 0
else
  echo "STATE: INCOMPLETE — controller ended with no success/failure marker."
  echo "  Usually the CONTROLLER job itself was killed (walltime/OOM). Check its"
  echo "  own SLURM output (snakemake_ctl.*.log) and raise -t / --mem in run_slurm.sh."
fi

# 3) Which rule(s) failed, and why (real error is in the rule's own logs)
echo; echo "=== failing rule block(s) from controller log ==="
grep -A 12 'Error in rule ' "$LOG" | sed 's/^/  /'

RULES="$(grep -oE 'Error in rule [A-Za-z0-9_]+' "$LOG" | awk '{print $4}' | sort -u)"
echo; echo "failed rule(s): ${RULES:-<none parsed>}"

# 4) Tail the SLURM per-rule stderr for each failed rule (find handles nesting)
for r in $RULES; do
  echo; echo "=== newest SLURM log for rule '$r' ==="
  f="$(find "$SLURMLOGS/rule_$r" -name '*.log' 2>/dev/null -print0 \
        | xargs -0 ls -t 2>/dev/null | head -1)"
  if [ -n "${f:-}" ]; then echo "$f"; tail -20 "$f"
  else echo "  (no SLURM log under $SLURMLOGS/rule_$r/ — check the rule's declared log: path above)"; fi
done
