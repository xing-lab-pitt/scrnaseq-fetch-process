#!/usr/bin/env bash
# =============================================================================
# run_batch.sh — ONE reconcile -> decide -> (re)launch cycle over the manifest.
#
# Keeps mechanics in a script; keeps judgment with the human. Not a daemon: it
# performs a single pass and exits. The thin agent loop calls it repeatedly,
# watching progress with snakemake_status.sh between calls, until the reconciler
# reports COMPLETE or only human-flagged failures (qc_fail / no_layers) remain.
#
#   ./run_batch.sh [MANIFEST]        MANIFEST defaults to config/manifest.tsv
#   DRY_RUN=1 ./run_batch.sh         reconcile + show what WOULD launch, no sbatch
#   MAX_RETRIES=2 ./run_batch.sh     cap auto-relaunches per study (default 2)
#   LEDGER=results/successful_samples.tsv   success logbook (append-only)
#
# Decision per study (from reconcile.py --manifest categories + `action`):
#   missing   -> RERUN: transient / not-run-yet. (Re)launch this study's Snakemake
#                via run_slurm.sh --configfile (idempotent; Snakemake skips done
#                work). Capped by MAX_RETRIES.
#   corrupt   -> RERUN: the .h5ad exists but is unreadable (a truncated / killed
#                write). QUARANTINE it (move to <workdir>/.quarantine/<stamp>/) so
#                Snakemake regenerates it, then relaunch. Also capped.
#   read_qc_fail / qc_fail / no_layers -> FLAG: a rerun of the same inputs can't
#                fix it (bad reads / failed gate / mis-wire). NEVER auto-retried;
#                listed for human review.
#
# Every cycle also records the samples that are DONE into the success logbook
# ($LEDGER, append-only + idempotent) via reconcile.py --ledger.
#
# Exit code: 0 when the manifest to-do list is empty (reconciler COMPLETE),
# 1 when work remains (launched and/or flagged), 2 on setup error.
# =============================================================================
set -uo pipefail

PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPE"

MANIFEST="${1:-config/manifest.tsv}"
MAX_RETRIES="${MAX_RETRIES:-2}"
DRY_RUN="${DRY_RUN:-0}"
BASE_CONFIG="${BASE_CONFIG:-config/config.yaml}"
STATE="${STATE:-.batch_state.json}"          # per-study relaunch counters
STUDIES_DIR="config/studies"
LEDGER="${LEDGER:-results/successful_samples.tsv}"   # success logbook (append-only)

PY="${SCRNASEQ_PY:-python}"
command -v "$PY" >/dev/null || { echo "ERROR: python ('$PY') not on PATH — activate the venv" >&2; exit 2; }
[ -f "$MANIFEST" ]     || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 2; }
[ -f "$BASE_CONFIG" ]  || { echo "ERROR: base config not found: $BASE_CONFIG" >&2; exit 2; }
mkdir -p "$STUDIES_DIR"

RECON_JSON="$(mktemp)"
trap 'rm -f "$RECON_JSON"' EXIT

# ---- 1. Reconcile the whole manifest -------------------------------------- #
# Also records DONE samples into the success logbook (append-only + idempotent).
echo "### Reconciling manifest: $MANIFEST"
"$PY" workflow/scripts/reconcile.py --manifest "$MANIFEST" --base-dir "$PIPE" \
      --json "$RECON_JSON" --ledger "$LEDGER"
RECON_EXIT=$?
if [ "$RECON_EXIT" -eq 2 ]; then
    echo "ERROR: reconcile failed (see above)" >&2; exit 2
fi
if [ "$RECON_EXIT" -eq 0 ]; then
    echo "### Manifest COMPLETE — nothing to launch."
    exit 0
fi

# ---- 2. Per study: decide launch vs flag ---------------------------------- #
# Emit one TAB line per study, parsed from the reconcile JSON (feature comes from
# the JSON now that reconcile carries it):
#   accession workdir samples_tsv feature n_missing n_corrupt n_readqcfail n_qcfail n_nolayers error
mapfile -t STUDY_LINES < <(
  "$PY" - "$RECON_JSON" <<'PY'
import json, sys
recon = json.load(open(sys.argv[1]))
for s in recon.get("studies", []):
    c = s["counts"]
    print("\t".join([
        s["accession"], s["workdir"], s["samples_tsv"], s.get("feature", ""),
        str(c["missing"]), str(c["corrupt"]), str(c["read_qc_fail"]),
        str(c["qc_fail"]), str(c["no_layers"]),
        (s.get("error") or "").replace("\t", " "),
    ]))
PY
)

# ---- state helpers (relaunch counters) ------------------------------------ #
get_count() { # get_count <accession>
  [ -f "$STATE" ] || { echo 0; return; }
  "$PY" - "$STATE" "$1" <<'PY'
import json, sys
try: d = json.load(open(sys.argv[1]))
except Exception: d = {}
print(int(d.get(sys.argv[2], 0)))
PY
}
bump_count() { # bump_count <accession>
  "$PY" - "$STATE" "$1" <<'PY'
import json, sys
p, acc = sys.argv[1], sys.argv[2]
try: d = json.load(open(p))
except Exception: d = {}
d[acc] = int(d.get(acc, 0)) + 1
json.dump(d, open(p, "w"), indent=2)
PY
}

# quarantine_corrupt <accession> <stamp> — move every `corrupt` sample's unreadable
# .h5ad for this study into <workdir>/.quarantine/<stamp>/ so Snakemake regenerates
# it on the next launch (we never delete; the bad file is preserved for inspection).
# Prints one line per moved file. Reads the corrupt rows straight from RECON_JSON.
quarantine_corrupt() {
  "$PY" - "$RECON_JSON" "$1" "$2" <<'PY'
import json, shutil, sys
from pathlib import Path
recon, acc, stamp = json.load(open(sys.argv[1])), sys.argv[2], sys.argv[3]
for s in recon.get("studies", []):
    if s["accession"] != acc:
        continue
    qdir = Path(s["workdir"]) / ".quarantine" / stamp
    for r in s["rows"]:
        if r.get("category") != "corrupt":
            continue
        src = Path(r["h5ad"])
        if not src.exists():
            continue
        qdir.mkdir(parents=True, exist_ok=True)
        dst = qdir / src.name
        shutil.move(str(src), str(dst))
        print(f"    quarantined {src} -> {dst}")
PY
}

LAUNCHED=(); FLAGGED=(); CAPPED=()

for line in "${STUDY_LINES[@]}"; do
  IFS=$'\t' read -r acc workdir samples_tsv feature \
      n_missing n_corrupt n_readqcfail n_qcfail n_nolayers err <<<"$line"

  if [ -n "$err" ]; then
    FLAGGED+=("$acc: setup error — $err")
    continue
  fi
  # FLAG categories: a rerun of the same inputs can't fix these — human review.
  if [ "$n_readqcfail" -gt 0 ] || [ "$n_qcfail" -gt 0 ] || [ "$n_nolayers" -gt 0 ]; then
    FLAGGED+=("$acc: read_qc_fail=$n_readqcfail qc_fail=$n_qcfail no_layers=$n_nolayers (needs human — not auto-retried)")
  fi

  # RERUN categories: missing (not produced yet) + corrupt (truncated write).
  n_relaunch=$((n_missing + n_corrupt))
  if [ "$n_relaunch" -eq 0 ]; then
    continue                              # nothing re-runnable for this study
  fi

  tries="$(get_count "$acc")"
  if [ "$tries" -ge "$MAX_RETRIES" ]; then
    CAPPED+=("$acc: $n_missing missing + $n_corrupt corrupt after $tries relaunch(es) — escalating to human")
    continue
  fi

  # --- generate the per-study config (override workdir/samples_tsv/feature) --
  study_cfg="$STUDIES_DIR/${acc}.yaml"
  "$PY" - "$BASE_CONFIG" "$study_cfg" "$workdir" "$samples_tsv" "$feature" <<'PY'
import sys, yaml
base, out, workdir, samples_tsv, feature = sys.argv[1:6]
cfg = yaml.safe_load(open(base))
cfg["workdir"] = workdir
cfg["samples_tsv"] = samples_tsv
if feature:
    cfg["feature"] = feature
yaml.safe_dump(cfg, open(out, "w"), sort_keys=False)
print(f"wrote {out}")
PY

  if [ "$DRY_RUN" = "1" ]; then
    if [ "$n_corrupt" -gt 0 ]; then
      echo "DRY_RUN: would quarantine $n_corrupt corrupt .h5ad for $acc into $workdir/.quarantine/<stamp>/"
    fi
    echo "DRY_RUN: would launch $acc: sbatch run_slurm.sh --configfile $study_cfg"
    LAUNCHED+=("$acc (dry-run)")
    continue
  fi

  # Move corrupt outputs aside FIRST so Snakemake regenerates them on this launch.
  if [ "$n_corrupt" -gt 0 ]; then
    stamp="$(date +%Y%m%dT%H%M%S)"
    echo "### Quarantining $n_corrupt corrupt .h5ad for $acc (stamp $stamp)"
    quarantine_corrupt "$acc" "$stamp"
  fi

  echo "### Launching $acc ($n_missing missing + $n_corrupt corrupt, attempt $((tries + 1))/$MAX_RETRIES)"
  if sbatch run_slurm.sh --configfile "$study_cfg"; then
    bump_count "$acc"
    LAUNCHED+=("$acc")
  else
    FLAGGED+=("$acc: sbatch submission failed")
  fi
done

# ---- 3. Report ------------------------------------------------------------ #
echo
echo "=== batch cycle summary ==="
printf 'launched : %s\n' "${LAUNCHED[*]:-none}"
printf 'flagged  : %s\n' "${FLAGGED[*]:-none}"
printf 'capped   : %s\n' "${CAPPED[*]:-none}"
echo "Manifest still INCOMPLETE — rerun this cycle after jobs finish"
echo "(watch: PIPE=$PIPE ./snakemake_status.sh ; or squeue -u \$USER)"
exit 1
