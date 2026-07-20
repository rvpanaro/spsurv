#!/usr/bin/env bash
# Wait for an in-flight degree study, then run remaining paper workflow steps.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
WAIT_PID="${1:-}"
if [[ -n "$WAIT_PID" ]] && kill -0 "$WAIT_PID" 2>/dev/null; then
  echo "Waiting for PID $WAIT_PID (degree study)..."
  while kill -0 "$WAIT_PID" 2>/dev/null; do
    sleep 60
  done
  echo "Degree study process finished."
fi
export SPSURV_ROOT="$ROOT"
export SPSURV_MC_REPLICATES=1000
export SPSURV_MC_NSIZES=50,100
export SPSURV_DEGREE_MC_REPS=1000
export SPSURV_DEGREE_MC_NSIZES=50,100
unset SPSURV_MC_FAST_BAYES || true
unset SPSURV_MC_SEQUENTIAL || true
LOG="paper/logs/workflow-mc-$(date +%Y%m%d-%H%M%S).log"
echo "Starting MC + render + PDF; log=$LOG"
nohup Rscript paper/run-full-workflow.R --from-step=mc_regression >>"$LOG" 2>&1 &
echo $! > paper/logs/workflow-mc.pid
echo "PID=$(cat paper/logs/workflow-mc.pid)"
