#!/usr/bin/env bash
# Production overnight run: preflight, simulation (parallel), figures, TeX tables.
# Skips PDF build (requires paper/spsurv.bbl).
#
# Usage (from package root):
#   ./paper/run-overnight.sh
#
# Monitor:
#   tail -f paper/logs/workflow-*.log
#   wc -l paper/simulation/output/results-*.txt

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mkdir -p paper/logs paper/logs/tmp paper/simulation/output

export SPSURV_ROOT="$ROOT"
export SPSURV_MC_REPLICATES=1000
export SPSURV_MC_NSIZES=50,100
export SPSURV_DEGREE_MC_REPS=1000
export SPSURV_DEGREE_MC_NSIZES=50,100
export SPSURV_MC_BAYES_CORES=1
unset SPSURV_MC_FAST_BAYES || true
unset SPSURV_MC_SEQUENTIAL || true

Rscript paper/preflight-overnight.R

LOG="paper/logs/workflow-$(date +%Y%m%d-%H%M%S).log"
echo "Starting overnight workflow; log=$LOG"
echo "Config: R=1000, n={50,100}, m=ceiling(n^0.5), parallel MC, Stan defaults, --skip-pdf"
nohup Rscript paper/run-full-workflow.R --skip-pdf >>"$LOG" 2>&1 &
echo $! > paper/logs/workflow.pid
echo "PID=$(cat paper/logs/workflow.pid)"
echo "Log=$LOG"
