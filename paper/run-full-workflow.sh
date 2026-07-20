#!/usr/bin/env bash
# Detached full paper workflow (simulation + figures + PDF).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p paper/logs paper/logs/tmp
export SPSURV_ROOT="$ROOT"
export SPSURV_MC_REPLICATES=1000
export SPSURV_MC_NSIZES=50,100
export SPSURV_DEGREE_MC_REPS=1000
export SPSURV_DEGREE_MC_NSIZES=50,100
unset SPSURV_MC_FAST_BAYES || true
unset SPSURV_MC_SEQUENTIAL || true
LOG="paper/logs/workflow-$(date +%Y%m%d-%H%M%S).log"
echo "Starting paper workflow; log=$LOG"
nohup Rscript paper/run-full-workflow.R >>"$LOG" 2>&1 &
echo $! > paper/logs/workflow.pid
echo "PID=$(cat paper/logs/workflow.pid)"
