#!/usr/bin/env bash
# MLE-only MC regression (Bayes kept in paper/archive/mc-bayes/).
# After completion: Rscript paper/simulation/combine-mc-bayes-mle.R
#
# Usage (from package root):
#   ./paper/run-mle-only.sh
#
# Monitor:
#   tail -f paper/logs/mle-only-*.log
#   wc -l paper/simulation/output/results-*.txt

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mkdir -p paper/logs paper/logs/tmp paper/simulation/output

export SPSURV_ROOT="$ROOT"
export SPSURV_MC_REPLICATES=1000
export SPSURV_MC_NSIZES=50,100
export SPSURV_MC_MLE_ONLY=1
export SPSURV_MC_BAYES_CORES=1
unset SPSURV_MC_FAST_BAYES || true
unset SPSURV_MC_SEQUENTIAL || true

rm -f paper/logs/workflow.pid
Rscript paper/preflight-overnight.R

LOG="paper/logs/mle-only-$(date +%Y%m%d-%H%M%S).log"
echo "Starting MLE-only MC; log=$LOG"
echo "LOG=$LOG" > paper/logs/mle-only.latest
Rscript -e 'source("paper/paths.R"); source("paper/paper-runtime.R"); paths <- source_paper_paths(); paper_configure_runtime(paths); devtools::load_all(".", quiet=TRUE); source("paper/simulation/monte_carlo_regression.R")' \
  >>"$LOG" 2>&1
echo "EXIT=$?" >>"$LOG"
