#!/usr/bin/env Rscript
# Preflight checks before a production overnight paper run.
#
# Usage (from package root):
#   Rscript paper/preflight-overnight.R

for (src in c("paper/paths.R", "paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
source(file.path("paper", "paper-runtime.R"))

paths <- source_paper_paths()
paper_configure_runtime(paths)

req_pkgs <- c(
  "devtools", "survival", "rsurv", "flexsurv", "rstan",
  "snowfall", "ggplot2", "patchwork", "dplyr"
)
miss <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop(
    "Install missing package(s) before overnight run: ",
    paste(miss, collapse = ", "),
    call. = FALSE
  )
}

pid_path <- file.path(paths$paper_dir, "logs", "workflow.pid")
if (file.exists(pid_path)) {
  pid_txt <- trimws(readLines(pid_path, warn = FALSE)[1L])
  pid_num <- suppressWarnings(as.integer(pid_txt))
  if (is.finite(pid_num)) {
    alive <- identical(
      suppressWarnings(system2("ps", c("-p", pid_num), stdout = FALSE, stderr = FALSE)),
      0L
    )
    if (alive) {
      stop(
        "Workflow already running (PID ", pid_num, "). Stop it before a new overnight run.",
        call. = FALSE
      )
    }
  }
}

workers <- max(1L, parallel::detectCores(logical = TRUE) - 2L)
workers_env <- suppressWarnings(as.integer(Sys.getenv("SPSURV_MC_WORKERS", "")))
if (is.finite(workers_env) && workers_env >= 1L) {
  workers <- workers_env
}

cat("=== Overnight preflight OK ===\n")
cat("Package root:", paths$pkg_root, "\n")
cat("R:", R.version.string, "\n")
cat("Cores (logical):", parallel::detectCores(logical = TRUE), "\n")
cat("MC workers (snowfall):", workers, "\n")
cat("Bayes cores per fit:", Sys.getenv("SPSURV_MC_BAYES_CORES", "1"), "\n")
cat("Replicates:", Sys.getenv("SPSURV_MC_REPLICATES", "1000"), "\n")
cat("Sample sizes:", Sys.getenv("SPSURV_MC_NSIZES", "50,100"), "\n")
cat("Degree rule: m = ceiling(n^0.5)\n")
cat("Fast Bayes smoke:", nzchar(Sys.getenv("SPSURV_MC_FAST_BAYES", "")), "\n")
cat("Estimated runtime: ~7-9 h (simulation + figures; PDF skipped without spsurv.bbl)\n")
cat("Monitor: tail -f paper/logs/workflow-*.log\n")
cat("MC progress: wc -l paper/simulation/output/results-*.txt\n")
