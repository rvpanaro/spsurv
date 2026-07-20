#!/usr/bin/env Rscript
# Full paper workflow: simulation -> summaries -> figures -> PDF.
#
# Usage (from package root):
#   Rscript paper/run-full-workflow.R
#   Rscript paper/run-full-workflow.R --skip-simulation
#   Rscript paper/run-full-workflow.R --simulation-only
#   Rscript paper/run-full-workflow.R --skip-pdf
#   Rscript paper/run-full-workflow.R --from-step=aggregate_summaries
#
# Production Bayes MC uses spbp/rstan defaults (no SPSURV_MC_FAST_BAYES).
# Temp files: paper/logs/tmp/

args <- commandArgs(trailingOnly = TRUE)
skip_simulation <- "--skip-simulation" %in% args
simulation_only <- "--simulation-only" %in% args
skip_pdf <- "--skip-pdf" %in% args
from_step <- NA_character_
if (length(args)) {
  fs <- grep("^--from-step=", args, value = TRUE)
  if (length(fs)) {
    from_step <- sub("^--from-step=", "", fs[[1L]])
  }
}

step_order <- c(
  "degree_study",
  "mc_regression",
  "aggregate_summaries",
  "render_figures",
  "build_pdf"
)
if (skip_simulation) {
  step_order <- setdiff(step_order, c("degree_study", "mc_regression", "aggregate_summaries"))
}
if (simulation_only) {
  step_order <- intersect(step_order, c("degree_study", "mc_regression", "aggregate_summaries"))
}
if (skip_pdf) {
  step_order <- setdiff(step_order, "build_pdf")
}
if (!is.na(from_step)) {
  idx <- match(from_step, step_order)
  if (is.na(idx)) {
    stop("Unknown --from-step: ", from_step, call. = FALSE)
  }
  step_order <- step_order[seq.int(idx, length(step_order))]
}

for (src in c("paper/paths.R", "paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
source(file.path("paper", "paper-runtime.R"))

paths <- source_paper_paths()
dirs <- paper_configure_runtime(paths)
pkg_root <- paths$pkg_root

run_id <- format(Sys.time(), "%Y%m%d-%H%M%S")
manifest <- list(
  run_id = run_id,
  started = as.character(Sys.time()),
  git_revision = paper_git_revision(pkg_root),
  r_version = R.version.string,
  platform = R.version$platform,
  stan_defaults = paper_spbp_stan_defaults(),
  simulation_cases = paper_simulation_cases(),
  steps = list(),
  status = "running"
)
paper_write_manifest(paths, manifest)

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required to load spsurv for simulation.", call. = FALSE)
}

run_step <- function(name, expr) {
  if (!name %in% step_order) {
    return(invisible(FALSE))
  }
  message("\n=== ", name, " ===\n")
  log_path <- paper_log_path(paths, name)
  message("Log: ", log_path)
  t0 <- Sys.time()
  res <- tryCatch(
    {
      con <- file(log_path, open = "wt")
      sink(con, split = TRUE)
      sink(con, type = "message")
      on.exit({
        sink(type = "message")
        sink()
        close(con)
      }, add = TRUE)
      force(expr)
      "ok"
    },
    error = function(e) {
      message("ERROR in ", name, ": ", conditionMessage(e))
      "error"
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  manifest$steps[[name]] <<- list(
    status = res,
    elapsed_secs = elapsed,
    finished = as.character(Sys.time()),
    log = log_path
  )
  paper_write_manifest(paths, manifest)
  if (identical(res, "error")) {
    stop("Step failed: ", name, call. = FALSE)
  }
  invisible(TRUE)
}

if ("degree_study" %in% step_order || "mc_regression" %in% step_order) {
  old_out <- Sys.glob(file.path(paths$sim_output_dir, "*.txt"))
  if (length(old_out) && "degree_study" %in% step_order) {
    archive_dir <- file.path(paths$sim_output_dir, "archive", run_id)
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    for (f in old_out) {
      file.rename(f, file.path(archive_dir, basename(f)))
    }
    message("Archived ", length(old_out), " prior simulation output file(s)")
  }
  stale_summaries <- c(paths$mcsim_rds, paths$degree_csv)
  stale_summaries <- stale_summaries[file.exists(stale_summaries)]
  if (length(stale_summaries) && "degree_study" %in% step_order) {
    archive_dir <- file.path(paths$paper_dir, "logs", "archive", run_id)
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    for (f in stale_summaries) {
      file.rename(f, file.path(archive_dir, basename(f)))
    }
    message("Archived ", length(stale_summaries), " prior summary file(s)")
  }
  paper_configure_simulation_env()
}

run_step("degree_study", {
  devtools::load_all(pkg_root, quiet = TRUE)
  source(file.path(paths$sim_dir, "build_degree_llph_bpph_csv.R"), local = FALSE)
})

run_step("mc_regression", {
  devtools::load_all(pkg_root, quiet = TRUE)
  source(file.path(paths$sim_dir, "monte_carlo_regression.R"), local = FALSE)
})

run_step("aggregate_summaries", {
  source(file.path(paths$sim_dir, "build_bp_mcsim_rds.R"), local = FALSE)
})

run_step("render_figures", {
  status <- system2(
    "Rscript",
    c(
      file.path(paths$render_dir, "render-paper-figures.R"),
      "--figures=all",
      "--require"
    )
  )
  if (!identical(status, 0L)) {
    stop("render-paper-figures.R failed", call. = FALSE)
  }
})

run_step("build_pdf", {
  tex_name <- "spsurv.TeX"
  build_dir <- dirs$build
  work_tex <- file.path(build_dir, tex_name)
  file.copy(paths$tex_path, work_tex, overwrite = TRUE)
  figs <- Sys.glob(file.path(paths$paper_dir, "[0-9][0-9][0-9]_*.pdf"))
  for (f in figs) {
    file.copy(f, file.path(build_dir, basename(f)), overwrite = TRUE)
  }
  bbl <- file.path(paths$paper_dir, "spsurv.bbl")
  if (file.exists(bbl)) {
    file.copy(bbl, file.path(build_dir, "spsurv.bbl"), overwrite = TRUE)
  }
  owd <- setwd(build_dir)
  on.exit(setwd(owd), add = TRUE)
  for (i in seq_len(2L)) {
    status <- system2("pdflatex", c("-interaction=nonstopmode", tex_name))
    if (!identical(status, 0L)) {
      warning("pdflatex pass ", i, " returned status ", status)
    }
  }
  pdf_out <- file.path(build_dir, "spsurv.pdf")
  if (!file.exists(pdf_out)) {
    stop("pdflatex did not produce spsurv.pdf (missing spsurv.bbl?)", call. = FALSE)
  }
  final_pdf <- file.path(paths$paper_dir, "spsurv.pdf")
  file.copy(pdf_out, final_pdf, overwrite = TRUE)
  message("Wrote ", final_pdf)
})

iso <- paper_check_isolation(pkg_root, paths$paper_dir)
manifest$isolation <- iso
manifest$status <- if (iso$ok) "completed" else "completed_with_isolation_warnings"
manifest$finished <- as.character(Sys.time())
paper_write_manifest(paths, manifest)

if (!iso$ok) {
  warning(
    "Artifacts found outside paper/: ",
    paste(iso$offenders, collapse = "; ")
  )
}

message("\nPaper workflow finished (run_id=", run_id, ").")
