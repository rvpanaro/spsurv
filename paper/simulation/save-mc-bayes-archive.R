#!/usr/bin/env Rscript
# Snapshot Bayes MC results (and split out MLE) from the latest production run.
#
# Usage (from package root):
#   Rscript paper/simulation/save-mc-bayes-archive.R
#
# Writes:
#   paper/archive/mc-bayes/<timestamp>/
#     results-bayes.txt, censoring-bayes.txt, errors-bayes.txt
#     bp-mcsim-results-bayes.rds
#     results-mle.txt, censoring-mle.txt, errors-mle.txt   (for MLE-only rework)
#     manifest.rds
#   paper/archive/mc-bayes/latest -> symlink to timestamp dir (if supported)

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$sim_dir, "R", "simulation_io.R"), local = TRUE)

output_dir <- paths$sim_output_dir
results_path <- find_latest_output("results", output_dir)
errors_path <- find_latest_output("errors", output_dir)
censoring_path <- find_latest_output("censoring", output_dir)
if (is.null(results_path)) {
  stop("No results-*.txt found under ", output_dir, call. = FALSE)
}

write_table <- function(df, path) {
  utils::write.table(df, path, row.names = FALSE, col.names = FALSE, quote = TRUE)
}

build_summary <- function(replicates, censoring_path) {
  event_pct <- NULL
  if (!is.null(censoring_path) && file.exists(censoring_path)) {
    cen <- read_censoring_table(censoring_path)
    event_pct <- aggregate(
      event_proportion ~ nsize + gdist + approach + model,
      data = cen,
      FUN = mean,
      na.rm = TRUE
    )
    names(event_pct)[names(event_pct) == "event_proportion"] <- "event_pct"
    event_pct$event_pct <- 100 * event_pct$event_pct
  }
  split_keys <- split(
    replicates,
    interaction(
      replicates$nsize,
      replicates$gdist,
      replicates$model,
      replicates$approach,
      replicates$par,
      drop = TRUE
    )
  )
  do.call(rbind, lapply(split_keys, function(dd) {
    out <- data.frame(
      nsize = unique(dd$nsize),
      gdist = unique(dd$gdist),
      model = unique(dd$model),
      approach = unique(dd$approach),
      par = unique(dd$par),
      real = mean(dd$real, na.rm = TRUE),
      mean = mean(dd$estimate, na.rm = TRUE),
      ase = mean(dd$se, na.rm = TRUE),
      asd = stats::sd(dd$estimate, na.rm = TRUE),
      rb = mean(dd$RB, na.rm = TRUE),
      CP = 100 * mean(dd$CP, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (!is.null(event_pct)) {
      merge(
        out,
        event_pct,
        by = c("nsize", "gdist", "approach", "model"),
        all.x = TRUE
      )
    } else {
      out$event_pct <- NA_real_
      out
    }
  }))
}

split_errors_censoring <- function(path, col_names) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  d <- utils::read.table(path, header = FALSE, quote = "\"", stringsAsFactors = FALSE)
  names(d) <- col_names
  d
}

replicates <- read_results_table(results_path)
errors <- split_errors_censoring(
  errors_path,
  c("nsize", "gdist", "approach", "model", "code", "rep")
)
censoring <- if (!is.null(censoring_path) && file.exists(censoring_path)) {
  read_censoring_table(censoring_path)
} else {
  NULL
}

bayes_rep <- replicates[replicates$approach == "bayes", , drop = FALSE]
mle_rep <- replicates[replicates$approach == "mle", , drop = FALSE]
if (!nrow(bayes_rep)) {
  stop("No Bayes rows in ", results_path, call. = FALSE)
}

bayes_err <- if (!is.null(errors)) errors[errors$approach == "bayes", , drop = FALSE] else NULL
mle_err <- if (!is.null(errors)) errors[errors$approach == "mle", , drop = FALSE] else NULL
bayes_cen <- if (!is.null(censoring)) {
  censoring[censoring$approach == "bayes", , drop = FALSE]
} else {
  NULL
}
mle_cen <- if (!is.null(censoring)) {
  censoring[censoring$approach == "mle", , drop = FALSE]
} else {
  NULL
}

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
archive_root <- file.path(paths$paper_dir, "archive", "mc-bayes")
archive_dir <- file.path(archive_root, stamp)
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

bayes_results_path <- file.path(archive_dir, "results-bayes.txt")
bayes_errors_path <- file.path(archive_dir, "errors-bayes.txt")
bayes_censoring_path <- file.path(archive_dir, "censoring-bayes.txt")
bayes_rds_path <- file.path(archive_dir, "bp-mcsim-results-bayes.rds")

mle_results_path <- file.path(archive_dir, "results-mle.txt")
mle_errors_path <- file.path(archive_dir, "errors-mle.txt")
mle_censoring_path <- file.path(archive_dir, "censoring-mle.txt")

write_table(bayes_rep, bayes_results_path)
write_table(mle_rep, mle_results_path)
if (!is.null(bayes_err) && nrow(bayes_err)) {
  write_table(bayes_err, bayes_errors_path)
}
if (!is.null(mle_err) && nrow(mle_err)) {
  write_table(mle_err, mle_errors_path)
}
if (!is.null(bayes_cen) && nrow(bayes_cen)) {
  write_table(bayes_cen, bayes_censoring_path)
}
if (!is.null(mle_cen) && nrow(mle_cen)) {
  write_table(mle_cen, mle_censoring_path)
}

bayes_summary <- build_summary(bayes_rep, bayes_censoring_path)
bayes_meta <- list(
  approach = "bayes",
  source_results = results_path,
  source_censoring = censoring_path,
  source_errors = errors_path,
  n_replicate_rows = nrow(bayes_rep),
  n_replicates = max(bayes_rep$rep, na.rm = TRUE),
  nsizes = sort(unique(bayes_rep$nsize)),
  models = sort(unique(bayes_rep$model)),
  gdists = sort(unique(bayes_rep$gdist)),
  archived = Sys.time()
)
saveRDS(
  list(replicates = bayes_rep, summary = bayes_summary, meta = bayes_meta),
  bayes_rds_path
)

# Stable copies at archive root for easy loading
stable_bayes_rds <- file.path(archive_root, "bp-mcsim-results-bayes.rds")
file.copy(bayes_rds_path, stable_bayes_rds, overwrite = TRUE)
file.copy(bayes_results_path, file.path(archive_root, "results-bayes.txt"), overwrite = TRUE)
if (file.exists(bayes_censoring_path)) {
  file.copy(bayes_censoring_path, file.path(archive_root, "censoring-bayes.txt"), overwrite = TRUE)
}
if (file.exists(bayes_errors_path)) {
  file.copy(bayes_errors_path, file.path(archive_root, "errors-bayes.txt"), overwrite = TRUE)
}

manifest <- list(
  archive_dir = archive_dir,
  stamp = stamp,
  source = list(
    results = results_path,
    censoring = censoring_path,
    errors = errors_path
  ),
  bayes = list(
    results = bayes_results_path,
    censoring = if (file.exists(bayes_censoring_path)) bayes_censoring_path else NA_character_,
    errors = if (file.exists(bayes_errors_path)) bayes_errors_path else NA_character_,
    rds = bayes_rds_path,
    stable_rds = stable_bayes_rds,
    n_rows = nrow(bayes_rep),
    n_summary_rows = nrow(bayes_summary)
  ),
  mle_split = list(
    results = mle_results_path,
    censoring = if (file.exists(mle_censoring_path)) mle_censoring_path else NA_character_,
    errors = if (file.exists(mle_errors_path)) mle_errors_path else NA_character_,
    n_rows = nrow(mle_rep),
    note = "Use for MLE-only rework; Bayes preserved separately."
  ),
  full_combined_preserved = list(
    results = results_path,
    note = "Original combined run left in simulation/output/."
  )
)
saveRDS(manifest, file.path(archive_dir, "manifest.rds"))
saveRDS(manifest, file.path(archive_root, "manifest-latest.rds"))

message("Bayes archive written to: ", archive_dir)
message("  results: ", bayes_results_path, " (", nrow(bayes_rep), " rows)")
message("  rds:     ", bayes_rds_path)
message("  stable:  ", stable_bayes_rds)
message("MLE split: ", mle_results_path, " (", nrow(mle_rep), " rows)")
