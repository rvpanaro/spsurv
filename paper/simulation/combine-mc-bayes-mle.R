#!/usr/bin/env Rscript
# Combine archived Bayes MC results with a fresh MLE-only run.
#
# Usage (from package root):
#   Rscript paper/simulation/combine-mc-bayes-mle.R
#   Rscript paper/simulation/combine-mc-bayes-mle.R \
#     --bayes-rds=paper/archive/mc-bayes/bp-mcsim-results-bayes.rds \
#     --mle-results=paper/simulation/output/results-<tag>.txt
#
# Writes:
#   paper/bp-mcsim-results.rds
#   paper/simulation/output/results-combined-<timestamp>.txt (+ censoring, errors)

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$sim_dir, "R", "simulation_io.R"), local = TRUE)

args <- commandArgs(trailingOnly = TRUE)
bayes_rds <- paths$paper_dir
bayes_rds <- file.path(bayes_rds, "archive", "mc-bayes", "bp-mcsim-results-bayes.rds")
mle_results <- NA_character_
if (length(args)) {
  br <- grep("^--bayes-rds=", args, value = TRUE)
  if (length(br)) {
    bayes_rds <- sub("^--bayes-rds=", "", br[[1L]])
  }
  mr <- grep("^--mle-results=", args, value = TRUE)
  if (length(mr)) {
    mle_results <- sub("^--mle-results=", "", mr[[1L]])
  }
}

if (!file.exists(bayes_rds)) {
  stop("Missing Bayes archive: ", bayes_rds, call. = FALSE)
}
if (!nzchar(mle_results)) {
  mle_results <- find_latest_output("results", paths$sim_output_dir)
}
if (is.null(mle_results) || !file.exists(mle_results)) {
  stop("No MLE results-*.txt found under ", paths$sim_output_dir, call. = FALSE)
}

mle_errors <- sub("^results-", "errors-", basename(mle_results))
mle_errors <- file.path(dirname(mle_results), mle_errors)
mle_censoring <- sub("^results-", "censoring-", basename(mle_results))
mle_censoring <- file.path(dirname(mle_results), mle_censoring)

read_optional_table <- function(path, col_names) {
  if (!file.exists(path)) {
    return(NULL)
  }
  d <- utils::read.table(path, header = FALSE, quote = "\"", stringsAsFactors = FALSE)
  names(d) <- col_names
  d
}

bayes_obj <- readRDS(bayes_rds)
bayes_rep <- bayes_obj$replicates
if (!is.data.frame(bayes_rep) || !nrow(bayes_rep)) {
  stop("No Bayes replicates in ", bayes_rds, call. = FALSE)
}

mle_rep <- read_results_table(mle_results)
mle_rep <- mle_rep[mle_rep$approach == "mle", , drop = FALSE]
if (!nrow(mle_rep)) {
  stop("No MLE rows in ", mle_results, call. = FALSE)
}

bayes_err <- read_optional_table(
  file.path(dirname(bayes_rds), "errors-bayes.txt"),
  c("nsize", "gdist", "approach", "model", "code", "rep")
)
if (is.null(bayes_err) || !nrow(bayes_err)) {
  manifest_path <- file.path(paths$paper_dir, "archive", "mc-bayes", "manifest-latest.rds")
  if (file.exists(manifest_path)) {
    man <- readRDS(manifest_path)
    ep <- man$bayes$errors
    if (is.character(ep) && file.exists(ep)) {
      bayes_err <- read_optional_table(
        ep,
        c("nsize", "gdist", "approach", "model", "code", "rep")
      )
    }
  }
}

mle_err <- read_optional_table(
  mle_errors,
  c("nsize", "gdist", "approach", "model", "code", "rep")
)
bayes_cen <- read_optional_table(
  file.path(dirname(bayes_rds), "censoring-bayes.txt"),
  c("nsize", "gdist", "approach", "model", "rep", "event_proportion")
)
if (is.null(bayes_cen) || !nrow(bayes_cen)) {
  manifest_path <- file.path(paths$paper_dir, "archive", "mc-bayes", "manifest-latest.rds")
  if (file.exists(manifest_path)) {
    man <- readRDS(manifest_path)
    cp <- man$bayes$censoring
    if (is.character(cp) && file.exists(cp)) {
      bayes_cen <- read_optional_table(
        cp,
        c("nsize", "gdist", "approach", "model", "rep", "event_proportion")
      )
    }
  }
}
mle_cen <- if (file.exists(mle_censoring)) {
  read_censoring_table(mle_censoring)
} else {
  NULL
}

replicates <- rbind(bayes_rep, mle_rep)
errors <- rbind(
  bayes_err[bayes_err$approach == "bayes", , drop = FALSE],
  mle_err[mle_err$approach == "mle", , drop = FALSE]
)
censoring <- rbind(
  bayes_cen[bayes_cen$approach == "bayes", , drop = FALSE],
  mle_cen[mle_cen$approach == "mle", , drop = FALSE]
)

build_summary <- function(rep_df, cen_df) {
  event_pct <- NULL
  if (!is.null(cen_df) && nrow(cen_df)) {
    event_pct <- aggregate(
      event_proportion ~ nsize + gdist + approach + model,
      data = cen_df,
      FUN = mean,
      na.rm = TRUE
    )
    names(event_pct)[names(event_pct) == "event_proportion"] <- "event_pct"
    event_pct$event_pct <- 100 * event_pct$event_pct
  }
  split_keys <- split(
    rep_df,
    interaction(
      rep_df$nsize,
      rep_df$gdist,
      rep_df$model,
      rep_df$approach,
      rep_df$par,
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

write_table <- function(df, path) {
  utils::write.table(df, path, row.names = FALSE, col.names = FALSE, quote = TRUE)
}

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
combined_results <- file.path(
  paths$sim_output_dir,
  paste0("results-combined-", stamp, ".txt")
)
combined_errors <- file.path(
  paths$sim_output_dir,
  paste0("errors-combined-", stamp, ".txt")
)
combined_censoring <- file.path(
  paths$sim_output_dir,
  paste0("censoring-combined-", stamp, ".txt")
)
dir.create(paths$sim_output_dir, recursive = TRUE, showWarnings = FALSE)
write_table(replicates, combined_results)
write_table(errors, combined_errors)
write_table(censoring, combined_censoring)

summary <- build_summary(replicates, censoring)
meta <- list(
  bayes_rds = bayes_rds,
  mle_results = mle_results,
  combined_results = combined_results,
  combined_censoring = combined_censoring,
  combined_errors = combined_errors,
  n_bayes_rows = nrow(bayes_rep),
  n_mle_rows = nrow(mle_rep),
  n_replicates = max(replicates$rep, na.rm = TRUE),
  nsizes = sort(unique(replicates$nsize)),
  built = Sys.time()
)
out_rds <- paths$mcsim_rds
saveRDS(
  list(replicates = replicates, summary = summary, meta = meta),
  out_rds
)

message("Combined ", nrow(bayes_rep), " Bayes + ", nrow(mle_rep), " MLE replicate rows")
message("Wrote ", out_rds)
message("Wrote ", combined_results)
