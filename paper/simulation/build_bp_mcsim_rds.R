#!/usr/bin/env Rscript
# Aggregate paper/simulation/output/results-*.txt into paper/bp-mcsim-results.rds.
#
# Usage (from package root):
#   Rscript paper/simulation/build_bp_mcsim_rds.R

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()

sim_dir <- paths$sim_dir
source(file.path(sim_dir, "R", "simulation_io.R"), local = TRUE)

output_dir <- paths$sim_output_dir
results_path <- find_latest_output("results", output_dir)
censoring_path <- find_latest_output("censoring", output_dir)
if (is.null(results_path)) {
  stop("No results-*.txt found under ", output_dir, call. = FALSE)
}

replicates <- read_results_table(results_path)
if (!nrow(replicates)) {
  stop("Empty results table: ", results_path, call. = FALSE)
}

event_pct <- NULL
if (!is.null(censoring_path)) {
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

summary <- do.call(rbind, lapply(split_keys, function(dd) {
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

paper_dir <- paths$paper_dir
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- paths$mcsim_rds

meta <- list(
  results_path = results_path,
  censoring_path = censoring_path,
  n_replicates = max(replicates$rep, na.rm = TRUE),
  nsizes = sort(unique(replicates$nsize)),
  built = Sys.time()
)

saveRDS(
  list(
    replicates = replicates,
    summary = summary,
    meta = meta
  ),
  out_rds
)

message("Wrote ", out_rds, " (", nrow(replicates), " replicate rows)")
