#!/usr/bin/env Rscript
# Import Fabio/anexos Monte Carlo CSVs into paper/bp-mcsim-results.rds
# (source of Table tab:bp-mcsim-abc) and a censoring table for tab:mc-design-events.
#
# Usage (from package root):
#   Rscript paper/simulation/import_anexos_mc.R
#   SPSURV_ANEXOS_DIR=/path/to/resultados Rscript paper/simulation/import_anexos_mc.R
#
# Main-text tables use n in {50, 100, 200} (n = 500 in the anexos files is ignored).

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$render_dir, "mc-paper-summary.R"))

anexos_dir <- Sys.getenv(
  "SPSURV_ANEXOS_DIR",
  unset = "/Users/rvpanaro/Downloads/anexos/resultados"
)
rep_path <- file.path(anexos_dir, "main_replications.csv")
ev_path <- file.path(anexos_dir, "main_event_rates_raw.csv")
err_path <- file.path(anexos_dir, "main_errors.csv")
if (!file.exists(rep_path) || !file.exists(ev_path)) {
  stop("Anexos CSVs not found under ", anexos_dir, call. = FALSE)
}

paper_nsizes <- c(50L, 100L, 200L)

rep <- utils::read.csv(rep_path, stringsAsFactors = FALSE)
ev <- utils::read.csv(ev_path, stringsAsFactors = FALSE)
rep <- rep[rep$nsize %in% paper_nsizes, , drop = FALSE]
ev <- ev[ev$nsize %in% paper_nsizes, , drop = FALSE]

need_rep <- c(
  "nsize", "gdist", "model", "approach", "rep", "par",
  "real", "estimate", "se", "rb", "lower", "upper", "coverage"
)
miss <- setdiff(need_rep, names(rep))
if (length(miss)) {
  stop("main_replications.csv missing: ", paste(miss, collapse = ", "), call. = FALSE)
}

replicates <- data.frame(
  nsize = as.integer(rep$nsize),
  par = as.character(rep$par),
  real = as.numeric(rep$real),
  estimate = as.numeric(rep$estimate),
  se = as.numeric(rep$se),
  RB = as.numeric(rep$rb),
  lwr = as.numeric(rep$lower),
  upr = as.numeric(rep$upper),
  CP = as.logical(rep$coverage),
  gdist = as.character(rep$gdist),
  approach = as.character(rep$approach),
  model = as.character(rep$model),
  rep = as.integer(rep$rep),
  degree = as.integer(rep$degree),
  stringsAsFactors = FALSE
)

ok <- is.finite(replicates$estimate) &
  is.finite(replicates$se) &
  is.finite(replicates$RB) &
  !is.na(replicates$CP)
replicates <- replicates[ok, , drop = FALSE]

censoring <- data.frame(
  nsize = as.integer(ev$nsize),
  gdist = as.character(ev$gdist),
  approach = "mle",
  model = as.character(ev$model),
  rep = as.integer(ev$rep),
  event_proportion = as.numeric(ev$event_rate),
  stringsAsFactors = FALSE
)

event_pct <- stats::aggregate(
  event_proportion ~ nsize + gdist + approach + model,
  data = censoring,
  FUN = mean,
  na.rm = TRUE
)
names(event_pct)[names(event_pct) == "event_proportion"] <- "event_pct"
event_pct$event_pct <- 100 * event_pct$event_pct

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
  merge(out, event_pct, by = c("nsize", "gdist", "approach", "model"), all.x = TRUE)
}))

dir.create(paths$sim_output_dir, recursive = TRUE, showWarnings = FALSE)
censoring_path <- file.path(paths$sim_output_dir, "censoring-anexos-main.txt")
utils::write.table(
  censoring[, c("nsize", "gdist", "approach", "model", "rep", "event_proportion")],
  censoring_path,
  row.names = FALSE,
  col.names = FALSE,
  quote = TRUE
)

meta <- list(
  results_path = normalizePath(rep_path, winslash = "/"),
  censoring_path = normalizePath(censoring_path, winslash = "/"),
  errors_path = if (file.exists(err_path)) {
    normalizePath(err_path, winslash = "/")
  } else {
    NA_character_
  },
  n_replicates = 1000L,
  nsizes = paper_nsizes,
  degree_rule = "m = ceiling(n^0.4)",
  source = "anexos/resultados",
  built = Sys.time()
)

saveRDS(
  list(replicates = replicates, summary = summary, meta = meta),
  paths$mcsim_rds
)

message(
  "Wrote ", paths$mcsim_rds, " (", nrow(replicates), " replicate rows; n=",
  paste(paper_nsizes, collapse = ","), ")"
)
message("Wrote ", censoring_path)
