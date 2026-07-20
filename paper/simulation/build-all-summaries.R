#!/usr/bin/env Rscript
# Aggregate raw Monte Carlo text output into paper/ summary files.
#
# Usage (from package root):
#   Rscript paper/simulation/build-all-summaries.R

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
sim_dir <- paths$sim_dir

status <- system2("Rscript", file.path(sim_dir, "build_bp_mcsim_rds.R"))
if (!identical(status, 0L)) {
  stop("build_bp_mcsim_rds.R failed", call. = FALSE)
}

status <- system2("Rscript", file.path(sim_dir, "build_degree_llph_bpph_csv.R"))
if (!identical(status, 0L)) {
  stop("build_degree_llph_bpph_csv.R failed", call. = FALSE)
}

message("Wrote summaries:")
message("  ", paths$mcsim_rds)
message("  ", paths$degree_csv)
