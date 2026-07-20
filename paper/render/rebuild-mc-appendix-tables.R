# Patch Monte Carlo tables in paper/spsurv.TeX from paper/ simulation summaries.
#
# Usage (from package root):
#   Rscript paper/render/rebuild-mc-appendix-tables.R

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$render_dir, "mc-paper-summary.R"))

tex_path <- paths$tex_path
tex <- readLines(tex_path, warn = FALSE)

replace_table_body <- function(lines, label, new_body) {
  cap_idx <- grep(sprintf("\\\\label\\{%s\\}", label), lines, fixed = FALSE)[1L]
  if (is.na(cap_idx)) {
    stop("Could not find table label: ", label, call. = FALSE)
  }
  begin_idx <- cap_idx
  while (begin_idx > 1L && !grepl("\\\\begin\\{table", lines[begin_idx])) {
    begin_idx <- begin_idx - 1L
  }
  mid_idx <- begin_idx +
    grep("\\\\midrule", lines[begin_idx:length(lines)], fixed = FALSE)[1L] - 1L
  bot_idx <- begin_idx +
    grep("\\\\bottomrule", lines[begin_idx:length(lines)], fixed = FALSE)[1L] - 1L
  if (is.na(mid_idx) || is.na(bot_idx) || bot_idx <= mid_idx) {
    stop("Could not locate table body for: ", label, call. = FALSE)
  }
  c(
    lines[seq_len(mid_idx)],
    new_body,
    lines[seq.int(bot_idx, length(lines))]
  )
}

mcsim_rds <- paths$mcsim_rds
degree_csv <- paths$degree_csv
source(file.path(paths$sim_dir, "R", "simulation_io.R"), local = TRUE)
censoring_path <- if (file.exists(mcsim_rds)) {
  meta_cen <- readRDS(mcsim_rds)$meta$censoring_path
  if (!is.null(meta_cen) && file.exists(meta_cen)) {
    meta_cen
  } else {
    find_latest_output("censoring", paths$sim_output_dir)
  }
} else {
  find_latest_output("censoring", paths$sim_output_dir)
}

updated <- FALSE

if (!is.null(censoring_path) && file.exists(censoring_path)) {
  cen <- read_censoring_table(censoring_path)
  event_df <- summarize_mc_design_events(cen)
  tex <- replace_table_body(
    tex,
    "tab:mc-design-events",
    format_mc_design_events_tex_body(event_df)
  )
  updated <- TRUE
  message("Updated table tab:mc-design-events from ", censoring_path)
} else {
  message("Skipping tab:mc-design-events (missing censoring-*.txt)")
}

if (file.exists(mcsim_rds)) {
  x <- readRDS(mcsim_rds)
  plot_df <- summarize_bp_mcsim_replicates(x$replicates)
  tex <- replace_table_body(
    tex,
    "tab:bp-mcsim-abc",
    format_bp_mcsim_tex_body(plot_df)
  )
  updated <- TRUE
  message("Updated table tab:bp-mcsim-abc from ", mcsim_rds)
} else {
  message("Skipping tab:bp-mcsim-abc (missing ", mcsim_rds, ")")
}

if (file.exists(degree_csv)) {
  deg_tbl <- utils::read.csv(degree_csv, stringsAsFactors = FALSE)
  tex <- replace_table_body(
    tex,
    "tab:degree-llph-bpph",
    format_degree_llph_tex_body(deg_tbl)
  )
  updated <- TRUE
  message("Updated appendix table tab:degree-llph-bpph from ", degree_csv)
} else {
  message("Skipping tab:degree-llph-bpph (missing ", degree_csv, ")")
}

if (updated) {
  writeLines(tex, tex_path)
  message("Patched Monte Carlo tables in ", tex_path)
} else {
  message("No Monte Carlo appendix tables updated.")
}
