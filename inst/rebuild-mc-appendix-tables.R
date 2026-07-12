# Patch appendix Monte Carlo tables in paper/spsurv.TeX from the same
# paper/ artifacts used for figures 007--008.
#
# Usage (from package root):
#   Rscript inst/rebuild-mc-appendix-tables.R

pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(".", winslash = "/")
}

source(file.path(pkg_root, "inst", "mc-paper-summary.R"))

tex_path <- file.path(pkg_root, "paper", "spsurv.TeX")
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

mcsim_rds <- file.path(pkg_root, "paper", "bp-mcsim-results.rds")
degree_csv <- file.path(pkg_root, "paper", "degree_llph_bpph.csv")

updated <- FALSE

if (file.exists(mcsim_rds)) {
  x <- readRDS(mcsim_rds)
  plot_df <- summarize_bp_mcsim_replicates(x$replicates)
  tex <- replace_table_body(
    tex,
    "tab:bp-mcsim-abc",
    format_bp_mcsim_tex_body(plot_df)
  )
  updated <- TRUE
  message("Updated appendix table tab:bp-mcsim-abc from ", mcsim_rds)
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
  message("Patched Monte Carlo appendix tables in ", tex_path)
} else {
  message("No Monte Carlo appendix tables updated.")
}
