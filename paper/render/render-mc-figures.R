# Render Monte Carlo figure 008 and patch simulation tables in spsurv.TeX.
#
# Usage (from package root):
#   Rscript paper/render/render-mc-figures.R
#   Rscript paper/render/render-mc-figures.R --require

args <- commandArgs(trailingOnly = TRUE)
require_summaries <- "--require" %in% args

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
render_dir <- paths$render_dir
figures <- paths$figures

render_one <- function(spec, render_fun, ...) {
  if (!file.exists(spec$summary_path)) {
    if (require_summaries) {
      stop(
        "Missing simulation summary for ", spec$pdf, ": ",
        spec$summary_path,
        call. = FALSE
      )
    }
    message("Skipping ", spec$pdf, " (missing ", basename(spec$summary_path), ")")
    return(invisible(FALSE))
  }
  render_fun(output = spec$path, ...)
  invisible(TRUE)
}

fig008 <- figures$fig_008
mcsim_rds <- paths$mcsim_rds

source(file.path(render_dir, "render-degree-llph-bpph.R"))
render_one(
  fig008,
  render_degree_llph_bpph,
  csv_path = fig008$summary_path
)

if (file.exists(mcsim_rds) || file.exists(fig008$summary_path)) {
  status <- system2(
    "Rscript",
    file.path(render_dir, "rebuild-mc-appendix-tables.R")
  )
  if (!identical(status, 0L)) {
    stop("rebuild-mc-appendix-tables.R failed", call. = FALSE)
  }
}

source(file.path(render_dir, "sync-tex-figures.R"))
illustration_specs <- paper_figures_by_kind(paths$figures, "illustration")
simulation_specs <- paper_figures_by_kind(paths$figures, "simulation")
paths_sync <- paths
paths_sync$figures <- c(illustration_specs, simulation_specs)
sync_tex_figure_paths(paths_sync, patch = TRUE, verify_files = require_summaries)

message("Monte Carlo figure pipeline complete.")
