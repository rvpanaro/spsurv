# Regenerate paper/spsurv.TeX figures into paper/.
#
# Illustrations (001--006):
#   Rscript paper/render/render-paper-figures.R
#   Rscript paper/render/render-paper-figures.R --figures=illustrations
#
# Monte Carlo (007--008) from paper/ simulation summaries:
#   Rscript paper/render/render-paper-figures.R --figures=mc
#   Rscript paper/render/render-paper-figures.R --figures=mc --require
#
# Full pipeline:
#   Rscript -e 'devtools::load_all(".", quiet=TRUE); source("paper/simulation/monte_carlo_regression.R")'
#   Rscript paper/simulation/build-all-summaries.R
#   Rscript paper/render/render-paper-figures.R --figures=all --require
#
# Smoke test:
#   Rscript paper/render/render-paper-figures.R --run-smoke-sim --figures=all

args <- commandArgs(trailingOnly = TRUE)
run_smoke <- "--run-smoke-sim" %in% args
require_mc <- "--require" %in% args

figures_arg <- "all"
if (length(args)) {
  fig_flags <- grep("^--figures=", args, value = TRUE)
  if (length(fig_flags)) {
    figures_arg <- sub("^--figures=", "", fig_flags[[1L]])
  }
}
figures_arg <- match.arg(
  figures_arg,
  c("all", "illustrations", "mc")
)

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
pkg_root <- paths$pkg_root
Sys.setenv(SPSURV_ROOT = pkg_root)

if (run_smoke) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required for --run-smoke-sim", call. = FALSE)
  }
  message("Running smoke Monte Carlo (R=1, n=100, fast Bayes)...")
  devtools::load_all(pkg_root, quiet = TRUE)
  Sys.setenv(
    SPSURV_MC_REPLICATES = "1",
    SPSURV_MC_NSIZES = "100",
    SPSURV_MC_FAST_BAYES = "1",
    SPSURV_MC_SEQUENTIAL = "1",
    SPSURV_DEGREE_MC_REPS = "2"
  )
  source(file.path(paths$sim_dir, "monte_carlo_regression.R"))
  status <- system2(
    "Rscript",
    file.path(paths$sim_dir, "build-all-summaries.R")
  )
  if (!identical(status, 0L)) {
    stop("build-all-summaries.R failed", call. = FALSE)
  }
}

run_script <- function(script, extra_args = character()) {
  status <- system2("Rscript", c(script, extra_args))
  if (!identical(status, 0L)) {
    stop(basename(script), " failed", call. = FALSE)
  }
}

if (figures_arg %in% c("all", "illustrations")) {
  run_script(file.path(paths$render_dir, "render-tex-assets.R"))
}

if (figures_arg %in% c("all", "mc")) {
  mc_args <- if (require_mc) "--require" else character()
  run_script(file.path(paths$render_dir, "render-mc-figures.R"), mc_args)
}

message("Paper figure pipeline complete (", figures_arg, ").")
