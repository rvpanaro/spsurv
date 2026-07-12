# Regenerate all paper/spsurv.TeX figures (001--008) into figures/.
#
# Quick render (figures 001--006 only unless paper/ artifacts exist):
#   Rscript inst/render-paper-figures.R
#
# Full Monte Carlo pipeline (long-running; produces paper/ artifacts for 007--008):
#   Rscript -e 'devtools::load_all(".", quiet=TRUE); source("inst/simulation/monte_carlo_regression.R")'
#   Rscript inst/simulation/build_bp_mcsim_rds.R
#   Rscript inst/simulation/build_degree_llph_bpph_csv.R
#   Rscript inst/render-paper-figures.R
#
# Figures 007--008 and appendix Tables tab:bp-mcsim-abc / tab:degree-llph-bpph
# are then regenerated from paper/bp-mcsim-results.rds and paper/degree_llph_bpph.csv.
#
# End-to-end smoke test (R=1 MC + R=2 degree study):
#   Rscript inst/render-paper-figures.R --run-smoke-sim

args <- commandArgs(trailingOnly = TRUE)
run_smoke <- "--run-smoke-sim" %in% args

pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkg_root <- normalizePath(".", winslash = "/")
}
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
  source(file.path(pkg_root, "inst", "simulation", "monte_carlo_regression.R"))
  status <- system2(
    "Rscript",
    file.path(pkg_root, "inst", "simulation", "build_bp_mcsim_rds.R")
  )
  if (!identical(status, 0L)) {
    stop("build_bp_mcsim_rds.R failed", call. = FALSE)
  }
  status <- system2(
    "Rscript",
    file.path(pkg_root, "inst", "simulation", "build_degree_llph_bpph_csv.R")
  )
  if (!identical(status, 0L)) {
    stop("build_degree_llph_bpph_csv.R failed", call. = FALSE)
  }
}

status <- system2(
  "Rscript",
  file.path(pkg_root, "inst", "render-tex-assets.R")
)
if (!identical(status, 0L)) {
  stop("render-tex-assets.R failed", call. = FALSE)
}

message("Paper figure pipeline complete.")
