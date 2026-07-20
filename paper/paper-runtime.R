# Paper-local runtime configuration (temp dirs, manifests, isolation checks).

paper_runtime_dirs <- function(paper_dir) {
  list(
    logs = file.path(paper_dir, "logs"),
    tmp = file.path(paper_dir, "logs", "tmp"),
    manifests = file.path(paper_dir, "logs", "manifests"),
    build = file.path(paper_dir, "build")
  )
}

paper_configure_runtime <- function(paths) {
  dirs <- paper_runtime_dirs(paths$paper_dir)
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  tmp <- dirs$tmp
  Sys.setenv(
    TMPDIR = tmp,
    TEMP = tmp,
    TMP = tmp,
    SPSURV_ROOT = paths$pkg_root
  )
  dirs
}

paper_git_revision <- function(pkg_root) {
  out <- suppressWarnings(
    system2(
      "git",
      c("-C", pkg_root, "rev-parse", "HEAD"),
      stdout = TRUE,
      stderr = FALSE
    )
  )
  if (length(out)) out[[1L]] else NA_character_
}

paper_spbp_stan_defaults <- function() {
  list(
    interface = "spbp.default() -> rstan::sampling()",
    chains = 4L,
    iter = 2000L,
    warmup = 1000L,
    thin = 1L,
    adapt_delta = 0.8,
    max_treedepth = 10L,
    cores = "spbp.default() -> .spbp_default_cores()",
    note = paste(
      "Production paper simulation passes no Bayes overrides to bpph/bppo/bpaft;",
      "only SPSURV_MC_FAST_BAYES=1 uses reduced MCMC for smoke tests."
    )
  )
}

paper_configure_simulation_env <- function() {
  Sys.unsetenv("SPSURV_MC_FAST_BAYES")
  Sys.unsetenv("SPSURV_MC_SEQUENTIAL")
  Sys.setenv(
    SPSURV_MC_REPLICATES = "1000",
    SPSURV_MC_NSIZES = "50,100",
    SPSURV_DEGREE_MC_REPS = "1000",
    SPSURV_DEGREE_MC_NSIZES = "50,100",
    SPSURV_MC_BAYES_CORES = "1"
  )
  invisible(NULL)
}

paper_simulation_cases <- function() {
  list(
    mc_regression = list(
      id = "mc_regression",
      script = "paper/simulation/monte_carlo_regression.R",
      replicates = 1000L,
      nsizes = c(50L, 100L),
      parallel = "snowfall (detectCores() - 2); opt-in sequential via SPSURV_MC_SEQUENTIAL=1",
      approaches = c("mle", "bayes"),
      gdists = c("weibull", "llogis"),
      models = c("ph", "po", "aft"),
      degree_rule = "m = ceiling(n^0.5)",
      seed_schedule = "set.seed(r) per replicate in simulate_dataset()",
      outputs = c(
        "paper/simulation/output/results-*.txt",
        "paper/simulation/output/censoring-*.txt",
        "paper/simulation/output/errors-*.txt",
        "paper/bp-mcsim-results.rds"
      )
    ),
    degree_llph_bpph = list(
      id = "degree_llph_bpph",
      script = "paper/simulation/build_degree_llph_bpph_csv.R",
      replicates = 1000L,
      nsizes = c(50L, 100L),
      approach = "mle",
      model = "bpph",
      generator = "llogis PH",
      degree_exponents = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
      seed_schedule = "set.seed(100000 + 1000*nsize + r) in generate_llph_dataset()",
      outputs = c("paper/degree_llph_bpph.csv")
    )
  )
}

paper_write_manifest <- function(paths, manifest) {
  dirs <- paper_runtime_dirs(paths$paper_dir)
  dir.create(dirs$manifests, recursive = TRUE, showWarnings = FALSE)
  run_id <- manifest$run_id %||% format(Sys.time(), "%Y%m%d-%H%M%S")
  manifest$run_id <- run_id
  manifest$updated <- as.character(Sys.time())
  path <- file.path(dirs$manifests, paste0("run-", run_id, ".rds"))
  saveRDS(manifest, path)
  message("Wrote manifest: ", path)
  invisible(path)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

paper_log_path <- function(paths, step) {
  dirs <- paper_runtime_dirs(paths$paper_dir)
  file.path(dirs$logs, paste0(step, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".log"))
}

paper_run_logged <- function(paths, step, expr) {
  log_path <- paper_log_path(paths, step)
  message("Logging ", step, " to ", log_path)
  con <- file(log_path, open = "wt")
  sink(con, split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)
  force(expr)
}

paper_check_isolation <- function(pkg_root, paper_dir) {
  pkg_root <- normalizePath(pkg_root, winslash = "/")
  paper_dir <- normalizePath(paper_dir, winslash = "/")
  offenders <- character()
  check_glob <- function(glob, label) {
    hits <- Sys.glob(file.path(pkg_root, glob))
    hits <- hits[!startsWith(normalizePath(hits, winslash = "/"), paper_dir)]
    if (length(hits)) {
      offenders <<- c(offenders, paste0(label, ": ", hits))
    }
  }
  check_glob("Rplots.pdf", "plot")
  check_glob("_spsurv_bayes_*.rds", "bayes cache")
  check_glob("results-*.txt", "mc raw")
  check_glob("censoring-*.txt", "mc raw")
  check_glob("degree_llph_bpph.csv", "degree csv")
  list(ok = !length(offenders), offenders = offenders)
}
