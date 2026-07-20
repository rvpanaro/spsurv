# Monte Carlo regression study for spsurv.TeX Table tab:bp-mcsim-abc.
#
# Usage (from package root):
#   Rscript -e 'devtools::load_all(".", quiet=TRUE); source("paper/simulation/monte_carlo_regression.R")'
#
# Production defaults: R=1000, n in {50, 100}, parallel replicates (snowfall).
# MLE only (Bayes archived separately):
#   SPSURV_MC_MLE_ONLY=1 Rscript -e '... source("paper/simulation/monte_carlo_regression.R")'
#   Then: Rscript paper/simulation/combine-mc-bayes-mle.R
# Sequential smoke test only:
#   SPSURV_MC_REPLICATES=2 SPSURV_MC_NSIZES=100 SPSURV_MC_FAST_BAYES=1 \
#   SPSURV_MC_SEQUENTIAL=1 Rscript -e 'devtools::load_all(".", quiet=TRUE); source("paper/simulation/monte_carlo_regression.R")'
#
# Production Bayes fits use \code{spbp}/\code{rstan} defaults (4 chains,
# \code{iter = 2000}, \code{warmup = 1000}); set \code{SPSURV_MC_FAST_BAYES=1}
# only for smoke tests.

if (!exists("pkg_root", inherits = FALSE)) {
  paths <- source_paper_paths()
  pkg_root <- paths$pkg_root
}

paths <- if (exists("paths", inherits = FALSE)) paths else source_paper_paths()
sim_dir <- paths$sim_dir
output_dir <- paths$sim_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

req_pkgs <- c("survival", "rsurv", "flexsurv")
miss <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop(
    "Missing package(s) for Monte Carlo simulation: ",
    paste(miss, collapse = ", "),
    call. = FALSE
  )
}

if (!"package:spsurv" %in% search()) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    stop("Load spsurv with devtools::load_all() before sourcing this script.", call. = FALSE)
  }
}

parse_int_vec <- function(x, default) {
  if (!nzchar(x)) {
    return(default)
  }
  as.integer(strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1L]])
}

R <- suppressWarnings(as.integer(Sys.getenv("SPSURV_MC_REPLICATES", "1000")))
if (!is.finite(R) || R < 1L) {
  R <- 1000L
}

nsizes <- parse_int_vec(Sys.getenv("SPSURV_MC_NSIZES", ""), c(50L, 100L))
nsizes <- nsizes[is.finite(nsizes) & nsizes > 0L]
if (!length(nsizes)) {
  nsizes <- c(50L, 100L)
}

approaches <- {
  if (nzchar(Sys.getenv("SPSURV_MC_MLE_ONLY", ""))) {
    c("mle")
  } else if (nzchar(Sys.getenv("SPSURV_MC_BAYES_ONLY", ""))) {
    c("bayes")
  } else {
    raw <- Sys.getenv("SPSURV_MC_APPROACHES", "mle,bayes")
    out <- strsplit(gsub("\\s+", "", raw), ",", fixed = TRUE)[[1L]]
    out <- out[out %in% c("mle", "bayes")]
    if (!length(out)) {
      c("mle", "bayes")
    } else {
      out
    }
  }
}
gdists <- c("weibull", "llogis")
models <- c("ph", "po", "aft")
if (nzchar(Sys.getenv("SPSURV_MC_CORRECT_SPECIFIED_ONLY", ""))) {
  models <- c("aft")
}

beta_dgp <- c(age = -2, sexm = 1)
shape <- 1.5
scale <- 1
tau <- 10

fast_bayes <- nzchar(Sys.getenv("SPSURV_MC_FAST_BAYES", ""))
sequential <- nzchar(Sys.getenv("SPSURV_MC_SEQUENTIAL", ""))
use_parallel <- !sequential &&
  requireNamespace("snowfall", quietly = TRUE) &&
  R > 1L

bayes_fit_args <- function(approach) {
  if (!identical(approach, "bayes")) {
    return(list())
  }
  if (fast_bayes) {
    return(list(
      chains = 1L,
      iter = 150L,
      warmup = 75L,
      cores = 1L,
      verbose = FALSE
    ))
  }
  cores_env <- suppressWarnings(as.integer(Sys.getenv("SPSURV_MC_BAYES_CORES", "")))
  if (is.finite(cores_env) && cores_env >= 1L) {
    return(list(cores = cores_env))
  }
  if (use_parallel) {
    return(list(cores = 1L))
  }
  list()
}

resume_results <- Sys.getenv("SPSURV_MC_RESUME_RESULTS", "")
if (nzchar(resume_results)) {
  results_path <- resume_results
  dirn <- dirname(results_path)
  base <- basename(results_path)
  errors_path <- file.path(dirn, sub("^results-", "errors-", base))
  censoring_path <- file.path(dirn, sub("^results-", "censoring-", base))
  tag <- sub("^results-", "", sub("\\.txt$", "", base))
} else {
  tag <- format(Sys.time(), "%Y%m%d-%H%M%S")
  tag <- paste0(tag, "-", Sys.getpid())
  results_path <- file.path(output_dir, paste0("results-", tag, ".txt"))
  errors_path <- file.path(output_dir, paste0("errors-", tag, ".txt"))
  censoring_path <- file.path(output_dir, paste0("censoring-", tag, ".txt"))
}

completed_scenarios <- character()
if (file.exists(results_path) && file.info(results_path)$size > 0) {
  # Keep this light: only read columns needed to identify completed cells.
  # Columns (no header): 1=nsize, 10=gdist, 11=approach, 12=model
  cc <- rep("NULL", 13L)
  cc[c(1L, 10L, 11L, 12L)] <- c("integer", "character", "character", "character")
  existing <- utils::read.table(
    results_path,
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    colClasses = cc
  )
  if (nrow(existing)) {
    completed_scenarios <- unique(
      paste(existing$V1, existing$V10, existing$V11, existing$V12)
    )
  }
}

truth_for_names <- function(par_names) {
  out <- rep(NA_real_, length(par_names))
  names(out) <- par_names
  if ("age" %in% names(out)) {
    out["age"] <- beta_dgp[["age"]]
  }
  sx <- grep("^sex", names(out), value = TRUE)
  if (length(sx)) {
    out[sx] <- beta_dgp[["sexm"]]
  }
  out
}

se_spbp <- function(object) {
  if (identical(object$call$approach, "mle")) {
    sqrt(diag(as.matrix(stats::vcov(object))))
  } else {
    apply(object$posterior$beta, 2L, stats::sd)
  }
}

rsurv_package <- function(gdist) {
  if (gdist == "weibull") {
    "stats"
  } else {
    "flexsurv"
  }
}

simulate_dataset <- function(r, nsize, gdist, model) {
  set.seed(r)
  dat <- data.frame(
    age = stats::rnorm(nsize),
    sex = factor(
      sample(c("f", "m"), size = nsize, replace = TRUE),
      levels = c("f", "m")
    )
  )
  beta_vec <- unname(c(beta_dgp[["age"]], beta_dgp[["sexm"]]))
  pkg <- rsurv_package(gdist)
  u <- stats::runif(nsize)
  gen_args <- list(
    u = u,
    formula = ~age + sex,
    beta = beta_vec,
    dist = gdist,
    shape = shape,
    scale = scale,
    package = pkg,
    data = dat
  )
  t_ev <- switch(
    model,
    ph = do.call(rsurv::rphreg, gen_args),
    po = do.call(rsurv::rporeg, gen_args),
    aft = do.call(rsurv::raftreg, gen_args),
    stop("Unknown model: ", model, call. = FALSE)
  )
  c_ev <- stats::runif(nsize, 0, tau)
  dat$time <- pmin(t_ev, c_ev)
  dat$status <- as.integer(t_ev <= c_ev)
  list(data = dat, event_proportion = mean(dat$status))
}

fit_spbp_row <- function(dat, nsize, gdist, approach, model, rep_id) {
  m <- as.integer(ceiling(nsize^0.5))
  fit_fun <- switch(
    model,
    ph = spsurv::bpph,
    po = spsurv::bppo,
    aft = spsurv::bpaft,
    stop("Unknown model: ", model, call. = FALSE)
  )
  fit <- try(
    do.call(
      fit_fun,
      c(
        list(
          formula = survival::Surv(time, status) ~ age + sex,
          data = dat,
          degree = m,
          approach = approach
        ),
        if (identical(approach, "bayes")) bayes_fit_args(approach) else list()
      )
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(NULL)
  }

  estimates <- stats::coef(fit)
  par <- names(estimates)
  truth <- truth_for_names(par)
  ci <- try(stats::confint(fit, parm = par), silent = TRUE)
  if (inherits(ci, "try-error") || is.null(ci)) {
    return(NULL)
  }
  SE <- se_spbp(fit)
  names(SE) <- par
  rb <- 100 * (estimates - truth) / pmax(abs(truth), .Machine$double.eps)
  cp <- truth > ci[, 1L] & truth < ci[, 2L]

  data.frame(
    nsize = nsize,
    par = par,
    real = as.numeric(truth[par]),
    estimate = as.numeric(estimates),
    se = as.numeric(SE[par]),
    RB = as.numeric(rb),
    lwr = as.numeric(ci[, 1L]),
    upr = as.numeric(ci[, 2L]),
    CP = as.logical(cp),
    gdist = gdist,
    approach = approach,
    model = model,
    rep = rep_id,
    stringsAsFactors = FALSE
  )
}

fit_flexsurv_row <- function(dat, nsize, gdist, rep_id) {
  distr <- if (gdist == "weibull") "weibull" else "llogis"
  fit <- try(
    flexsurv::flexsurvreg(
      survival::Surv(time, status) ~ age + sex,
      data = dat,
      dist = distr
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  co <- stats::coef(fit)
  par <- names(co)
  truth <- truth_for_names(par)
  se <- sqrt(diag(stats::vcov(fit)))[par]
  ci <- cbind(co - 1.96 * se, co + 1.96 * se)
  rb <- 100 * (co - truth) / pmax(abs(truth), .Machine$double.eps)
  cp <- truth > ci[, 1L] & truth < ci[, 2L]
  data.frame(
    nsize = nsize,
    par = par,
    real = as.numeric(truth[par]),
    estimate = as.numeric(co),
    se = as.numeric(se),
    RB = as.numeric(rb),
    lwr = as.numeric(ci[, 1L]),
    upr = as.numeric(ci[, 2L]),
    CP = as.logical(cp),
    gdist = gdist,
    approach = "mle",
    model = "flexsurv",
    rep = rep_id,
    stringsAsFactors = FALSE
  )
}

run_monte_carlo <- function(r, nsize, gdist, approach, model) {
  sim <- simulate_dataset(r, nsize, gdist, model)
  rows <- fit_spbp_row(sim$data, nsize, gdist, approach, model, r)
  list(
    rows = rows,
    censoring = data.frame(
      nsize = nsize,
      gdist = gdist,
      approach = approach,
      model = model,
      rep = r,
      event_proportion = sim$event_proportion,
      stringsAsFactors = FALSE
    ),
    error = data.frame(
      nsize = nsize,
      gdist = gdist,
      approach = approach,
      model = model,
      code = if (is.null(rows)) 1L else 0L,
      rep = r,
      stringsAsFactors = FALSE
    )
  )
}

run_flexsurv_mc <- function(r, nsize, gdist) {
  sim <- simulate_dataset(r, nsize, gdist, "aft")
  rows <- fit_flexsurv_row(sim$data, nsize, gdist, r)
  list(
    rows = rows,
    censoring = data.frame(
      nsize = nsize,
      gdist = gdist,
      approach = "mle",
      model = "flexsurv",
      rep = r,
      event_proportion = sim$event_proportion,
      stringsAsFactors = FALSE
    ),
    error = data.frame(
      nsize = nsize,
      gdist = gdist,
      approach = "mle",
      model = "flexsurv",
      code = if (is.null(rows)) 1L else 0L,
      rep = r,
      stringsAsFactors = FALSE
    )
  )
}

append_block <- function(block, results_path, errors_path, censoring_path) {
  if (!is.null(block$rows) && nrow(block$rows)) {
    utils::write.table(
      block$rows,
      file = results_path,
      append = file.exists(results_path) && file.info(results_path)$size > 0,
      row.names = FALSE,
      col.names = FALSE,
      quote = TRUE
    )
  }
  utils::write.table(
    block$error,
    file = errors_path,
    append = file.exists(errors_path) && file.info(errors_path)$size > 0,
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
  utils::write.table(
    block$censoring,
    file = censoring_path,
    append = file.exists(censoring_path) && file.info(censoring_path)$size > 0,
    row.names = FALSE,
    col.names = FALSE,
    quote = TRUE
  )
}

run_block <- function(reps, nsize, gdist, approach, model) {
  lapply(reps, function(r) {
    if (identical(model, "flexsurv")) {
      run_flexsurv_mc(r, nsize, gdist)
    } else {
      run_monte_carlo(r, nsize, gdist, approach, model)
    }
  })
}

message(
  "Monte Carlo: R=", R,
  ", nsizes=", paste(nsizes, collapse = ","),
  if (fast_bayes) " (fast Bayes smoke)" else " (spbp/rstan default Bayes MCMC)",
  if (sequential) " (sequential)" else if (use_parallel) " (parallel)" else " (sequential)"
)
message("Writing: ", results_path)
if (length(completed_scenarios)) {
  message(
    "Resume: found ", length(completed_scenarios),
    " completed scenario cell(s) in existing results; will skip them."
  )
}

reps <- seq_len(R)

run_parallel_reps <- function(reps, worker_fun, export_vars) {
  workers_env <- suppressWarnings(as.integer(Sys.getenv("SPSURV_MC_WORKERS", "")))
  ncpus <- if (is.finite(workers_env) && workers_env >= 1L) {
    workers_env
  } else {
    max(1L, parallel::detectCores(logical = TRUE) - 2L)
  }
  snowfall::sfInit(parallel = TRUE, cpus = ncpus)
  on.exit(snowfall::sfStop(), add = TRUE)
  # sfLibrary(spsurv) loads the *installed* package and ignores devtools::load_all()
  # in the parent process. Reload package source on each worker instead.
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required for parallel MC with load_all() source.", call. = FALSE)
  }
  snowfall::sfExport("pkg_root")
  snowfall::sfClusterEval(devtools::load_all(pkg_root, quiet = TRUE, reset = FALSE))
  snowfall::sfLibrary(survival)
  snowfall::sfLibrary(rsurv)
  snowfall::sfLibrary(flexsurv)
  do.call(snowfall::sfExport, as.list(export_vars))
  snowfall::sfLapply(reps, worker_fun)
}

for (nsize in nsizes) {
  for (approach in approaches) {
    for (gdist in gdists) {
      for (model in models) {
        scenario_key <- paste(nsize, gdist, approach, model)
        if (scenario_key %in% completed_scenarios) {
          message(
            "Skipping completed scenario n=", nsize, " ", gdist, " ", approach, " ", model,
            " ..."
          )
          next
        }
        message(
          "Scenario n=", nsize, " ", gdist, " ", approach, " ", model,
          " ..."
        )
        if (use_parallel) {
          blocks <- run_parallel_reps(
            reps,
            function(r) run_monte_carlo(r, nsize, gdist, approach, model),
            c(
              "run_monte_carlo", "simulate_dataset", "fit_spbp_row",
              "truth_for_names", "se_spbp", "beta_dgp", "shape", "scale",
              "tau", "bayes_fit_args", "fast_bayes", "rsurv_package", "nsize", "gdist",
              "approach", "model"
            )
          )
        } else {
          blocks <- run_block(reps, nsize, gdist, approach, model)
        }
        for (block in blocks) {
          append_block(block, results_path, errors_path, censoring_path)
        }
      }
      if (identical(approach, "mle") && "aft" %in% models) {
        scenario_key <- paste(nsize, gdist, "mle", "flexsurv")
        if (scenario_key %in% completed_scenarios) {
          message("Skipping completed scenario n=", nsize, " ", gdist, " mle flexsurv ...")
          next
        }
        message("Scenario n=", nsize, " ", gdist, " mle flexsurv ...")
        if (use_parallel) {
          blocks <- run_parallel_reps(
            reps,
            function(r) run_flexsurv_mc(r, nsize, gdist),
            c(
              "run_flexsurv_mc", "simulate_dataset", "fit_flexsurv_row",
              "truth_for_names", "beta_dgp", "shape", "scale", "tau",
              "rsurv_package", "nsize", "gdist"
            )
          )
        } else {
          blocks <- lapply(reps, function(r) run_flexsurv_mc(r, nsize, gdist))
        }
        for (block in blocks) {
          append_block(block, results_path, errors_path, censoring_path)
        }
      }
    }
  }
}

message("Done. Results: ", results_path)
message("Censoring: ", censoring_path)
message("Errors: ", errors_path)
if (file.exists(file.path(paths$paper_dir, "paper-runtime.R"))) {
  source(file.path(paths$paper_dir, "paper-runtime.R"), local = FALSE)
  manifest_path <- paper_write_manifest(
    paths,
    list(
      run_id = tag,
      step = "mc_regression_complete",
      results_path = results_path,
      censoring_path = censoring_path,
      errors_path = errors_path,
      replicates = R,
      nsizes = nsizes,
      fast_bayes = fast_bayes,
      sequential = sequential,
      stan_defaults = if (!fast_bayes) paper_spbp_stan_defaults() else "smoke_fast_bayes"
    )
  )
  message("Manifest: ", manifest_path)
}
message("Next: Rscript paper/simulation/build_bp_mcsim_rds.R")
