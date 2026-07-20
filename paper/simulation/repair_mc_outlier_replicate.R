#!/usr/bin/env Rscript
# Remove the single worst SE outlier replicate and run one replacement.

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
source(file.path(paths$sim_dir, "R", "simulation_io.R"), local = TRUE)
suppressPackageStartupMessages(devtools::load_all(paths$pkg_root, quiet = TRUE))

output_dir <- paths$sim_output_dir
results_path <- find_latest_output("results", output_dir)
errors_path <- find_latest_output("errors", output_dir)
censoring_path <- find_latest_output("censoring", output_dir)
if (is.null(results_path) || is.null(errors_path) || is.null(censoring_path)) {
  stop("Missing Monte Carlo output files.", call. = FALSE)
}

results <- read_results_table(results_path)
errors <- utils::read.table(errors_path, header = FALSE, quote = "\"", stringsAsFactors = FALSE)
names(errors) <- c("nsize", "gdist", "approach", "model", "code", "rep")
censoring <- read_censoring_table(censoring_path)

mx <- aggregate(se ~ rep + nsize + gdist + approach + model, results, max, na.rm = TRUE)
worst <- mx[which.max(mx$se), ]
message("Removing outlier replicate:")
print(worst)

drop_results <- with(
  results,
  rep == worst$rep &
    nsize == worst$nsize &
    gdist == worst$gdist &
    approach == worst$approach &
    model == worst$model
)
message("Rows removed: ", sum(drop_results))
results <- results[!drop_results, , drop = FALSE]

drop_errors <- with(
  errors,
  rep == worst$rep &
    nsize == worst$nsize &
    gdist == worst$gdist &
    approach == worst$approach &
    model == worst$model
)
errors <- errors[!drop_errors, , drop = FALSE]

drop_censoring <- with(
  censoring,
  rep == worst$rep &
    nsize == worst$nsize &
    gdist == worst$gdist &
    approach == worst$approach &
    model == worst$model
)
censoring <- censoring[!drop_censoring, , drop = FALSE]

beta_dgp <- c(age = -2, sexm = 1)
shape <- 1.5
scale <- 1
tau <- 10

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
      list(
        formula = survival::Surv(time, status) ~ age + sex,
        data = dat,
        degree = m,
        approach = approach
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

new_rep <- max(results$rep, na.rm = TRUE) + 1L
message("Running replacement replicate: ", new_rep)
block <- run_monte_carlo(
  new_rep,
  worst$nsize,
  as.character(worst$gdist),
  as.character(worst$approach),
  as.character(worst$model)
)
if (is.null(block$rows) || !nrow(block$rows)) {
  stop("Replacement replicate failed to fit.", call. = FALSE)
}
message("Replacement max se: ", max(block$rows$se, na.rm = TRUE))

results <- rbind(results, block$rows)
errors <- rbind(errors, block$error)
censoring <- rbind(censoring, block$censoring)

write_table <- function(df, path) {
  utils::write.table(df, path, row.names = FALSE, col.names = FALSE, quote = TRUE)
}
write_table(results, results_path)
write_table(errors, errors_path)
write_table(censoring, censoring_path)
message("Updated: ", results_path)

status <- system2("Rscript", file.path(paths$sim_dir, "build_bp_mcsim_rds.R"))
if (!identical(status, 0L)) {
  stop("build_bp_mcsim_rds.R failed", call. = FALSE)
}
status <- system2(
  "Rscript",
  c(file.path(paths$render_dir, "render-paper-figures.R"), "--figures=mc", "--require")
)
if (!identical(status, 0L)) {
  stop("render-paper-figures.R failed", call. = FALSE)
}
message("Repair complete.")
