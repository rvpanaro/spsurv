#!/usr/bin/env Rscript
# Degree-sensitivity study for spsurv.TeX figure 008 (LLPH -> BPPH, MLE).
#
# Output: paper/degree_llph_bpph.csv
#
# Usage (from package root):
#   Rscript paper/simulation/build_degree_llph_bpph_csv.R
#
# Smoke test:
#   SPSURV_DEGREE_MC_REPS=5 SPSURV_DEGREE_MC_NSIZES=50 Rscript paper/simulation/build_degree_llph_bpph_csv.R

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}
paths <- source_paper_paths()
pkg_root <- paths$pkg_root

suppressPackageStartupMessages({
  if (!requireNamespace("rsurv", quietly = TRUE)) {
    stop("Install rsurv before running the degree study.", call. = FALSE)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Install survival before running the degree study.", call. = FALSE)
  }
  if (!"package:spsurv" %in% search()) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(pkg_root, quiet = TRUE)
    } else {
      stop("Load spsurv with devtools::load_all() first.", call. = FALSE)
    }
  }
})

parse_int_vec <- function(x, default) {
  if (!nzchar(x)) {
    return(default)
  }
  as.integer(strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1L]])
}

R <- suppressWarnings(as.integer(Sys.getenv("SPSURV_DEGREE_MC_REPS", "1000")))
if (!is.finite(R) || R < 1L) {
  R <- 1000L
}

nsizes <- parse_int_vec(Sys.getenv("SPSURV_DEGREE_MC_NSIZES", ""), c(50L, 100L))
nsizes <- nsizes[is.finite(nsizes) & nsizes > 0L]
if (!length(nsizes)) {
  nsizes <- c(50L, 100L)
}
degree_exponents <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
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

degree_rule_label <- function(degree, nsize) {
  for (p in degree_exponents) {
    if (degree == as.integer(ceiling(nsize^p))) {
      return(paste0("n^", p))
    }
  }
  paste0("m=", degree)
}

generate_llph_dataset <- function(r, nsize) {
  set.seed(100000L + 1000L * nsize + r)
  dat <- data.frame(
    age = stats::rnorm(nsize),
    sex = factor(
      sample(c("f", "m"), size = nsize, replace = TRUE),
      levels = c("f", "m")
    )
  )
  t_ev <- rsurv::rphreg(
    stats::runif(nsize),
    ~age + sex,
    beta = unname(c(beta_dgp[["age"]], beta_dgp[["sexm"]])),
    dist = "llogis",
    shape = shape,
    scale = scale,
    package = "flexsurv",
    data = dat
  )
  c_ev <- stats::runif(nsize, 0, tau)
  dat$time <- pmin(t_ev, c_ev)
  dat$status <- as.integer(t_ev <= c_ev)
  dat
}

fit_degree_rows <- function(dat, nsize, degree, rep_id) {
  fit <- try(
    spsurv::bpph(
      survival::Surv(time, status) ~ age + sex,
      data = dat,
      degree = degree,
      approach = "mle"
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  est <- stats::coef(fit)
  par <- names(est)
  truth <- truth_for_names(par)
  ci <- try(stats::confint(fit, parm = par), silent = TRUE)
  if (inherits(ci, "try-error") || is.null(ci)) {
    return(NULL)
  }
  se <- sqrt(diag(as.matrix(stats::vcov(fit))))[par]
  rb <- 100 * (est - truth) / pmax(abs(truth), .Machine$double.eps)
  cp <- truth > ci[, 1L] & truth < ci[, 2L]
  data.frame(
    rep = rep_id,
    nsize = nsize,
    degree = degree,
    par = par,
    estimate = as.numeric(est),
    se = as.numeric(se),
    rb = as.numeric(rb),
    cp = as.numeric(cp),
    stringsAsFactors = FALSE
  )
}

all_rows <- list()
idx <- 1L
for (nsize in nsizes) {
  deg_grid <- unique(vapply(
    degree_exponents,
    function(p) as.integer(ceiling(nsize^p)),
    integer(1L)
  ))
  for (r in seq_len(R)) {
    dat <- generate_llph_dataset(r, nsize)
    for (deg in deg_grid) {
      out <- fit_degree_rows(dat, nsize, deg, r)
      if (!is.null(out)) {
        all_rows[[idx]] <- out
        idx <- idx + 1L
      }
    }
  }
}

res <- do.call(rbind, all_rows)
if (is.null(res) || !nrow(res)) {
  stop("No successful replicates in degree study.", call. = FALSE)
}

split_keys <- split(res, list(res$nsize, res$degree, res$par), drop = TRUE)
summ <- do.call(rbind, lapply(split_keys, function(dd) {
  par_raw <- unique(dd$par)
  par_label <- ifelse(par_raw == "sexm", "binary", par_raw)
  data.frame(
    nsize = unique(dd$nsize),
    degree = unique(dd$degree),
    degree_rule = degree_rule_label(unique(dd$degree), unique(dd$nsize)),
    parameter = par_label,
    coverage = mean(dd$cp, na.rm = TRUE),
    bias = mean(dd$rb, na.rm = TRUE),
    se_ratio = mean(dd$se, na.rm = TRUE) / stats::sd(dd$estimate, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

summ$parameter <- factor(summ$parameter, levels = c("age", "binary"))
summ <- summ[order(summ$nsize, summ$degree, summ$parameter), ]
summ$parameter <- as.character(summ$parameter)
summ$coverage <- round(summ$coverage, 3)
summ$bias <- round(summ$bias, 3)
summ$se_ratio <- round(summ$se_ratio, 3)

paper_dir <- paths$paper_dir
dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- paths$degree_csv
utils::write.csv(summ, out_csv, row.names = FALSE)
message("Wrote ", out_csv)
if (file.exists(file.path(paths$paper_dir, "paper-runtime.R"))) {
  source(file.path(paths$paper_dir, "paper-runtime.R"), local = FALSE)
  manifest_path <- paper_write_manifest(
    paths,
    list(
      run_id = format(Sys.time(), "%Y%m%d-%H%M%S"),
      step = "degree_llph_bpph_complete",
      output_csv = out_csv,
      replicates = R,
      nsizes = nsizes,
      nrow_summary = nrow(summ),
      approach = "mle"
    )
  )
  message("Manifest: ", manifest_path)
}
