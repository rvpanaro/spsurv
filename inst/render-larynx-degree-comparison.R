# Larynx BPPH: compare m = 3 MLE, default m = 10 MLE, and m = 10 Bayes.
# Illustrates delta-method instability vs lower degree vs posterior bands.
#
# Usage (from package root):
#   Rscript inst/render-larynx-degree-comparison.R
#   Rscript inst/render-larynx-degree-comparison.R /path/to/output.pdf

render_larynx_degree_comparison <- function(
    output = NULL,
    age_ref = 65,
    degree_low = 3L,
    degree_high = NULL,
    width = 9,
    height = 10.5) {
  pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkg_root <- normalizePath(".", winslash = "/")
  }

  suppressPackageStartupMessages({
    if (!"package:spsurv" %in% search() &&
      requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(pkg_root, quiet = TRUE)
    } else if (!"package:spsurv" %in% search()) {
      library(spsurv)
    }
    library(KMsurv)
    library(ggplot2)
    library(patchwork)
  })

  data("larynx", package = "KMsurv")
  larynx$stage <- factor(larynx$stage)
  if (is.null(degree_high)) {
    degree_high <- as.integer(ceiling(nrow(larynx)^(0.5)))
  }

  newdata <- data.frame(
    age = age_ref,
    stage = factor(levels(larynx$stage), levels = levels(larynx$stage))
  )
  plot_times <- seq(0, max(larynx$time), length.out = 121)
  curve_labs <- paste("Stage", levels(newdata$stage))

  .mle_panel_label <- function(fit) {
    diag <- spsurv:::.spbp_gamma_information_diagnostics(fit)
    status <- if (isTRUE(diag$stable)) {
      "stable gamma information"
    } else {
      "ill-conditioned gamma block (bands may be unreliable)"
    }
    sprintf(
      "m = %d  |  kappa(gamma) = %s\n%s",
      length(fit$bp.param),
      format(signif(diag$kappa_gamma, 3), scientific = TRUE),
      status
    )
  }

  .predict_long <- function(fit, label, interval.type = NULL) {
    args <- list(
      object = fit,
      newdata = newdata,
      times = plot_times,
      interval = 0.95
    )
    if (!is.null(interval.type)) {
      args$interval.type <- interval.type
    }
    pr <- do.call(predict, args)
    pr$series <- curve_labs[match(as.character(pr$id), as.character(seq_len(nrow(newdata))))]
    pr$panel <- label
    pr$has_bands <- is.finite(pr$lower) & is.finite(pr$upper)
    pr
  }

  .survival_panel <- function(pr, title, subtitle) {
    ggplot(pr, aes(x = time, y = surv, color = series, fill = series)) +
      geom_ribbon(
        data = subset(pr, has_bands),
        aes(ymin = lower, ymax = upper),
        alpha = 0.22,
        colour = NA
      ) +
      geom_line(linewidth = 0.55) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
      labs(
        x = "Time (years)",
        y = "Survival probability",
        title = title,
        subtitle = subtitle,
        color = NULL,
        fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  }

  fit_low <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "mle",
    degree = degree_low,
    init = 0
  )
  fit_high <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "mle",
    degree = degree_high,
    init = 0
  )
  set.seed(1)
  fit_bayes <- bpph(
    Surv(time, delta) ~ age + stage,
    data = larynx,
    approach = "bayes",
    degree = degree_high
  )

  pr_low <- .predict_long(fit_low, .mle_panel_label(fit_low))
  pr_high <- .predict_long(fit_high, .mle_panel_label(fit_high))
  pr_bayes <- .predict_long(
    fit_bayes,
    sprintf("m = %d  |  Bayes BPPH\nposterior 95%% HPD bands", degree_high),
    interval.type = "hpd"
  )

  p_low <- .survival_panel(
    pr_low,
    title = sprintf("MLE, m = %d", degree_low),
    subtitle = "Stable gamma information: pointwise 95% delta-method bands"
  )
  p_high <- .survival_panel(
    pr_high,
    title = sprintf("MLE, m = %d (default)", degree_high),
    subtitle = "Ill-conditioned gamma block: delta-method bands with survfit() warning"
  )
  p_bayes <- .survival_panel(
    pr_bayes,
    title = sprintf("Bayes, m = %d (default)", degree_high),
    subtitle = "Posterior HPD bands; no delta-method or gamma-information warning"
  )

  p <- p_low / p_high / p_bayes +
    plot_annotation(
      title = sprintf(
        "Larynx BPPH at age %d years: degree and uncertainty route",
        age_ref
      ),
      caption = paste0(
        "n = ", nrow(larynx), ", events = ", sum(larynx$delta),
        ". Same covariate profiles in all panels; compare stable MLE bands at ",
        "low degree, unreliable delta-method bands at default degree, and ",
        "Bayesian posterior bands at the same default degree."
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (is.null(output)) {
    output <- file.path(pkg_root, "figures", "006_larynx_degree_comparison.pdf")
  }
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(output, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  message("Wrote ", output)
  invisible(list(
    plot = p,
    data = rbind(pr_low, pr_high, pr_bayes),
    fits = list(low = fit_low, high = fit_high, bayes = fit_bayes)
  ))
}

if (identical(sys.nframe(), 0L)) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[[1L]] else NULL
  render_larynx_degree_comparison(output = out)
}
