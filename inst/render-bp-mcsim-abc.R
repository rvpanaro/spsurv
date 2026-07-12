# Monte Carlo summary figure (007_bp_mcsim_abc.pdf): coverage, bias, SE calibration.
#
# Usage (from package root):
#   Rscript inst/render-bp-mcsim-abc.R
#   Rscript inst/render-bp-mcsim-abc.R /path/to/output.pdf

render_bp_mcsim_abc <- function(
    output = NULL,
    results_rds = NULL,
    width = 10.5,
    height = 9) {
  pkg_root <- Sys.getenv("SPSURV_ROOT", unset = normalizePath("..", winslash = "/"))
  if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkg_root <- normalizePath(".", winslash = "/")
  }

  if (is.null(results_rds)) {
    results_rds <- file.path(pkg_root, "paper", "bp-mcsim-results.rds")
  }
  if (!file.exists(results_rds)) {
    stop("Missing simulation results: ", results_rds, call. = FALSE)
  }

  req <- c("ggplot2", "patchwork", "dplyr")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing package(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }

  source(file.path(pkg_root, "inst", "mc-paper-summary.R"))

  x <- readRDS(results_rds)
  rep_df <- x$replicates
  if (!is.data.frame(rep_df) || !nrow(rep_df)) {
    stop("No replicates data in: ", results_rds, call. = FALSE)
  }

  plot_df <- summarize_bp_mcsim_replicates(rep_df)

  offset <- 0.14
  plot_df$x <- as.numeric(plot_df$model) +
    ifelse(plot_df$approach == "MLE", -offset, offset)

  make_segments <- function(df, ycol) {
    byv <- c("model", "par", "generator_cell", "nsize")
    a <- df[df$approach == "MLE", c(byv, "x"), drop = FALSE]
    b <- df[df$approach == "Bayes", c(byv, "x"), drop = FALSE]
    names(a)[names(a) == "x"] <- "x_mle"
    names(b)[names(b) == "x"] <- "x_bayes"
    m <- merge(a, b, by = byv, all = FALSE)
    a_y <- df[df$approach == "MLE", c(byv, ycol), drop = FALSE]
    b_y <- df[df$approach == "Bayes", c(byv, ycol), drop = FALSE]
    names(a_y)[names(a_y) == ycol] <- "y_mle"
    names(b_y)[names(b_y) == ycol] <- "y_bayes"
    m <- merge(m, a_y, by = byv, all = FALSE)
    merge(m, b_y, by = byv, all = FALSE)
  }

  segA <- make_segments(plot_df, "rb")
  segB <- make_segments(plot_df, "calib")
  segC <- make_segments(plot_df, "cov")

  shape_map <- c(continuous = 16, binary = 17)
  panel_aspect <- 0.5
  nsize_lab <- function(x) paste0("n = ", x)

  base_theme <- ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = ggplot2::unit(10, "pt"),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio = panel_aspect
    )

  facet_ns <- ggplot2::facet_wrap(
    ~nsize,
    nrow = 1L,
    scales = "fixed",
    labeller = ggplot2::labeller(nsize = nsize_lab)
  )

  shared_guides <- ggplot2::guides(
    shape = ggplot2::guide_legend(order = 1),
    colour = ggplot2::guide_legend(order = 2),
    alpha = ggplot2::guide_legend(order = 3)
  )

  pA <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segA,
      ggplot2::aes(
        x = x_mle, xend = x_bayes, y = y_mle, yend = y_bayes,
        colour = generator_cell
      ),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = x, y = rb, colour = generator_cell, shape = par, alpha = approach
      ),
      size = 3.5
    ) +
    ggplot2::scale_shape_manual(values = shape_map, name = "Parameter") +
    ggplot2::scale_alpha_manual(values = c(MLE = 0.55, Bayes = 1), name = "Approach") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(plot_df$model)),
      labels = levels(plot_df$model)
    ) +
    ggplot2::labs(x = NULL, y = "Relative bias (%)", colour = "Generator cell") +
    facet_ns +
    shared_guides +
    base_theme

  pB <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segB,
      ggplot2::aes(
        x = x_mle, xend = x_bayes, y = y_mle, yend = y_bayes,
        colour = generator_cell
      ),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 1, colour = "grey40") +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = x, y = calib, colour = generator_cell, shape = par, alpha = approach
      ),
      size = 3.5
    ) +
    ggplot2::scale_shape_manual(values = shape_map, name = "Parameter") +
    ggplot2::scale_alpha_manual(values = c(MLE = 0.55, Bayes = 1), name = "Approach") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(plot_df$model)),
      labels = levels(plot_df$model)
    ) +
    ggplot2::labs(x = NULL, y = "SE calibration ratio", colour = "Generator cell") +
    facet_ns +
    shared_guides +
    base_theme

  pC <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = segC,
      ggplot2::aes(
        x = x_mle, xend = x_bayes, y = y_mle, yend = y_bayes,
        colour = generator_cell
      ),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 95, colour = "grey40") +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        x = x, y = cov, colour = generator_cell, shape = par, alpha = approach
      ),
      size = 3.5
    ) +
    ggplot2::scale_shape_manual(values = shape_map, name = "Parameter") +
    ggplot2::scale_alpha_manual(values = c(MLE = 0.55, Bayes = 1), name = "Approach") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(plot_df$model)),
      labels = levels(plot_df$model)
    ) +
    ggplot2::labs(x = NULL, y = "Coverage (%)", colour = "Generator cell") +
    facet_ns +
    shared_guides +
    base_theme

  combined <- pC / pA / pB +
    patchwork::plot_layout(guides = "collect", heights = c(1, 1, 1)) &
    ggplot2::theme(legend.position = "bottom")

  if (is.null(output)) {
    output <- file.path(pkg_root, "figures", "007_bp_mcsim_abc.pdf")
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = output,
    plot = combined,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  message("Wrote ", output)
  invisible(list(plot = combined, data = plot_df))
}

if (identical(sys.nframe(), 0L)) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[[1L]] else NULL
  render_bp_mcsim_abc(output = out)
}
