# Polynomial-degree sensitivity figure (008_degree_llph_bpph.pdf).
#
# Usage (from package root):
#   Rscript paper/render/render-degree-llph-bpph.R

for (src in c("paper/paths.R", "../paths.R")) {
  if (file.exists(src)) {
    source(src, local = FALSE)
    break
  }
}

render_degree_llph_bpph <- function(
    output = NULL,
    csv_path = NULL,
    width = 9.5,
    height = 5.2,
    base_size = 16) {
  paths <- source_paper_paths()

  if (is.null(csv_path)) {
    csv_path <- paths$figures$fig_008$summary_path
  }
  if (!file.exists(csv_path)) {
    stop("Missing degree sensitivity table: ", csv_path, call. = FALSE)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to render figure 008.", call. = FALSE)
  }

  deg_tbl <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  need_cols <- c(
    "degree_rule", "nsize", "parameter", "coverage", "bias", "se_ratio"
  )
  miss_cols <- setdiff(need_cols, names(deg_tbl))
  if (length(miss_cols)) {
    stop(
      "Input table missing column(s): ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  deg_tbl$degree_rule <- factor(
    deg_tbl$degree_rule,
    levels = c("n^0.2", "n^0.3", "n^0.4", "n^0.5", "n^0.6", "n^0.7", "n^0.8")
  )
  param_label <- ifelse(
    deg_tbl$parameter %in% c("sex", "sexm", "binary"),
    "binary",
    ifelse(deg_tbl$parameter == "age", "continuous", deg_tbl$parameter)
  )
  deg_tbl$line_group <- paste0("n=", deg_tbl$nsize, ", ", param_label)

  d_long <- rbind(
    data.frame(
      degree_rule = deg_tbl$degree_rule,
      line_group = deg_tbl$line_group,
      metric = "Coverage",
      value = deg_tbl$coverage
    ),
    data.frame(
      degree_rule = deg_tbl$degree_rule,
      line_group = deg_tbl$line_group,
      metric = "Bias (%)",
      value = deg_tbl$bias
    ),
    data.frame(
      degree_rule = deg_tbl$degree_rule,
      line_group = deg_tbl$line_group,
      metric = "SE ratio",
      value = deg_tbl$se_ratio
    )
  )
  d_long$metric <- factor(
    d_long$metric,
    levels = c("Coverage", "Bias (%)", "SE ratio")
  )

  se_ylim <- c(0.5, 1.2)
  se_rows <- d_long$metric == "SE ratio"
  n_se_total <- sum(se_rows)
  n_se_omitted <- sum(
    se_rows & (d_long$value < se_ylim[1] | d_long$value > se_ylim[2]),
    na.rm = TRUE
  )
  d_long$value[se_rows & (d_long$value < se_ylim[1] | d_long$value > se_ylim[2])] <- NA

  # Keep the SE-ratio panel on a readable scale when extreme values are omitted.
  d_anchor <- data.frame(
    degree_rule = factor(deg_tbl$degree_rule[1], levels = levels(deg_tbl$degree_rule)),
    line_group = deg_tbl$line_group[1],
    metric = factor("SE ratio", levels = levels(d_long$metric)),
    value = se_ylim,
    anchor = TRUE,
    stringsAsFactors = FALSE
  )
  d_long$anchor <- FALSE
  d_long <- rbind(d_long, d_anchor)

  ref_df <- data.frame(
    metric = c("Coverage", "Bias (%)", "SE ratio"),
    yint = c(0.95, 0, 1),
    stringsAsFactors = FALSE
  )

  if (is.null(output)) {
    output <- paths$figures$fig_008$path
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

  p <- ggplot2::ggplot(d_long, ggplot2::aes(
    x = degree_rule, y = value, group = line_group, colour = line_group
  )) +
    ggplot2::geom_hline(
      data = ref_df,
      ggplot2::aes(yintercept = yint),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.45
    ) +
    ggplot2::geom_line(
      data = subset(d_long, !anchor),
      na.rm = TRUE,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = subset(d_long, !anchor & !is.na(value)),
      size = 2.6
    ) +
    ggplot2::geom_point(
      data = subset(d_long, anchor),
      alpha = 0,
      size = 0
    ) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1L) +
    ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) +
    ggplot2::labs(x = "Degree rule", y = NULL, colour = "(nsize, parameter)") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = ggplot2::unit(10, "pt"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )

  ggplot2::ggsave(
    filename = output,
    plot = p,
    device = "pdf",
    width = width,
    height = height,
    units = "in"
  )
  message(
    "SE ratio panel ylim = [", se_ylim[1], ", ", se_ylim[2],
    "]; omitted ", n_se_omitted, " of ", n_se_total, " summaries"
  )
  message("Wrote ", output)
  invisible(list(
    plot = p,
    data = d_long,
    se_ylim = se_ylim,
    n_se_omitted = n_se_omitted,
    n_se_total = n_se_total
  ))
}

if (identical(sys.nframe(), 0L)) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[[1L]] else NULL
  render_degree_llph_bpph(output = out)
}
