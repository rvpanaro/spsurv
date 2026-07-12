# Shared Monte Carlo summaries for paper figures 007--008 and appendix tables.

summarize_bp_mcsim_replicates <- function(rep_df) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr for Monte Carlo summaries.", call. = FALSE)
  }

  need_cols <- c(
    "par", "estimate", "se", "real", "RB", "CP",
    "gdist", "model", "rep", "nsize", "approach"
  )
  miss_cols <- setdiff(need_cols, names(rep_df))
  if (length(miss_cols)) {
    stop(
      "replicates missing column(s): ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  rep_df <- rep_df[rep_df$model %in% c("ph", "po", "aft"), , drop = FALSE]
  rep_df <- rep_df[rep_df$approach %in% c("mle", "bayes"), , drop = FALSE]

  model_lab <- c(ph = "BPPH", po = "BPPO", aft = "BPAFT")
  gdist_lab <- c(weibull = "WAFT", loglogistic = "LLAFT", llogis = "LLAFT")
  approach_lab <- c(mle = "MLE", bayes = "Bayes")
  par_lab <- function(p) {
    ifelse(p == "sexm", "binary", ifelse(p == "age", "continuous", p))
  }
  gen_cell_label <- function(family_label, model_label) {
    fam_prefix <- ifelse(family_label == "WAFT", "W", "LL")
    cls <- ifelse(
      model_label == "BPPH",
      "PH",
      ifelse(model_label == "BPPO", "PO", "AFT")
    )
    paste0(fam_prefix, cls)
  }

  rep_df$model_f <- factor(
    model_lab[as.character(rep_df$model)],
    levels = c("BPPH", "BPPO", "BPAFT")
  )
  gd <- as.character(rep_df$gdist)
  gd[gd %in% names(gdist_lab)] <- gdist_lab[gd[gd %in% names(gdist_lab)]]
  rep_df$gdist_f <- factor(gd, levels = c("WAFT", "LLAFT"))
  rep_df$par_f <- factor(
    par_lab(as.character(rep_df$par)),
    levels = c("continuous", "binary")
  )
  rep_df$approach_f <- factor(
    approach_lab[as.character(rep_df$approach)],
    levels = c("MLE", "Bayes")
  )
  rep_df$nsize_f <- factor(as.character(rep_df$nsize), levels = c("100", "300"))
  rep_df$gen_cell_f <- factor(
    dplyr::recode(
      gen_cell_label(as.character(rep_df$gdist_f), as.character(rep_df$model_f)),
      "WPO" = "WPO",
      "LLPH" = "LLPH"
    ),
    levels = c("WPH", "WPO", "WAFT", "LLPH", "LLPO", "LLAFT")
  )
  rep_df$CP_num <- suppressWarnings(as.numeric(rep_df$CP))

  split_key <- interaction(
    rep_df$model_f,
    rep_df$par_f,
    rep_df$gen_cell_f,
    rep_df$approach_f,
    rep_df$nsize_f,
    drop = TRUE
  )
  rb_m <- tapply(rep_df$RB, split_key, mean, na.rm = TRUE)
  sde_m <- tapply(rep_df$estimate, split_key, stats::sd, na.rm = TRUE)
  se_mean_m <- tapply(rep_df$se, split_key, mean, na.rm = TRUE)
  calib_m <- se_mean_m / sde_m
  cov_m <- tapply(rep_df$CP_num, split_key, mean, na.rm = TRUE) * 100

  key_df <- do.call(rbind, strsplit(names(rb_m), split = "\\.", fixed = FALSE))
  if (ncol(key_df) != 5L) {
    stop(
      "Unexpected interaction key format (expected model.par.generator.approach.nsize): ncol=",
      ncol(key_df),
      call. = FALSE
    )
  }
  colnames(key_df) <- c("model", "par", "generator_cell", "approach", "nsize")
  data.frame(
    model = factor(key_df[, 1], levels = c("BPPH", "BPPO", "BPAFT")),
    par = factor(key_df[, 2], levels = c("continuous", "binary")),
    generator_cell = factor(
      key_df[, 3],
      levels = c("WPH", "WPO", "WAFT", "LLPH", "LLPO", "LLAFT")
    ),
    approach = factor(key_df[, 4], levels = c("MLE", "Bayes")),
    nsize = factor(key_df[, 5], levels = c("100", "300")),
    rb = as.numeric(rb_m),
    sde = as.numeric(sde_m),
    calib = as.numeric(calib_m),
    cov = as.numeric(cov_m),
    stringsAsFactors = FALSE
  )
}

format_bp_mcsim_tex_body <- function(plot_df) {
  par_tex <- c(
    continuous = "Cont. (age)",
    binary = "Bin. (sex)"
  )
  gen_order <- c("WPH", "WPO", "WAFT", "LLPH", "LLPO", "LLAFT")

  format_row <- function(df, nsize_label) {
    rows <- character(0)
    for (gen in gen_order) {
      sub <- df[df$generator_cell == gen & df$nsize == nsize_label, , drop = FALSE]
      if (!nrow(sub)) {
        next
      }
      fit <- as.character(unique(sub$model))[1L]
      for (par in c("continuous", "binary")) {
        par_rows <- sub[sub$par == par, , drop = FALSE]
        if (!nrow(par_rows)) {
          next
        }
        mle <- par_rows[par_rows$approach == "MLE", , drop = FALSE]
        bayes <- par_rows[par_rows$approach == "Bayes", , drop = FALSE]
        if (!nrow(mle) || !nrow(bayes)) {
          next
        }
        rows <- c(
          rows,
          sprintf(
            paste0(
              "%s   & %s  & %s & %.1f & %.1f & % .1f & % .1f & %.2f & %.2f \\\\"
            ),
            gen,
            fit,
            par_tex[[par]],
            mle$cov,
            bayes$cov,
            mle$rb,
            bayes$rb,
            mle$calib,
            bayes$calib
          )
        )
      }
    }
    c(
      sprintf("\\multicolumn{9}{@{}l}{\\textit{\\(n=%s\\)}} \\\\", nsize_label),
      "\\cmidrule(lr){1-9}",
      rows
    )
  }

  c(
    format_row(plot_df, "100"),
    "\\addlinespace[4pt]",
    format_row(plot_df, "300")
  )
}

format_degree_llph_tex_body <- function(deg_tbl) {
  need_cols <- c(
    "degree_rule", "degree", "nsize", "parameter", "coverage", "bias", "se_ratio"
  )
  miss_cols <- setdiff(need_cols, names(deg_tbl))
  if (length(miss_cols)) {
    stop(
      "degree table missing column(s): ",
      paste(miss_cols, collapse = ", "),
      call. = FALSE
    )
  }

  deg_tbl$degree_rule <- factor(
    deg_tbl$degree_rule,
    levels = c("n^0.2", "n^0.3", "n^0.4", "n^0.5", "n^0.6", "n^0.7", "n^0.8")
  )
  param_tex <- function(p) {
    ifelse(p %in% c("age", "continuous"), "Continuous (age)", "Binary (sex)")
  }

  format_block <- function(nsize_label) {
    sub <- deg_tbl[deg_tbl$nsize == as.integer(nsize_label), , drop = FALSE]
    sub <- sub[order(sub$degree_rule, sub$parameter), , drop = FALSE]
    rows <- apply(sub, 1L, function(row) {
      cov_pct <- if (max(sub$coverage, na.rm = TRUE) <= 1) {
        100 * as.numeric(row[["coverage"]])
      } else {
        as.numeric(row[["coverage"]])
      }
      sprintf(
        paste0(
          "\\(%s\\) & %2d & %s & %.1f & % .1f & %.2f \\\\"
        ),
        row[["degree_rule"]],
        as.integer(row[["degree"]]),
        param_tex(row[["parameter"]]),
        cov_pct,
        as.numeric(row[["bias"]]),
        as.numeric(row[["se_ratio"]])
      )
    })
    c(
      sprintf("\\multicolumn{6}{@{}l}{\\textit{\\(n=%s\\)}} \\\\", nsize_label),
      "\\cmidrule(lr){1-6}",
      rows
    )
  }

  c(
    format_block("100"),
    "\\addlinespace[4pt]",
    format_block("300")
  )
}
