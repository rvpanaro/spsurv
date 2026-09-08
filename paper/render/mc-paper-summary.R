# Shared Monte Carlo summaries for paper figures 007--008 and appendix tables.

mc_paper_nsize_levels <- function(nsize) {
  as.character(sort(unique(as.integer(as.character(nsize)))))
}

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
  nsize_levels <- mc_paper_nsize_levels(rep_df$nsize)
  rep_df$nsize_f <- factor(as.character(rep_df$nsize), levels = nsize_levels)
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
    nsize = factor(key_df[, 5], levels = nsize_levels),
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
  n_levels <- mc_paper_nsize_levels(plot_df$nsize)

  pick <- function(df, gen, par, nsize, approach) {
    row <- df[
      df$generator_cell == gen &
        df$par == par &
        as.character(df$nsize) == as.character(nsize) &
        df$approach == approach,
      ,
      drop = FALSE
    ]
    if (!nrow(row)) {
      return(NULL)
    }
    row[1L, , drop = FALSE]
  }

  fmt_triple <- function(vals, pattern) {
    paste(sprintf(pattern, vals), collapse = " & ")
  }

  fmt_rb_triple <- function(vals) {
    r <- round(as.numeric(vals), 1)
    r[abs(r) < 0.05] <- 0
    paste(sprintf("% .1f", r), collapse = " & ")
  }

  rows <- character(0)
  for (gen in gen_order) {
    sub <- plot_df[plot_df$generator_cell == gen, , drop = FALSE]
    if (!nrow(sub)) {
      next
    }
    fit <- as.character(unique(sub$model))[1L]
    if (gen %in% c("LLPH") && length(rows)) {
      rows <- c(rows, "\\addlinespace[2pt]")
    }
    for (par in c("continuous", "binary")) {
      cov_mle <- cov_bayes <- rb_mle <- rb_bayes <- se_mle <- se_bayes <- numeric(0)
      ok <- TRUE
      for (nsize in n_levels) {
        mle <- pick(sub, gen, par, nsize, "MLE")
        bayes <- pick(sub, gen, par, nsize, "Bayes")
        if (is.null(mle) || is.null(bayes)) {
          ok <- FALSE
          break
        }
        cov_mle <- c(cov_mle, mle$cov)
        cov_bayes <- c(cov_bayes, bayes$cov)
        rb_mle <- c(rb_mle, mle$rb)
        rb_bayes <- c(rb_bayes, bayes$rb)
        se_mle <- c(se_mle, mle$calib)
        se_bayes <- c(se_bayes, bayes$calib)
      }
      if (!ok) {
        next
      }
      rows <- c(
        rows,
        sprintf(
          "%s  & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          gen,
          fit,
          par_tex[[par]],
          fmt_triple(cov_mle, "%.1f"),
          fmt_triple(cov_bayes, "%.1f"),
          fmt_rb_triple(rb_mle),
          fmt_rb_triple(rb_bayes),
          fmt_triple(se_mle, "%.2f"),
          fmt_triple(se_bayes, "%.2f")
        )
      )
    }
  }
  rows
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
  n_levels <- mc_paper_nsize_levels(deg_tbl$nsize)
  param_tex <- c(
    continuous = "Cont. (age)",
    binary = "Bin. (sex)"
  )
  param_key <- function(p) {
    ifelse(p %in% c("age", "continuous"), "continuous", "binary")
  }
  rule_tex <- function(rule) {
    sub("^n\\^(.*)$", "n^{\\1}", as.character(rule))
  }
  pick <- function(df, rule, par, nsize) {
    row <- df[
      as.character(df$degree_rule) == as.character(rule) &
        param_key(df$parameter) == par &
        as.character(df$nsize) == as.character(nsize),
      ,
      drop = FALSE
    ]
    if (!nrow(row)) {
      return(NULL)
    }
    row[1L, , drop = FALSE]
  }
  fmt_n <- function(vals, pattern) {
    paste(sprintf(pattern, vals), collapse = " & ")
  }
  coverage_as_pct <- max(deg_tbl$coverage, na.rm = TRUE) <= 1

  rows <- character(0)
  for (rule in levels(deg_tbl$degree_rule)) {
    if (!any(as.character(deg_tbl$degree_rule) == as.character(rule))) {
      next
    }
    for (par in c("continuous", "binary")) {
      deg <- cov <- rb <- se <- numeric(0)
      ok <- TRUE
      for (nsize in n_levels) {
        row <- pick(deg_tbl, rule, par, nsize)
        if (is.null(row)) {
          ok <- FALSE
          break
        }
        deg <- c(deg, as.integer(row$degree))
        cov_raw <- as.numeric(row$coverage)
        cov <- c(cov, if (coverage_as_pct) 100 * cov_raw else cov_raw)
        rb <- c(rb, as.numeric(row$bias))
        se <- c(se, as.numeric(row$se_ratio))
      }
      if (!ok) {
        next
      }
      rows <- c(
        rows,
        sprintf(
          "\\(%s\\) & %s & %s & %s & %s & %s \\\\",
          rule_tex(rule),
          param_tex[[par]],
          fmt_n(deg, "%2d"),
          fmt_n(cov, "%.1f"),
          fmt_n(rb, "% .1f"),
          fmt_n(se, "%.2f")
        )
      )
    }
  }
  rows
}

#' Mean event percentages for Table tab:mc-design-events.
#'
#' Uses one approach only (MLE by default); MLE and Bayes share realisations.
summarize_mc_design_events <- function(
    censoring,
    approach = "mle",
    shape = 1.5,
    scale = 1) {
  if ("event_rate" %in% names(censoring) && !"event_proportion" %in% names(censoring)) {
    censoring$event_proportion <- censoring$event_rate
  }
  need <- c("nsize", "gdist", "model", "event_proportion")
  miss <- setdiff(need, names(censoring))
  if (length(miss)) {
    stop(
      "censoring table missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  if ("approach" %in% names(censoring)) {
    cen <- censoring[
      censoring$approach == approach &
        censoring$model %in% c("ph", "po", "aft"),
      ,
      drop = FALSE
    ]
  } else {
    cen <- censoring[censoring$model %in% c("ph", "po", "aft"), , drop = FALSE]
  }
  if (!nrow(cen)) {
    stop("No censoring rows for approach=", approach, call. = FALSE)
  }
  agg <- stats::aggregate(
    event_proportion ~ nsize + gdist + model,
    data = cen,
    FUN = mean,
    na.rm = TRUE
  )
  agg$event_pct <- 100 * agg$event_proportion
  agg$shape <- shape
  agg$scale <- scale
  agg
}

format_mc_design_events_tex_body <- function(event_df) {
  gdist_tex <- c(weibull = "Weibull", llogis = "Log-logistic")
  model_tex <- c(ph = "PH", po = "PO", aft = "AFT")
  n_levels <- as.integer(mc_paper_nsize_levels(event_df$nsize))
  if (!length(n_levels)) {
    n_levels <- c(50L, 100L, 200L)
  }
  rows <- character(0)
  for (g in c("weibull", "llogis")) {
    for (m in c("ph", "po", "aft")) {
      pcts <- numeric(0)
      shape <- scale <- NA_real_
      ok <- TRUE
      for (n in n_levels) {
        r <- event_df[
          event_df$gdist == g & event_df$model == m & event_df$nsize == n,
          ,
          drop = FALSE
        ]
        if (!nrow(r)) {
          ok <- FALSE
          break
        }
        pcts <- c(pcts, r$event_pct[[1L]])
        shape <- r$shape[[1L]]
        scale <- r$scale[[1L]]
      }
      if (!ok) {
        next
      }
      rows <- c(
        rows,
        sprintf(
          "%s      & %s  & %g & %g & %s \\\\",
          gdist_tex[[g]],
          model_tex[[m]],
          shape,
          scale,
          paste(sprintf("%.1f", pcts), collapse = " & ")
        )
      )
    }
  }
  rows
}
