#' Print a coefficient table after dropping all-NA columns.
#'
#' \code{printCoefmat()} treats the last column as a p-value when
#' \code{has.Pvalue} is TRUE. That errors if the p-value column was removed
#' because every standard error was \code{NA}.
#'
#' @noRd
.spbp_print_coefmat <- function(x, digits, signif.stars, has_pvalue = FALSE, ...) {
  x <- .spbp_drop_all_na_cols(x)
  if (is.null(x) || !is.matrix(x) || ncol(x) == 0L || nrow(x) == 0L) {
    return(invisible(NULL))
  }
  cn <- colnames(x)
  nc <- ncol(x)
  has_pvalue <- isTRUE(has_pvalue) &&
    nc >= 2L &&
    length(cn) > 0L &&
    substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v")
  if (!has_pvalue && nc == 1L) {
    stats::printCoefmat(
      x,
      digits = digits,
      signif.stars = signif.stars,
      P.values = FALSE,
      has.Pvalue = FALSE,
      cs.ind = 1L,
      tst.ind = integer(0),
      ...
    )
  } else {
    stats::printCoefmat(
      x,
      digits = digits,
      signif.stars = signif.stars,
      P.values = has_pvalue,
      has.Pvalue = has_pvalue,
      ...
    )
  }
}

#' @keywords internal
#' @noRd
.print.summary.spbp <- function(x,
                                digits = 2,
                                signif.stars = getOption("show.signif.stars"),
                                approach = c("mle", "bayes"),
                                ...) {
  approach <- match.arg(approach)
  controls <- x$controls
  if (is.null(controls)) {
    controls <- if (approach == "mle") {
      list(
        compact = TRUE,
        show_call = TRUE,
        show_intervals = TRUE,
        show_baseline = FALSE,
        mle_test = "lr"
      )
    } else {
      list(
        compact = TRUE,
        show_call = TRUE,
        show_intervals = TRUE,
        show_baseline = FALSE,
        bayes_criterion = "waic"
      )
    }
  }

  if (!identical(controls$show_call, FALSE) && !is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  cat(.spbp_summary_message(x))

  if (is.null(x$coefficients) || nrow(x$coefficients) == 0) {
    cat("Null model\n")
    return(invisible())
  }

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  coef_tab <- x$coefficients
  coef_int_cols <- NULL
  if (isTRUE(controls$show_intervals) && !is.null(x$coef_interval)) {
    coef_int_tab <- x$coef_interval
    lower_col <- grep("^lower", colnames(coef_int_tab), value = TRUE)[1]
    upper_col <- grep("^upper", colnames(coef_int_tab), value = TRUE)[1]
    if (!is.na(lower_col) && !is.na(upper_col)) {
      coef_int_cols <- coef_int_tab[, c(lower_col, upper_col), drop = FALSE]
      colnames(coef_int_cols) <- .spbp_interval_colnames(colnames(coef_int_tab))
    }
  }

  if (approach == "mle") {
    reg_tab <- cbind(
      Estimate = coef_tab[, "coef"],
      coef_int_cols,
      `Std. Error` = coef_tab[, "se(coef)"],
      `z value` = coef_tab[, "z"],
      `Pr(>|z|)` = coef_tab[, "Pr(>|z|)"]
    )
    has_pvalue <- TRUE
    estimate_col <- "exp(coef)"
  } else {
    reg_tab <- cbind(
      Estimate = coef_tab[, "mean(coef)"],
      coef_int_cols,
      `Std. Error` = coef_tab[, "sd(coef)"]
    )
    has_pvalue <- FALSE
    estimate_col <- "mean(exp(coef))"
  }

  cat("Regression coefficients:\n")
  .spbp_print_coefmat(
    reg_tab,
    digits = digits,
    signif.stars = signif.stars,
    has_pvalue = has_pvalue,
    ...
  )
  cat("\n")

  if (isTRUE(controls$show_intervals) && !is.null(x$interval)) {
    int_tab <- x$interval
    lower_col <- grep("^lower", colnames(int_tab), value = TRUE)[1]
    upper_col <- grep("^upper", colnames(int_tab), value = TRUE)[1]
    if (!is.na(lower_col) && !is.na(upper_col)) {
      exp_tab <- cbind(
        Estimate = int_tab[, estimate_col],
        int_tab[, c(lower_col, upper_col), drop = FALSE]
      )
      colnames(exp_tab)[2:3] <- .spbp_interval_colnames(colnames(int_tab))
      cat("Exponentiated coefficients:\n")
      .spbp_print_coefmat(
        exp_tab,
        digits = digits,
        signif.stars = FALSE,
        has_pvalue = FALSE
      )
      cat("\n")
    }
  }

  if (isTRUE(controls$show_baseline) && !is.null(x$baseline)) {
    base_tab <- cbind(
      Estimate = x$baseline[, "gamma"],
      `Std. Error` = x$baseline[, "se_log_gamma"],
      `z value` = x$baseline[, "z"],
      `Pr(>|z|)` = x$baseline[, "Pr(>|z|)"]
    )
    cat("Bernstein baseline coefficients:\n")
    .spbp_print_coefmat(
      base_tab,
      digits = digits,
      signif.stars = signif.stars,
      has_pvalue = TRUE,
      ...
    )
    cat("\n")
  }

  cat("--- \n")
  if (approach == "mle") {
    loglik <- if (!is.null(x$loglik)) unname(x$loglik[2]) else NA_real_
    nparams <- if (!is.null(x$nparams)) x$nparams else nrow(x$coefficients)
    aic <- if (!is.null(x$loglik)) -2 * unname(x$loglik[2]) + 2 * nparams else NA_real_
    cat("loglik =", loglik, " ", "AIC =", aic, "\n")

    if (!isTRUE(controls$compact)) {
      pdig <- max(1, getOption("digits") - 4)
      show_test <- controls$mle_test
      if (identical(show_test, "wald")) {
        cat(
          "Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
          x$waldtest["df"], " df,", "   p=",
          format.pval(x$waldtest["pvalue"], digits = pdig),
          "\n",
          sep = ""
        )
      } else {
        cat(
          "Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
          x$logtest["df"], " df,", "   p=",
          format.pval(x$logtest["pvalue"], digits = pdig),
          "\n",
          sep = ""
        )
      }
      if (!identical(show_test, "wald")) {
        cat(
          "Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
          x$waldtest["df"], " df,", "   p=",
          format.pval(x$waldtest["pvalue"], digits = pdig),
          "\n",
          sep = ""
        )
      }
    }
  } else if (!isTRUE(controls$compact)) {
    cat("DIC =", x$dic, "  WAIC =", x$waic, "  LPML =", x$lpml, "\n")
  } else {
    crit <- controls$bayes_criterion
    if (identical(crit, "lpml")) {
      cat("DIC =", x$dic, "  LPML =", x$lpml, "\n")
    } else if (identical(crit, "dic")) {
      cat("DIC =", x$dic, "\n")
    } else {
      cat("DIC =", x$dic, "  WAIC =", x$waic, "\n")
    }
  }

  invisible()
}
