#' @keywords internal
#' @noRd
.spbp_is_bernstein_spec <- function(x) {
  is.list(x) &&
    identical(x$baseline, "bernstein") &&
    !is.null(x$m)
}

#' @keywords internal
#' @noRd
.spbp_check_parametric_baseline <- function(x, arg_name = "dist") {
  if (is.character(x) && length(x) == 1L && nzchar(x)) {
    if (!identical(x, "bernstein")) {
      stop(
        "Parametric baseline '", x, "' is not supported in spsurv. ",
        "Use ", arg_name, " = bernstein(m) or baseline = bernstein(m).",
        call. = FALSE
      )
    }
  }
}

#' Resolve Bernstein polynomial degree from degree / dist / baseline
#' @keywords internal
#' @noRd
.spbp_resolve_degree <- function(degree, dist = NULL, baseline = NULL) {
  if (!missing(degree) && !is.null(degree)) {
    if (!(degree %% 1 == 0)) {
      stop("Polynomial degree must be integer.")
    }
    return(as.integer(degree))
  }

  spec <- dist
  if (is.null(spec)) {
    spec <- baseline
  } else if (!is.null(baseline)) {
    warning(
      "Both 'dist' and 'baseline' were supplied; using 'dist'.",
      call. = FALSE
    )
  }

  if (is.null(spec)) {
    return(NULL)
  }

  .spbp_check_parametric_baseline(spec, "dist")

  if (.spbp_is_bernstein_spec(spec)) {
    if (is.null(spec$m)) {
      stop("bernstein(m) requires a non-NULL degree m.", call. = FALSE)
    }
    return(as.integer(spec$m))
  }

  if (is.character(spec) && identical(spec, "bernstein")) {
    stop(
      "Use bernstein(m) to specify the Bernstein polynomial degree.",
      call. = FALSE
    )
  }

  stop(
    "Unrecognized baseline/dist specification. ",
    "Use bernstein(m) or set 'degree' directly.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.spbp_check_mle <- function(object) {
  if (!inherits(object, "spbp")) {
    stop("Object must inherit from class 'spbp'.", call. = FALSE)
  }
  if (!identical(object$call$approach, "mle")) {
    stop(
      "This method is only defined for MLE fits (approach = 'mle').",
      call. = FALSE
    )
  }
  invisible(object)
}

#' @keywords internal
#' @noRd
.spbp_model_loglik <- function(object) {
  .spbp_check_mle(object)
  if (is.null(object$loglik) || length(object$loglik) < 2L) {
    return(NA_real_)
  }
  unname(object$loglik[2])
}

#' @keywords internal
#' @noRd
.spbp_nparams <- function(object) {
  n_gamma <- length(object$bp.param)
  if (is.null(object$coefficients)) {
    return(as.integer(n_gamma))
  }
  as.integer(length(object$coefficients) + n_gamma)
}

#' @keywords internal
#' @noRd
.spbp_model_label <- function(object) {
  model <- object$call$model
  if (is.null(model)) {
    model <- "spbp"
  }
  fml <- object$formula
  if (is.null(fml) && !is.null(object$call$formula)) {
    fml <- object$call$formula
  }
  fml_chr <- if (!is.null(fml)) {
    paste(deparse(fml), collapse = " ")
  } else {
    "formula"
  }
  paste0("bernstein(", model, "): ", fml_chr)
}

#' @keywords internal
#' @noRd
.spbp_aic <- function(object) {
  ll <- .spbp_model_loglik(object)
  k <- .spbp_nparams(object)
  -2 * ll + 2 * k
}

#' Residual degrees of freedom (n - total parameters)
#' @keywords internal
#' @noRd
.spbp_resid_df <- function(object) {
  as.integer(object$n - .spbp_nparams(object))
}

#' Short formula label for anova tables
#' @keywords internal
#' @noRd
.spbp_anova_terms_label <- function(object) {
  fml <- object$formula
  if (is.null(fml) && !is.null(object$call$formula)) {
    fml <- object$call$formula
  }
  if (is.null(fml)) {
    return("model")
  }
  lbl <- attr(stats::terms(fml), "term.labels")
  if (length(lbl) == 0L) {
    return("1")
  }
  paste(lbl, collapse = " + ")
}

#' Label for the Test column in pairwise anova tables (survreg-style)
#' @keywords internal
#' @noRd
.spbp_anova_test_label <- function(prev_fit, next_fit, df) {
  if (!is.finite(df) || df <= 1L) {
    return("")
  }
  prev_lbl <- attr(stats::terms(prev_fit$formula), "term.labels")
  next_lbl <- attr(stats::terms(next_fit$formula), "term.labels")
  added <- setdiff(next_lbl, prev_lbl)
  if (length(added) != 1L) {
    return("")
  }
  paste0("+", added)
}

#' Pairwise anova table for two or more nested fits (survreg column layout)
#' @keywords internal
#' @noRd
.spbp_anova_pairwise <- function(fits) {
  n <- length(fits)
  logliks <- vapply(fits, function(f) as.numeric(logLik(f)), numeric(1))
  npars <- vapply(fits, .spbp_nparams, integer(1))
  resid_df <- vapply(fits, .spbp_resid_df, integer(1))
  terms <- vapply(fits, .spbp_anova_terms_label, character(1))

  deviance <- rep(NA_real_, n)
  df <- rep(NA_integer_, n)
  test <- rep("", n)
  pval <- rep(NA_real_, n)

  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      j <- i + 1L
      deviance[j] <- -2 * (logliks[i] - logliks[j])
      df[j] <- npars[j] - npars[i]
      test[j] <- .spbp_anova_test_label(fits[[i]], fits[[j]], df[j])
      if (is.finite(deviance[j]) && df[j] > 0L) {
        pval[j] <- stats::pchisq(deviance[j], df[j], lower.tail = FALSE)
      }
    }
  }

  data.frame(
    Terms = terms,
    `Resid. Df` = resid_df,
    `-2*LL` = -2 * logliks,
    Test = test,
    Df = df,
    Deviance = deviance,
    `Pr(>Chi)` = pval,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Heading vector for print.anova
#' @keywords internal
#' @noRd
.spbp_anova_heading <- function(object) {
  model <- object$call$model
  if (is.null(model)) {
    model <- "ph"
  }
  c(
    "Analysis of Deviance Table",
    "",
    paste0(
      "Bernstein ", toupper(model),
      " model (MLE), n = ", object$n,
      ", degree = ", object$degree
    )
  )
}

#' Attach anova class and heading for print.anova
#' @keywords internal
#' @noRd
.spbp_as_anova_table <- function(tab, heading) {
  structure(
    tab,
    class = c("anova", "data.frame"),
    heading = heading
  )
}

#' Response expression as a character string
#' @keywords internal
#' @noRd
.spbp_formula_response_chr <- function(formula) {
  paste(deparse(stats::formula(formula)[[2L]], width.cutoff = 500L), collapse = " ")
}

#' Nested formulas from intercept-only through the full model
#' @keywords internal
#' @noRd
.spbp_nested_formulas <- function(formula) {
  formula <- stats::formula(formula)
  lbl <- attr(stats::terms(formula), "term.labels")
  resp <- .spbp_formula_response_chr(formula)
  out <- vector("list", length(lbl) + 1L)
  out[[1L]] <- stats::as.formula(paste(resp, "~ 1"))
  if (length(lbl) > 0L) {
    for (i in seq_along(lbl)) {
      out[[i + 1L]] <- stats::reformulate(lbl[seq_len(i)], response = resp)
    }
  }
  out
}

#' Refit an spbp model with a new formula using the original call settings
#' @keywords internal
#' @noRd
.spbp_safe_call_arg <- function(cl, name, default = NULL) {
  if (is.null(cl[[name]])) {
    return(default)
  }
  expr <- cl[[name]]
  if (is.name(expr) && identical(as.character(expr), name)) {
    return(default)
  }
  tryCatch(eval(expr, envir = parent.frame()), error = function(e) default)
}

#' @keywords internal
#' @noRd
.spbp_refit <- function(object, formula) {
  cl <- object$call
  fn <- cl[[1L]]
  model <- object$call$model
  if (is.null(model)) {
    if (identical(fn, quote(bpph))) {
      model <- "ph"
    } else if (identical(fn, quote(bppo))) {
      model <- "po"
    } else if (identical(fn, quote(bpaft))) {
      model <- "aft"
    } else {
      model <- .spbp_safe_call_arg(cl, "model", "ph")
    }
  }
  args <- list(
    formula = formula,
    degree = object$degree,
    data = object$data,
    approach = object$call$approach,
    model = model,
    scale = .spbp_safe_call_arg(cl, "scale", TRUE),
    verbose = .spbp_safe_call_arg(cl, "verbose", FALSE),
    chains = .spbp_safe_call_arg(cl, "chains", 4L)
  )
  for (nm in c("dist", "baseline", "priors", "init", "cores")) {
    val <- .spbp_safe_call_arg(cl, nm, default = "__missing__")
    if (!identical(val, "__missing__")) {
      args[[nm]] <- val
    }
  }
  fit <- do.call(spbp.default, args)
  if (!inherits(fit, "spbp")) {
    stop("Refit did not return an spbp object.", call. = FALSE)
  }
  fit
}

#' Sequential analysis-of-deviance table (anova single-model)
#' @keywords internal
#' @noRd
.spbp_anova_sequential <- function(object) {
  fmls <- .spbp_nested_formulas(object$formula)
  if (length(fmls) < 2L) {
    stop("Cannot perform LR test on a null model without covariates.", call. = FALSE)
  }

  fits <- lapply(fmls, function(fml) .spbp_refit(object, fml))
  logliks <- vapply(fits, function(f) as.numeric(logLik(f)), numeric(1))
  npars <- vapply(fits, .spbp_nparams, integer(1))
  resid_df <- vapply(fits, .spbp_resid_df, integer(1))
  term_labels <- attr(stats::terms(object$formula), "term.labels")
  row_names <- c("NULL", term_labels)
  n <- length(fits)

  deviance <- rep(NA_real_, n)
  df <- rep(NA_integer_, n)
  pval <- rep(NA_real_, n)

  for (i in seq_len(n - 1L)) {
    j <- i + 1L
    deviance[j] <- -2 * (logliks[i] - logliks[j])
    df[j] <- npars[j] - npars[i]
    if (is.finite(deviance[j]) && df[j] > 0L) {
      pval[j] <- stats::pchisq(deviance[j], df[j], lower.tail = FALSE)
    }
  }

  tab <- data.frame(
    Df = df,
    Deviance = deviance,
    `Resid. Df` = resid_df,
    `-2*LL` = -2 * logliks,
    `Pr(>Chi)` = pval,
    row.names = row_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  heading <- c(
    .spbp_anova_heading(object),
    "",
    "Terms added sequentially (first to last)"
  )
  .spbp_as_anova_table(tab, heading)
}
