#' Glance at a fitted spbp model
#'
#' @description
#' Return a one-row broom/generics-style summary of model-level statistics.
#'
#' @param x A fitted \code{"spbp"} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with model-level statistics. Columns
#'   that are entirely \code{NA} for the fitted object are omitted.
#' @export
#' @method glance spbp
#' @importFrom generics glance
glance.spbp <- function(x, ...) {
  nevent <- if (!is.null(x$nevent)) x$nevent else NA_integer_

  if (is.null(x$coefficients) && is.null(x$bp.param)) {
    out <- data.frame(
      n = x$n,
      nevent = nevent,
      logLik = NA_real_,
      approach = x$call$approach,
      model = x$call$model,
      df = 0L,
      statistic = NA_real_,
      p.value = NA_real_,
      rsq = NA_real_,
      max.rsq = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      waic = NA_real_,
      dic = NA_real_,
      lpml = NA_real_,
      stringsAsFactors = FALSE
    )
    return(.spbp_drop_all_na_cols(out))
  }

  sm <- summary(x)

  if (identical(x$call$approach, "mle")) {
    ll_obj <- logLik(x)
    logLik_val <- if (!is.null(ll_obj)) as.numeric(ll_obj) else NA_real_
    df <- if (!is.null(ll_obj)) attr(ll_obj, "df") else .spbp_nparams(x)
    statistic <- if (!is.null(x$coefficients)) unname(sm$logtest["test"]) else NA_real_
    p.value <- if (!is.null(x$coefficients)) unname(sm$logtest["pvalue"]) else NA_real_
    rsq <- if (!is.null(x$coefficients)) unname(sm$rsq["rsq"]) else NA_real_
    max.rsq <- if (!is.null(x$coefficients)) unname(sm$rsq["maxrsq"]) else NA_real_
    AIC <- if (is.finite(logLik_val)) -2 * logLik_val + 2 * df else NA_real_
    BIC <- if (is.finite(logLik_val)) -2 * logLik_val + log(x$n) * df else NA_real_
    waic <- NA_real_
    dic <- NA_real_
    lpml <- NA_real_
  } else {
    logLik_val <- if (!is.null(x$loglik)) sum(x$loglik) else NA_real_
    df <- if (is.null(x$coefficients)) {
      length(x$bp.param)
    } else {
      length(x$coefficients)
    }
    statistic <- NA_real_
    p.value <- NA_real_
    rsq <- NA_real_
    max.rsq <- NA_real_
    AIC <- NA_real_
    BIC <- NA_real_
    waic <- if (!is.null(sm$waic)) sm$waic else NA_real_
    dic <- if (!is.null(sm$dic)) sm$dic else NA_real_
    lpml <- if (!is.null(sm$lpml)) sm$lpml else NA_real_
  }

  out <- data.frame(
    n = x$n,
    nevent = nevent,
    logLik = logLik_val,
    approach = x$call$approach,
    model = x$call$model,
    df = df,
    statistic = statistic,
    p.value = p.value,
    rsq = rsq,
    max.rsq = max.rsq,
    AIC = AIC,
    BIC = BIC,
    waic = waic,
    dic = dic,
    lpml = lpml,
    stringsAsFactors = FALSE
  )
  .spbp_drop_all_na_cols(out)
}
