#' Proportional odds regression model specification
#'
#' @description
#' Proportional odds survival regression with a Bernstein polynomial baseline
#' (\code{\link{bppo}}). This model type is registered by \pkg{spsurv} for use
#' with tidymodels; it is not part of the \pkg{censored} package.
#'
#' @inheritParams parsnip::proportional_hazards
#' @return A model specification.
#' @export
proportional_odds <- function(mode = "censored regression", engine = "spsurv") {
  parsnip::new_model_spec(
    "proportional_odds",
    args = list(),
    eng_args = NULL,
    mode = mode,
    user_specified_mode = !missing(mode),
    method = NULL,
    engine = engine,
    user_specified_engine = !missing(engine)
  )
}

#' Bernstein survival regression specification
#'
#' @description Convenience constructor mapping \code{family} to the appropriate
#'   parsnip model (\code{proportional_hazards}, \code{proportional_odds}, or
#'   \code{survival_reg}).
#'
#' @param family Model family: \code{"ph"}, \code{"po"}, or \code{"aft"}.
#' @param mode Model mode (default \code{"censored regression"}).
#' @param engine parsnip engine (default \code{"spsurv"} for MLE).
#' @return A parsnip model specification.
#' @export
bp_survival_reg <- function(family = c("ph", "po", "aft"),
                            mode = "censored regression",
                            engine = "spsurv") {
  family <- match.arg(family)
  switch(
    family,
    ph = parsnip::proportional_hazards(mode = mode, engine = engine),
    po = proportional_odds(mode = mode, engine = engine),
    aft = parsnip::survival_reg(mode = mode, engine = engine)
  )
}

#' Register spsurv engines with parsnip
#'
#' @description Called on package load when \pkg{parsnip} is available.
#' @keywords internal
#' @noRd
spsurv_register_parsnip <- function() {
  if (!requireNamespace("parsnip", quietly = TRUE)) {
    return(invisible(NULL))
  }

  po_exists <- tryCatch({
    parsnip::show_engines("proportional_odds")
    TRUE
  }, error = function(e) FALSE)
  if (!po_exists) {
    tryCatch({
      parsnip::set_new_model("proportional_odds")
      parsnip::set_model_mode("proportional_odds", mode = "censored regression")
    }, error = function(e) invisible(NULL))
  }

  engines <- list(
    spsurv = list(approach = "mle"),
    spsurv_bayes = list(approach = "bayes", iter = 1000L, chains = 2L, cores = 1L)
  )

  for (eng in names(engines)) {
    .spsurv_register_ph_engine(eng, engines[[eng]])
    .spsurv_register_po_engine(eng, engines[[eng]])
    .spsurv_register_aft_engine(eng, engines[[eng]])
  }

  invisible(NULL)
}

#' @keywords internal
#' @noRd
.spsurv_register_common_fit <- function(model, eng, fit_fun, defaults) {
  parsnip::set_model_engine(model, mode = "censored regression", eng = eng)
  parsnip::set_dependency(model, eng = eng, pkg = "spsurv", mode = "censored regression")
  parsnip::set_fit(
    model = model,
    eng = eng,
    mode = "censored regression",
    value = list(
      interface = "formula",
      protect = c("formula", "data", "weights"),
      func = c(pkg = "spsurv", fun = fit_fun),
      defaults = defaults
    )
  )
  parsnip::set_encoding(
    model = model,
    eng = eng,
    mode = "censored regression",
    options = list(
      predictor_indicators = "traditional",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )
}

#' @keywords internal
#' @noRd
.spsurv_register_common_pred <- function(model, eng) {
  for (pred_type in c("survival", "time", "linear_pred")) {
    fun <- switch(
      pred_type,
      survival = "spsurv_pred_survival",
      time = "spsurv_pred_time",
      linear_pred = "spsurv_pred_linear_pred"
    )
    args <- list(
      object = quote(object),
      new_data = quote(new_data)
    )
    if (identical(pred_type, "survival")) {
      args$eval_time <- quote(eval_time)
    }
    parsnip::set_pred(
      model = model,
      eng = eng,
      mode = "censored regression",
      type = pred_type,
      value = list(
        pre = NULL,
        post = NULL,
        func = c(pkg = "spsurv", fun = fun),
        args = args
      )
    )
  }
}

#' @keywords internal
#' @noRd
.spsurv_register_ph_engine <- function(eng, extra_defaults) {
  defaults <- utils::modifyList(
    list(scale = FALSE, init = 0, verbose = FALSE),
    extra_defaults
  )
  .spsurv_register_common_fit(
    "proportional_hazards",
    eng,
    "spsurv_fit_proportional_hazards",
    defaults
  )
  .spsurv_register_common_pred("proportional_hazards", eng)
}

#' @keywords internal
#' @noRd
.spsurv_register_po_engine <- function(eng, extra_defaults) {
  defaults <- utils::modifyList(
    list(scale = FALSE, init = 0, verbose = FALSE),
    extra_defaults
  )
  .spsurv_register_common_fit(
    "proportional_odds",
    eng,
    "spsurv_fit_proportional_odds",
    defaults
  )
  .spsurv_register_common_pred("proportional_odds", eng)
}

#' @keywords internal
#' @noRd
.spsurv_register_aft_engine <- function(eng, extra_defaults) {
  defaults <- utils::modifyList(
    list(scale = FALSE, init = 0, verbose = FALSE),
    extra_defaults
  )
  .spsurv_register_common_fit(
    "survival_reg",
    eng,
    "spsurv_fit_survival_reg",
    defaults
  )
  .spsurv_register_common_pred("survival_reg", eng)
}
