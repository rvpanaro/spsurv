utils::globalVariables(c("means", "sdv", "p", "rand", "signif.stars", "x", "y"))

#' @importFrom stats AIC confint predict residuals
NULL

.onLoad <- function(libname, pkgname) {
  spsurv_register_parsnip()
  spsurv_register_tidybayes_s3()
}

#' @keywords internal
#' @noRd
spsurv_register_tidybayes_s3 <- function() {
  if (requireNamespace("posterior", quietly = TRUE)) {
    registerS3method(
      "as_draws_df", "spbp", as_draws_df.spbp,
      envir = asNamespace("posterior")
    )
    registerS3method(
      "as_draws", "spbp",
      function(x, ...) as_draws_df.spbp(x, ...),
      envir = asNamespace("posterior")
    )
  }
  if (requireNamespace("tidybayes", quietly = TRUE)) {
    tb_env <- asNamespace("tidybayes")
    registerS3method("tidy_draws", "spbp", tidy_draws.spbp, envir = tb_env)
    registerS3method("spread_draws", "spbp", spread_draws.spbp, envir = tb_env)
    registerS3method("gather_draws", "spbp", gather_draws.spbp, envir = tb_env)
  }
  invisible(NULL)
}
