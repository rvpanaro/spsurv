#' The 'spsurv' package.
#'
#' @description A set of reliable routines to ease semiparametric
#' survival regression modeling based on Bernstein Polynomials.
#' \code{spsurv::spbp()} includes Proportional Hazards, Proportional Odds and
#' Accelerated Failure Time frameworks for right censored data, it
#' also includes bayesian frailty modeling.
#'
#' @docType package
#' @name spsurv-package
#' @aliases spsurv
#' @useDynLib spsurv, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import survival
#' @importFrom rstan stan sampling optimizing
#'
#' @references
#'
#' Collett, D. (2003). Modelling survival data in medical research. Chapman and Hall/CRC.
#'
#' Colosimo, E., & Giolo S. R. (2006). Análise de sobrevivência aplicada. Editora Edgard Blucher Ltda.
#'
#' Ibrahim, J. G., Chen, M. H., & Sinha, D. (2001). Bayesian Survival Analysis. Wiley StatsRef: Statistics Reference Online.
#'
#' Klein, J. P., & Moeschberger, M. L. (1997). Survival analysis: techniques for censored and truncated data. Springer Science & Business Media.
#'
#' Lorentz, G. G. (1953). Bernstein polynomials. American Mathematical Society.
#'
#' Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for right-censored data using Bernstein polynomials. Computational Statistics & Data Analysis, 56(3), 559-573.
#'
NULL
