#' Bernstein polynomial baseline specification
#'
#' Helper for specifying a Bernstein polynomial baseline via the
#' \code{dist} / \code{baseline} arguments accepted by \code{\link{spbp}}.
#'
#' @param m Bernstein polynomial degree (number of basis coefficients).
#' @return A list with components \code{baseline} and \code{m}.
#' @export
#' @examples
#' bernstein(5)
bernstein <- function(m = NULL) {
  list(baseline = "bernstein", m = m)
}
