#' Power basis polynomials calculations
#' @export
#' @param degree Bernstein polynomial degree
#' @return A list containing matrices g and G corresponding BP basis and corresponding tau value used to compute them.

pw.basis <- function(degree) {
  if (degree < 0) {
    stop("polynomial degree must be positive.")
  }
  if (!(degree %% 1 == 0)) {
    stop("polynomial degree must be integer.")
  }
  sapply(1:degree, function(k) {
    sapply(1:degree, function(j) {
      ((-1)^(j - k)) * choose(degree - k, j - k) / beta(k, degree - k + 1)
    })
  }) |> as.matrix()

}
