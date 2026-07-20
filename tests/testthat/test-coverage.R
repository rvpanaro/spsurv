# Enforce minimum line coverage (covr, tests only).
# Run explicitly: NOT_CRAN=true SPSURV_RUN_COVERAGE=true Rscript -e 'testthat::test_file("tests/testthat/test-coverage.R")'

test_that("code coverage is maintained", {
  skip_if_not_installed("covr")
  skip_if_not(
    identical(Sys.getenv("SPSURV_RUN_COVERAGE"), "true"),
    "Set SPSURV_RUN_COVERAGE=true to run coverage check"
  )
  skip_on_cran()

  pkg_root <- testthat::test_path("..", "..")
  cov <- covr::package_coverage(
    path = pkg_root,
    type = "tests",
    line_exclusions = c("R/stanmodels.R")
  )
  pct <- covr::percent_coverage(cov)
  expect_gt(pct, 80)
})
