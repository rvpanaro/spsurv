# Unit tests for stanmodels (R/stanmodels.R)

context("stanmodels")

test_that("stanmodels object exists after load", {
  expect_true(exists("stanmodels", where = asNamespace("spsurv"), inherits = FALSE))
})

test_that("stanmodels includes spbp", {
  sm <- get("stanmodels", envir = asNamespace("spsurv"))
  has_spbp <- if (is.character(sm)) "spbp" %in% sm else "spbp" %in% names(sm)
  expect_true(has_spbp)
})
