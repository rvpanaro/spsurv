# Tests for validation and setup logic inlined in spbp.default (R/spbp.R)
# Handlers were inlined; these tests exercise the same behavior via spbp().

library(spsurv)

context("spbp validation and setup (inlined handlers)")

# --- Prior parsing (formerly .handler3): default and custom priors via spbp ---
test_that("spbp bayes uses default priors when priors not specified", {
  skip_if_not_installed("rstan")
  fit <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    model = "ph",
    chains = 1,
    iter = 100,
    refresh = 0
  )
  expect_true(is.list(fit$standata))
  expect_named(fit$standata, c(
    "n", "m", "p", "tau", "approach", "rand", "M", "status", "id", "z",
    "time", "X", "g", "G", "P",
    "priordist_beta", "location_beta", "scale_beta",
    "priordist_gamma", "location_gamma", "scale_gamma",
    "priordist_frailty", "par1_frailty", "par2_frailty",
    "means", "sdv"
  ), ignore.order = TRUE)
})

test_that("spbp bayes accepts custom priors", {
  skip_if_not_installed("rstan")
  fit <- spbp(
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    model = "ph",
    priors = list(
      beta = c("normal(0,4)"),
      gamma = c("lognormal(0,4)"),
      frailty = c("gamma(0.01,0.01)")
    ),
    chains = 1,
    iter = 100,
    refresh = 0
  )
  expect_equal(fit$standata$location_beta, 0)
  expect_equal(fit$standata$scale_beta, 4)
  expect_equal(fit$standata$location_gamma, 0)
  expect_equal(fit$standata$scale_gamma, 4)
})

# --- Model frame / response validation (formerly .handler4) via spbp ---
test_that("spbp errors when response is not Surv", {
  expect_error(
    spbp(time ~ karno + factor(celltype), data = veteran),
    "Response must be a survival object"
  )
})

test_that("spbp errors when variable on both sides of formula", {
  expect_error(
    spbp(Surv(time, status) ~ status + karno, data = veteran),
    "both the left and right|variable appears"
  )
})

test_that("spbp runs for valid right-censored model", {
  expect_s3_class(
    spbp(Surv(time, status) ~ karno, data = veteran),
    "spbp"
  )
})

test_that("spbp errors on non-integer degree", {
  expect_error(
    spbp(Surv(time, status) ~ karno + factor(celltype), degree = 1 / 2, data = veteran),
    "Polynomial degree must be integer"
  )
})
