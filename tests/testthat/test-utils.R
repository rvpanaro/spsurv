# Unit tests for utils (R/utils.R): terms.inner, read_prior, hmean, LPML, DIC, WAIC

library(spsurv)

context("utils")

test_that("terms.inner extracts terms from formula", {
  f <- Surv(time, status) ~ karno + factor(celltype)
  out <- terms.inner(f)
  expect_true(is.character(out))
  expect_true(length(out) >= 2L)
})

test_that("read_prior parses prior string", {
  out <- read_prior("normal(0,5)")
  expect_equal(out[1], "normal")
  expect_equal(out[2], "0")
  expect_equal(out[3], "5")
})

test_that("read_prior parses lognormal", {
  out <- read_prior("lognormal(0,4)")
  expect_equal(out[1], "lognormal")
  expect_equal(length(out), 3L)
})

test_that("hmean computes harmonic mean", {
  x <- c(1, 2, 4)
  expect_equal(hmean(x), 1 / mean(1 / x))
  expect_equal(hmean(c(2, 2)), 2)
})

test_that("LPML returns length-2 vector", {
  loglik <- matrix(rnorm(20, -2, 0.5), nrow = 10, ncol = 2)
  out <- LPML(loglik)
  expect_true(is.vector(out))
  expect_equal(length(out), 2L)
  expect_true(is.numeric(out))
})

test_that("DIC returns matrix with 2 columns", {
  loglik <- matrix(rnorm(20, -2, 0.5), nrow = 10, ncol = 2)
  out <- DIC(loglik)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(nrow(out), 1L)
})

test_that("WAIC returns matrix with 2 columns", {
  loglik <- matrix(rnorm(20, -2, 0.5), nrow = 10, ncol = 2)
  out <- WAIC(loglik)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(nrow(out), 1L)
})

# .survfit_confint (lines 71-120 in utils.R): lines 75-89 = selow, logse, conf.type "plain"
confint_surv <- function(...) spsurv:::.survfit_confint(...)

test_that(".survfit_confint with selow missing uses scale 1 (line 73-74)", {
  out <- confint_surv(p = 0.8, se = 0.1, conf.type = "plain", conf.int = 0.95)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower <= 0.8)
  expect_true(out$upper >= 0.8)
})

test_that(".survfit_confint with selow = 0 uses scale 1 (line 76)", {
  out <- confint_surv(p = 0.8, se = 0.1, conf.type = "plain", conf.int = 0.95, selow = 0)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower <= 0.8)
  expect_true(out$upper >= 0.8)
})

test_that(".survfit_confint with selow non-zero uses scale = selow/se (lines 75-76)", {
  out <- confint_surv(p = 0.8, se = 0.1, conf.type = "plain", conf.int = 0.95, selow = 0.05)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower <= 0.8)
  expect_true(out$upper >= 0.8)
})

test_that(".survfit_confint with logse = FALSE rescales se (lines 77-79)", {
  out <- confint_surv(p = 0.8, se = 0.1, logse = FALSE, conf.type = "plain", conf.int = 0.95)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower <= 0.8)
  expect_true(out$upper >= 0.8)
})

test_that(".survfit_confint conf.type plain with ulimit TRUE (lines 80-84)", {
  out <- confint_surv(p = 0.5, se = 0.1, conf.type = "plain", conf.int = 0.95, ulimit = TRUE)
  expect_true(all(out$lower >= 0))
  expect_true(all(out$upper <= 1))
  expect_true(out$lower < 0.5)
  expect_true(out$upper > 0.5)
})

test_that(".survfit_confint conf.type plain with ulimit FALSE (lines 85-86)", {
  out <- confint_surv(p = 0.5, se = 0.1, conf.type = "plain", conf.int = 0.95, ulimit = FALSE)
  expect_true(out$lower >= 0)
  expect_true(out$upper > 0.5)
  expect_true(out$lower < 0.5)
})

# Lines 97-120: conf.type "log-log", "logit", "arcsin", invalid type
test_that(".survfit_confint conf.type log-log (lines 99-104)", {
  out <- confint_surv(p = 0.5, se = 0.1, conf.type = "log-log", conf.int = 0.95)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower < 0.5)
  expect_true(out$upper > 0.5)
  expect_true(out$lower >= 0)
  expect_true(out$upper <= 1)
})

test_that(".survfit_confint conf.type logit (lines 105-110)", {
  out <- confint_surv(p = 0.5, se = 0.1, conf.type = "logit", conf.int = 0.95)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower < 0.5)
  expect_true(out$upper > 0.5)
  expect_true(out$lower > 0)
  expect_true(out$upper < 1)
})

test_that(".survfit_confint conf.type arcsin (lines 111-117)", {
  out <- confint_surv(p = 0.5, se = 0.1, conf.type = "arcsin", conf.int = 0.95)
  expect_named(out, c("lower", "upper"))
  expect_true(out$lower < 0.5)
  expect_true(out$upper > 0.5)
  expect_true(out$lower >= 0)
  expect_true(out$upper <= 1)
})

test_that(".survfit_confint invalid conf.type stops (lines 118-119)", {
  expect_error(
    confint_surv(p = 0.5, se = 0.1, conf.type = "invalid", conf.int = 0.95),
    "invalid conf.int type"
  )
})
