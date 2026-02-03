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
