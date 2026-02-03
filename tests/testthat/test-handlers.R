# Unit tests for .handler1, .handler2, .handler3, .handler4 (R/handlers.R)

library(spsurv)

# Call internal handlers explicitly (they now take arguments and return values)
handler1 <- function(...) spsurv:::.handler1(...)
handler2 <- function(...) spsurv:::.handler2(...)
handler3 <- function(...) spsurv:::.handler3(...)
handler4 <- function(...) spsurv:::.handler4(...)

context("handlers")

# --- .handler1(Call) ---
test_that(".handler1 returns indices for formula and data in Call", {
  # Simulate Call as in spbp.default: match.call() inside a function
  dummy <- function(formula, data, ...) match.call()
  Call <- dummy(Surv(time, status) ~ karno, data = veteran)
  aux <- handler1(Call)
  expect_equal(length(aux), 2L)
  expect_equal(names(Call)[aux], c("formula", "data"))
})

test_that(".handler1 returns nomatch=0 for missing names", {
  Call <- as.call(list(quote(f), x = 1))
  aux <- handler1(Call)
  expect_equal(aux, c(0L, 0L))
})

# --- .handler2(temp, formula) ---
test_that(".handler2 returns rand=0 when no frailty in formula", {
  f <- Surv(time, status) ~ karno + factor(celltype)
  temp <- list(formula = terms(f, data = veteran))
  attr(temp$formula, "specials") <- list()
  out <- handler2(temp, f)
  expect_named(out, c("rand", "id", "frailty_idx"))
  expect_equal(out$rand, 0)
  expect_null(out$id)
  expect_null(out$frailty_idx)
})

test_that(".handler2 returns list of correct structure", {
  f <- Surv(time, status) ~ karno
  temp <- list(formula = terms(f, data = veteran))
  attr(temp$formula, "specials") <- list()
  out <- handler2(temp, f)
  expect_type(out, "list")
  expect_equal(names(out), c("rand", "id", "frailty_idx"))
})

# --- .handler3(priors) ---
test_that(".handler3 returns default priors when priors empty", {
  priors <- list(beta = list(), gamma = list(), frailty = list())
  out <- handler3(priors)
  expect_named(out, c(
    "priordist_beta", "location_beta", "scale_beta",
    "priordist_gamma", "location_gamma", "scale_gamma",
    "priordist_frailty", "par1_frailty", "par2_frailty"
  ))
  expect_equal(out$priordist_beta, "normal")
  expect_equal(out$location_beta, "0")
  expect_equal(out$scale_beta, "5")
  expect_equal(out$priordist_gamma, "lognormal")
  expect_equal(out$priordist_frailty, "gamma")
  expect_equal(out$par1_frailty, "1")
  expect_equal(out$par2_frailty, "1")
})

test_that(".handler3 parses custom beta and gamma priors", {
  priors <- list(
    beta = c("normal(0,4)"),
    gamma = c("lognormal(0,4)"),
    frailty = c("gamma(0.01,0.01)")
  )
  out <- handler3(priors)
  expect_equal(out$priordist_beta, "normal")
  expect_equal(out$location_beta, "0")
  expect_equal(out$scale_beta, "4")
  expect_equal(out$priordist_gamma, "lognormal")
  expect_equal(out$location_gamma, "0")
  expect_equal(out$scale_gamma, "4")
  expect_equal(out$priordist_frailty, "gamma")
  expect_equal(out$par1_frailty, "0.01")
  expect_equal(out$par2_frailty, "0.01")
})

# --- .handler4(mf, Y, type, Terms, formula) ---
test_that(".handler4 stops when response is not Surv", {
  mf <- model.frame(Surv(time, status) ~ karno, data = veteran)
  Y <- model.extract(mf, "response")
  Terms <- terms(mf)
  f <- Surv(time, status) ~ karno
  # Pass non-Surv as Y
  expect_error(handler4(mf, veteran$time, "right", Terms, f), "Response must be a survival object")
})

test_that(".handler4 stops when type is not right or counting", {
  mf <- model.frame(Surv(time, status) ~ karno, data = veteran)
  Y <- model.extract(mf, "response")
  Terms <- terms(mf)
  f <- Surv(time, status) ~ karno
  expect_error(
    handler4(mf, Y, "left", Terms, f),
    "spsurv doesn't support"
  )
})

test_that(".handler4 stops when variable on both sides of formula", {
  mf <- model.frame(Surv(time, status) ~ status + karno, data = veteran)
  Y <- model.extract(mf, "response")
  Terms <- terms(mf)
  f <- Surv(time, status) ~ status + karno
  expect_error(
    handler4(mf, Y, "right", Terms, f),
    "both the left and right|variable appears"
  )
})

test_that(".handler4 runs without error for valid right-censored model", {
  mf <- model.frame(Surv(time, status) ~ karno, data = veteran)
  Y <- model.extract(mf, "response")
  Terms <- terms(mf)
  f <- Surv(time, status) ~ karno
  expect_silent(handler4(mf, Y, "right", Terms, f))
})

# --- Integration: spbp still triggers same validations (handler4) ---
test_that("spbp errors when variable on both sides of formula", {
  expect_error(
    spbp(Surv(time, status) ~ status + karno, data = veteran),
    "both the left and right|variable appears"
  )
})

# test_that("spbp errors on invalid formula (not Surv)", {
#   expect_error(
#     spbp(time ~ karno + factor(celltype), data = veteran),
#     "Response must be a survival object|formula"
#   )
# })

# test_that("spbp errors on unsupported Surv type", {
#   time2 <- veteran$time + 1
#   expect_error(
#     spbp(Surv(time = veteran$time, time2 = time2, status) ~ karno, data = veteran),
#     "doesn't support|counting|left"
#   )
# })

test_that("spbp errors on non-integer degree", {
  expect_error(
    spbp(Surv(time, status) ~ karno + factor(celltype), degree = 1 / 2, data = veteran),
    "Polynomial degree must be integer"
  )
})
