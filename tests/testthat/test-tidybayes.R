# tidybayes / posterior draws tests

test_that("as_draws_df.spbp includes chain metadata", {
  skip_if_not_installed("posterior")
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 20L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  dr <- spsurv::as_draws_df.spbp(fit)
  expect_true(all(c(".chain", ".iteration", ".draw") %in% colnames(dr)))
  beta_cols <- grep("^beta\\[", colnames(dr), value = TRUE)
  expect_length(beta_cols, 1L)
})

test_that("spread_draws.spbp works for beta", {
  skip_if_not_installed("posterior")
  skip_if_not_installed("tidybayes")
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 20L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  long <- tidybayes::spread_draws(fit, `beta[karno]`)
  expect_true(all(c(".draw", "beta[karno]") %in% names(long)))
  expect_equal(nrow(long), nrow(fit$posterior$beta))
})

test_that("spread_surv_draws returns long survival draws", {
  fit <- quick_bayes(
    bpph,
    Surv(time, status) ~ karno,
    data = veteran,
    approach = "bayes",
    iter = 10L,
    chains = 1L,
    cores = 1L,
    init = 0
  )
  times <- c(50, 100)
  nd <- veteran[1:2, , drop = FALSE]
  long <- spsurv::spread_surv_draws.spbp(fit, times = times, newdata = nd)
  expect_equal(nrow(long), nrow(fit$posterior$beta) * length(times) * nrow(nd))
  expect_true(all(c(".draw", "time", "surv") %in% names(long)))
})
