# parsnip / workflow integration tests

skip_if_not_installed <- function(pkg) {
  testthat::skip_if_not_installed(pkg)
}

test_that("spsurv registers parsnip engines on load", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  eng <- parsnip::show_engines("proportional_hazards")
  expect_true("spsurv" %in% eng$engine)
  expect_true("spsurv_bayes" %in% eng$engine)
})

test_that("parsnip MLE fit and censored predict", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  spec <- parsnip::proportional_hazards() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  fit <- parsnip::fit(spec, Surv(time, status) ~ karno, data = veteran)
  pr <- predict(fit, veteran[1:2, ], type = "survival", eval_time = c(100, 200))
  expect_true(".pred" %in% names(pr))
})

test_that("proportional_odds model fits via parsnip", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  spec <- spsurv::proportional_odds() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  fit <- parsnip::fit(spec, Surv(time, status) ~ karno, data = veteran)
  expect_s3_class(fit$fit, "spbp")
})

test_that("workflow with spsurv engine", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  skip_if_not_installed("workflows")
  spec <- parsnip::proportional_hazards() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  wf <- workflows::workflow() |>
    workflows::add_formula(Surv(time, status) ~ karno + celltype) |>
    workflows::add_model(spec)
  fit <- workflows::fit(wf, data = veteran)
  pr <- predict(fit, veteran[1:2, ], type = "survival", eval_time = 100)
  expect_true(".pred" %in% names(pr))
})

test_that("workflow resampling smoke test with spsurv engine", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  skip_if_not_installed("workflows")
  skip_if_not_installed("rsample")
  set.seed(1)
  split <- rsample::initial_split(veteran)
  train <- rsample::training(split)
  spec <- parsnip::proportional_hazards() |>
    parsnip::set_engine("spsurv", degree = 4L, scale = FALSE, init = 0)
  wf <- workflows::workflow() |>
    workflows::add_formula(Surv(time, status) ~ karno) |>
    workflows::add_model(spec)
  fit <- workflows::fit(wf, data = train)
  expect_s3_class(parsnip::extract_fit_engine(fit), "spbp")
})

test_that("bayes parsnip engine fits and predicts", {
  skip_if_not_installed("parsnip")
  skip_if_not_installed("censored")
  spec <- parsnip::proportional_hazards() |>
    parsnip::set_engine(
      "spsurv_bayes",
      degree = 3L,
      scale = FALSE,
      iter = 50L,
      chains = 1L,
      cores = 1L,
      init = 0
    )
  fit <- suppressWarnings(parsnip::fit(spec, Surv(time, status) ~ karno, data = veteran))
  pr <- predict(fit, veteran[1:2, ], type = "linear_pred")
  expect_true(".pred_linear_pred" %in% names(pr))
})
