# spsurv 1.1.0

* **CRAN resubmission:** Maintainer contact updated to `rvpanaro@gmail.com`;
  `Authors@R` roles reviewed (single maintainer; thesis advisors listed as
  `ths`, with `rev` for code review feedback); package `Description` citation
  updated to `<doi:10.48550/arXiv.2003.10548>`; set
  `Config/rstantools/auto_config: FALSE` and ship committed Stan C++ exports
  (no install-time `configure`/`rstan_config()`, avoiding Windows R-devel
  `rstan.dll` load failures); added `RcppParallel` and synced
  `src/Makevars.win` TBB flags with Unix `Makevars`; restored `covr`,
  `rsample`, and `rsurv` in `Suggests` for optional tests referenced in
  `tests/`.
* **Inference robustness:** `.spbp_gamma_information_diagnostics()` no longer
  calls `qr()` on non-finite gamma information matrices; `vcov()`/`survfit()`
  return NA uncertainty with a warning when that block is unavailable.
* **tidymodels:** parsnip engines `spsurv` (MLE) and `spsurv_bayes` for `proportional_hazards()`, new `proportional_odds()`, and `survival_reg()`; `bp_survival_reg()` convenience constructor; `workflow` / `predict(type = "survival")` support.
* **tidybayes:** `as_draws_df.spbp()`, `spread_surv_draws.spbp()`, and S3 hooks for `tidy_draws` / `spread_draws` / `gather_draws` on Bayes fits; posterior draw indices (`.chain`, `.iteration`, `.draw`).
* **Prediction:** censored-style `predict.spbp()` types (`survival`, `time`, `linear_pred`); `augment.spbp()`.
* **Vignette:** `vignette("tidymodels")`.
* **Inference:** `vcov.spbp()` no longer applies a diagonal ridge (`polish` removed); beta standard errors are unchanged when the Bernstein block is stable and no longer inflated when it is ill-conditioned.
* **Summary:** `summary.spbp()` now includes a `coef_interval` component (Wald for MLE, HPD for Bayes) and prints it when `show_intervals = TRUE`.
* **Repository:** paper workflow moved from `inst/` to `paper/` (excluded from the package tarball).
* **Dependencies:** trimmed `Suggests` to packages used in vignettes, tests, or optional integrations; dev tools (`devtools`, `roxygen2`, `covr`) and unused tidymodels/rsurv entries removed (install transitively or only when needed locally).

# spsurv 1.0.2

* **Vignettes:** Expanded with EDA plots, visual checks, and decision guidance across all seven articles.
* **tidy/broom:** `tidy()`, `glance()`, `logLik()`, `AIC()`, `anova()`, `estimates()`, and `se()` methods for `spbp` objects; `bernstein()` baseline specification.
* **Prediction:** `predict()` and improved `survfit()` (tidy output, newdata covariates, Bayesian monotone bands).
* **Residuals:** `residuals()` supports martingale, deviance, and Cox-Snell types; martingale diagnostics fixed in examples.
* **Handlers:** Removed `handlers.R`; validation and setup logic inlined in `spbp.default`.
* **vcov:** Block-wise inversion of the Hessian (regression vs Bernstein polynomial blocks) in `vcov.spbp`.
* **Stan:** Slim `transformed parameters` (only `alpha` stored); pointwise `log_lik` moved to `generated quantities`; likelihood via `bp_pointwise_log_lik()` in the model block.
* **Stan (AFT):** No clamping of `u` to avoid boundary bias; minimum feasible divisor (`min_range`) only when range is degenerate.
* **MLE:** Initial values for the optimizer drawn from prior distributions; retries limited to 3 attempts.
* **`spbp` objects:** Component `degree` on fits and in `call$degree`; default Bernstein degree `ceiling(sqrt(n))` when omitted (documented in help, README, and vignettes).
* **`anova.spbp`:** Pairwise and sequential analysis-of-deviance tables aligned with `survival::survreg` (`Terms`, `Test`, `Df`, `Deviance`, etc.).
* **`bpph` / `bppo` / `bpaft`:** User-facing `match.call()` without unevaluated `dist`/`baseline` symbols; optional `degree` argument.
* **Documentation:** Expanded `spbp` / `spbp.default` help (`@details`, return components, vignette links); README usage and model-comparison examples updated.
* **Examples:** `data("veteran", package = "survival")` in docs and examples.
* **pkgdown:** Reference index updated to match current package; Bootstrap 5 / default template; logo size and alt-text; README aligned with site style.
* **Repository:** CRAN submission hygiene (`.gitignore`, `.Rbuildignore`, `LICENSE`).

# spsurv 1.0.1

* [pkgdown](https://github.com/r-lib/pkgdown) reference manual [https://rvpanaro.github.io/spsurv/reference/index.html](https://rvpanaro.github.io/spsurv/reference/index.html).

# spsurv 1.0.0

* This is the first release of spsurv.
