# Changelog

## spsurv 1.1.0

- **CRAN resubmission:** Maintainer contact updated to
  `rvpanaro@gmail.com`; `Authors@R` roles reviewed (single maintainer;
  thesis advisors listed as `ths`, with `rev` for code review feedback);
  package `Description` citation updated to
  `<doi:10.48550/arXiv.2003.10548>`; set
  `Config/rstantools/auto_config: FALSE` and ship committed Stan C++
  exports (no install-time `configure`/`rstan_config()`, avoiding
  Windows R-devel `rstan.dll` load failures); added `RcppParallel` and
  synced `src/Makevars.win` TBB flags with Unix `Makevars`; restored
  `covr`, `rsample`, and `rsurv` in `Suggests` for optional tests
  referenced in `tests/`.
- **Inference robustness:** `.spbp_gamma_information_diagnostics()` no
  longer calls [`qr()`](https://rdrr.io/r/base/qr.html) on non-finite
  gamma information matrices;
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
  return NA uncertainty with a warning when that block is unavailable.
- **tidymodels:** parsnip engines `spsurv` (MLE) and `spsurv_bayes` for
  [`proportional_hazards()`](https://parsnip.tidymodels.org/reference/proportional_hazards.html),
  new
  [`proportional_odds()`](https://rvpanaro.github.io/spsurv/reference/proportional_odds.md),
  and
  [`survival_reg()`](https://parsnip.tidymodels.org/reference/survival_reg.html);
  [`bp_survival_reg()`](https://rvpanaro.github.io/spsurv/reference/bp_survival_reg.md)
  convenience constructor; `workflow` / `predict(type = "survival")`
  support.
- **tidybayes:**
  [`as_draws_df.spbp()`](https://rvpanaro.github.io/spsurv/reference/as_draws_df.spbp.md),
  [`spread_surv_draws.spbp()`](https://rvpanaro.github.io/spsurv/reference/spread_surv_draws.spbp.md),
  and S3 hooks for `tidy_draws` / `spread_draws` / `gather_draws` on
  Bayes fits; posterior draw indices (`.chain`, `.iteration`, `.draw`).
- **Prediction:** censored-style
  [`predict.spbp()`](https://rvpanaro.github.io/spsurv/reference/predict.spbp.md)
  types (`survival`, `time`, `linear_pred`);
  [`augment.spbp()`](https://rvpanaro.github.io/spsurv/reference/augment.spbp.md).
- **Vignette:**
  [`vignette("tidymodels")`](https://rvpanaro.github.io/spsurv/articles/tidymodels.md).
- **Inference:**
  [`vcov.spbp()`](https://rvpanaro.github.io/spsurv/reference/vcov.spbp.md)
  no longer applies a diagonal ridge (`polish` removed); beta standard
  errors are unchanged when the Bernstein block is stable and no longer
  inflated when it is ill-conditioned.
- **Summary:**
  [`summary.spbp()`](https://rvpanaro.github.io/spsurv/reference/summary.spbp.md)
  now includes a `coef_interval` component (Wald for MLE, HPD for Bayes)
  and prints it when `show_intervals = TRUE`.
- **Repository:** paper workflow moved from `inst/` to `paper/`
  (excluded from the package tarball).
- **Dependencies:** trimmed `Suggests` to packages used in vignettes,
  tests, or optional integrations; dev tools (`devtools`, `roxygen2`,
  `covr`) and unused tidymodels/rsurv entries removed (install
  transitively or only when needed locally).

## spsurv 1.0.2

- **Vignettes:** Expanded with EDA plots, visual checks, and decision
  guidance across all seven articles.
- **tidy/broom:**
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`anova()`](https://rdrr.io/r/stats/anova.html),
  [`estimates()`](https://rvpanaro.github.io/spsurv/reference/estimates.md),
  and [`se()`](https://rvpanaro.github.io/spsurv/reference/se.md)
  methods for `spbp` objects;
  [`bernstein()`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)
  baseline specification.
- **Prediction:** [`predict()`](https://rdrr.io/r/stats/predict.html)
  and improved
  [`survfit()`](https://rdrr.io/pkg/survival/man/survfit.html) (tidy
  output, newdata covariates, Bayesian monotone bands).
- **Residuals:** [`residuals()`](https://rdrr.io/r/stats/residuals.html)
  supports martingale, deviance, and Cox-Snell types; martingale
  diagnostics fixed in examples.
- **Handlers:** Removed `handlers.R`; validation and setup logic inlined
  in `spbp.default`.
- **vcov:** Block-wise inversion of the Hessian (regression vs Bernstein
  polynomial blocks) in `vcov.spbp`.
- **Stan:** Slim `transformed parameters` (only `alpha` stored);
  pointwise `log_lik` moved to `generated quantities`; likelihood via
  `bp_pointwise_log_lik()` in the model block.
- **Stan (AFT):** No clamping of `u` to avoid boundary bias; minimum
  feasible divisor (`min_range`) only when range is degenerate.
- **MLE:** Initial values for the optimizer drawn from prior
  distributions; retries limited to 3 attempts.
- **`spbp` objects:** Component `degree` on fits and in `call$degree`;
  default Bernstein degree `ceiling(sqrt(n))` when omitted (documented
  in help, README, and vignettes).
- **`anova.spbp`:** Pairwise and sequential analysis-of-deviance tables
  aligned with
  [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  (`Terms`, `Test`, `Df`, `Deviance`, etc.).
- **`bpph` / `bppo` / `bpaft`:** User-facing
  [`match.call()`](https://rdrr.io/r/base/match.call.html) without
  unevaluated `dist`/`baseline` symbols; optional `degree` argument.
- **Documentation:** Expanded `spbp` / `spbp.default` help (`@details`,
  return components, vignette links); README usage and model-comparison
  examples updated.
- **Examples:** `data("veteran", package = "survival")` in docs and
  examples.
- **pkgdown:** Reference index updated to match current package;
  Bootstrap 5 / default template; logo size and alt-text; README aligned
  with site style.
- **Repository:** CRAN submission hygiene (`.gitignore`,
  `.Rbuildignore`, `LICENSE`).

## spsurv 1.0.1

- [pkgdown](https://github.com/r-lib/pkgdown) reference manual
  <https://rvpanaro.github.io/spsurv/reference/index.html>.

## spsurv 1.0.0

CRAN release: 2020-03-31

- This is the first release of spsurv.
