## Test environments

* local: R 4.5.x (macOS, aarch64), `R CMD check --as-cran`
* rhub: run `rhub::check_for_cran()` before upload
* win-builder: run before upload

## testthat (local, spsurv 1.1.0)

```
576 tests, 0 failures, 39 warnings (Stan / ill-conditioned gamma diagnostics)
```

Run: `testthat::test_dir("tests/testthat")` after `pkgload::load_all()`.

## Code coverage (local, covr)

```
Overall line coverage: 81.18%  (type = "tests"; excludes R/stanmodels.R)
```

| File | Coverage |
|------|----------|
| R/vcov.R | 95.9% |
| R/summary.spbp.R | 98.6% |
| R/survfit.R | 85.2% |
| R/print.summary.spbp.R | 80.4% |
| R/utils.R | 78.5% |

Run: `covr::package_coverage(type = "tests", line_exclusions = "R/stanmodels.R")`.
Threshold enforced in `tests/testthat/test-coverage.R` (> 80%; opt-in via `SPSURV_RUN_COVERAGE=true`, skipped on CRAN).

## R CMD check results (spsurv 1.1.0)

```
Status: (run R CMD check --as-cran before upload)
```

## revdepcheck

No reverse dependencies known at time of release.

## Submission notes

* Minor release 1.1.0: tidymodels integration (parsnip engines `spsurv` and
  `spsurv_bayes`, `proportional_odds()`, `survival_reg()`, `bp_survival_reg()`),
  tidybayes/posterior draw helpers, censored-style `predict()` types
  (`survival`, `time`, `linear_pred`), and `augment.spbp()`.
* New vignette `vignette("tidymodels")`.
* Installed size is dominated by compiled Stan code in `libs/` (expected).
* Bayesian tests use reduced `iter`/`chains` for speed.
