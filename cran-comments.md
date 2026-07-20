## Test environments

* local: R 4.5.x (macOS, aarch64), `R CMD check --as-cran`
* rhub: run `rhub::check_for_cran()` before upload
* win-builder: run before upload

## testthat (local, spsurv 1.1.0, 2026-07-20)

```
devtools::test(): 576 passed, 0 failed, 39 warnings
(Stan / ill-conditioned gamma diagnostics)
```

## Code coverage (local, covr, 2026-07-20)

```
Overall line coverage: 81.18%  (type = "tests"; excludes R/stanmodels.R)
R/ mean line coverage: 81.92%
```

| File | Coverage |
|------|----------|
| R/vcov.R | 95.9% |
| R/summary.spbp.R | 96.4% |
| R/survfit.R | 85.5% |
| R/print.summary.spbp.R | 89.0% |
| R/utils.R | 80.9% |

Run: `NOT_CRAN=true SPSURV_RUN_COVERAGE=true Rscript -e 'testthat::test_file("tests/testthat/test-coverage.R")'`.
Install `covr` locally or in CI when running that opt-in check (`covr` is no longer in `Suggests`).
Threshold enforced in `tests/testthat/test-coverage.R` (> 80%; opt-in via `SPSURV_RUN_COVERAGE=true`, skipped on CRAN).

## Suggests (spsurv 1.1.0)

Required for CRAN check vignettes/tests: `testthat`, `knitr`, `rmarkdown`, `ggplot2`,
`KMsurv`, `parsnip`, `censored`, `workflows`, `tidybayes`, `posterior`.

Removed from `Suggests` (dev-only, transitive, or optional skipped tests):
`devtools`, `roxygen2`, `covr`, `ggsurvfit`, `recipes`, `tune`, `yardstick`,
`rsample`, `tidyr`, `dplyr`, `ggdist`, `hardhat`, `rsurv`.

## R CMD check results (spsurv 1.1.0)

```
Status: 0 errors, 0 test failures (local --as-cran, 2026-07-20)
Run full R CMD check --as-cran with vignettes before upload.
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
