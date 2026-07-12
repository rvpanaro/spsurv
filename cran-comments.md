## Test environments

* local: R 4.5.x (macOS, aarch64), `R CMD check --as-cran`
* rhub: run `rhub::check_for_cran()` before upload
* win-builder: run before upload

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
