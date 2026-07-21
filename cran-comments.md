## Resubmission context

This is a resubmission of **spsurv 1.1.0** after the package was archived on
CRAN because the former maintainer email address was unreachable.

* **Maintainer:** Renato Panaro (`rvpanaro@gmail.com`) is the package maintainer
  and the only person with the `cre` role in `Authors@R`.
* **Version:** increased from 1.0.0 (archived CRAN release) to 1.1.0 for this
  resubmission.
* **Archived-package issues addressed:** updated maintainer contact; corrected
  package citation in `Description` (`<doi:10.48550/arXiv.2003.10548>`);
  reviewed `Authors@R` (single maintainer; thesis advisors as `ths`/`rev` only);
  restored `rstantools` in `Imports` with explicit `NAMESPACE` import for
  `configure`/`configure.win`.

## Test environments

* local: R 4.5.x (macOS, aarch64), `R CMD check --as-cran`
* win-builder: run before upload
* macOS builder: run before upload if available

## R CMD check results (spsurv 1.1.0, local macOS aarch64, R 4.5.1, 2026-07-21)

```
R CMD build .  -> spsurv_1.1.0.tar.gz (OK)
R CMD check --as-cran spsurv_1.1.0.tar.gz

Status: 3 NOTEs
0 errors | 0 warnings | 0 failures
```

* **NOTE:** `checking CRAN incoming feasibility` — new resubmission; package was
  archived; maintainer `rvpanaro@gmail.com` (no invalid URL findings after
  README/`DESCRIPTION` URL canonicalisation).
* **NOTE:** `checking for future file timestamps` — environmental (sandbox/time).
* **NOTE:** `checking HTML version of manual` — local `tidy` too old; not a
  package issue.
* **INFO:** `GNU make is a SystemRequirements` — expected for Stan.
* **Resolved before final check:** unstated test dependencies (`covr`, `rsample`,
  `rsurv`) — added to `Suggests` for optional/skipped tests.

## testthat (local, spsurv 1.1.0, 2026-07-21)

```
devtools::test(): 665 passed, 0 failed, 1 skipped, 43 warnings
(Stan / ill-conditioned gamma diagnostics)
```

## revdepcheck

Not run for this resubmission.

## Submission notes

* Minor release 1.1.0: tidymodels integration (parsnip engines `spsurv` and
  `spsurv_bayes`, `proportional_odds()`, `survival_reg()`, `bp_survival_reg()`),
  tidybayes/posterior draw helpers, censored-style `predict()` types
  (`survival`, `time`, `linear_pred`), and `augment.spbp()`.
* New vignette `vignette("tidymodels")`.
* Installed size is dominated by compiled Stan code in `libs/` (expected).
* Bayesian tests use reduced `iter`/`chains` for speed.
