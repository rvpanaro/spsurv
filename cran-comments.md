## Test environments

* local: R 4.5.x (macOS, aarch64), `R CMD check --as-cran`
* rhub: run `rhub::check_for_cran()` before upload
* win-builder: run before upload

## R CMD check results (spsurv 1.0.2)

```
Status: 2 NOTEs, 0 warnings, 0 errors
```

Local check on 2026-06-29:

* **INFO** GNU make is a SystemRequirements (expected; Stan build).
* **NOTE** HTML manual validation skipped (no recent `tidy` on this machine).

All other `--as-cran` checks passed, including examples and tests.

## revdepcheck

No reverse dependencies known at time of release.

## Submission notes

* Resubmission of 1.0.2 with tidy/glance/AIC helpers, `bernstein()` baseline
  syntax, and improved `survfit`/`residuals`.
* Installed size is dominated by compiled Stan code in `libs/` (expected).
* Bayesian tests use reduced `iter`/`chains` for speed.
