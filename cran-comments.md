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
  disabled install-time Stan auto-config (`Config/rstantools/auto_config: FALSE`)
  so Windows builds use committed Stan C++ exports and do not call
  `rstantools::rstan_config()` during `configure`; added `RcppParallel` and
  aligned `src/Makevars.win` with Unix TBB flags (fixes missing
  `tbb/tbb_stddef.h`). R-devel Win-builder can still fail later when loading
  Import `rstan.dll` (guest environment); R-release Windows check is clean.

## Test environments

* local: R 4.6.1 (macOS Sequoia, Apple Silicon), `R CMD check --as-cran`
* win-builder R-release: R 4.6.1 (ucrt), x86_64-w64-mingw32
  (https://win-builder.r-project.org/6bFfd06dsW09) — Status: 1 NOTE
* win-builder R-devel (R 4.7.x guest): compile succeeds after `RcppParallel` /
  `Makevars.win` sync, but install fails at lazy-load when loading dependency
  `rstan.dll` (`LoadLibrary` / “Das angegebene Modul wurde nicht gefunden”).
  This is an R-devel Win-builder `rstan` binary/environment issue, not an
  `spsurv` compile or NAMESPACE error; R-release Windows check is clean.
* macOS builder: not yet run

## R CMD check results (spsurv 1.1.0, local macOS aarch64, R 4.6.1, 2026-07-21)

```
R CMD build .  -> spsurv_1.1.0.tar.gz (OK)
R CMD check --as-cran spsurv_1.1.0.tar.gz

Status: 2 NOTEs
0 errors | 0 warnings | 0 failures
```

Installation, compilation, examples, tests, vignettes, and PDF manual: OK.

* **NOTE:** `checking CRAN incoming feasibility` — archived-package resubmission;
  updated maintainer address `rvpanaro@gmail.com`; `GPL-3 + file LICENSE`; no
  invalid URL findings.
* **NOTE:** `checking HTML version of manual` — local HTML Tidy was not recent
  enough to validate the manual HTML; the optional V8-based math-rendering check
  was unavailable locally.
* **INFO:** `GNU make is a SystemRequirements` — expected for Stan.

## R CMD check results (win-builder R-release, R 4.6.1 ucrt, 2026-07-23)

```
https://win-builder.r-project.org/6bFfd06dsW09
Installation time: ~164s | Check time: ~321s

Status: 1 NOTE
0 errors | 0 warnings | 0 failures
```

Install, examples, tests, vignettes, PDF/HTML manuals: OK. Binary:
`spsurv_1.1.0.zip`.

* **NOTE:** `checking CRAN incoming feasibility` — New submission; package was
  archived on CRAN; `GPL-3 + file LICENSE`; maintainer
  `Renato Panaro <rvpanaro@gmail.com>`; CRAN comment notes archival on
  2026-06-11 due to undeliverable maintainer email (addressed by this
  resubmission).
* **INFO:** `checking C++ specification` — specified C++17 (expected for Stan).

## R CMD check results (win-builder R-devel, R 4.7.x guest)

Compile of `spsurv` succeeds. Install then fails at byte-compile / lazy-load:

```
unable to load shared object '.../rstan/libs/x64/rstan.dll':
  LoadLibrary failure: Das angegebene Modul wurde nicht gefunden.
ERROR: lazy loading failed for package 'spsurv'
```

Cause: loading Import `rstan` on that guest fails inside `rstan.dll` (missing
dependent module). Same class of failure previously seen at `configure` time;
with `auto_config: FALSE` it now appears only when R loads Imports. Not
actionable in `spsurv` without dropping `rstan` from Imports (incompatible with
Stan MLE/Bayes engines). R-release Win-builder check above is the Windows
reference result for this submission.

## testthat (local, spsurv 1.1.0, R 4.6.1)

```
devtools::test(): 678 passed, 0 failed, 1 skipped
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
* Stan C++ exports are committed; regenerate with `rstantools::rstan_config()`
  only when `inst/stan` models change (developers: install Suggests
  `rstantools`).
