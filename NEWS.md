# spsurv 1.0.2

* **Handlers:** Removed `handlers.R`; validation and setup logic inlined in `spbp.default`.
* **vcov:** Block-wise inversion of the Hessian (regression vs Bernstein polynomial blocks) in `vcov.spbp`.
* **Stan (AFT):** No clamping of `u` to avoid boundary bias; minimum feasible divisor (`min_range`) only when range is degenerate.
* **MLE:** Initial values for the optimizer drawn from prior distributions; retries limited to 3 attempts.
* **Examples:** `data("veteran", package = "survival")` in docs and examples.
* **pkgdown:** Reference index updated to match current package; Bootstrap 5 / default template; logo size and alt-text; README aligned with site style; itsamp removed from index.
* **Git:** `.gitignore` updated (rj/, .Rcheck, etc.); GitHub Action for pkgdown deploy to gh-pages.

# spsurv 1.0.1

* [pkgdown](https://github.com/r-lib/pkgdown) reference manual [https://rvpanaro.github.io/spsurv/reference/index.html](https://rvpanaro.github.io/spsurv/reference/index.html).

# spsurv 1.0.0

* This is the first release of spsurv.
