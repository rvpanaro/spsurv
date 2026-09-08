# The 'spsurv' package.

A set of flexible routines to allow semiparametric survival regression
modeling based on Bernstein polynomial, including Bernstein PH model
(BPPH), Bernstein PO model (BPPO), and Bernstein AFT model (BPAFT) for
right-censored data.

## Value

none

## Details

[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md) fits
semi-parametric models for time-to-event survival data. Non-informative
right-censoring assumption is available. Any user-defined Bernstein
polynomial can be user-defined using an arbitrary degree, i.e. highest
basis polynomials order.

The framework takes advantage of fully likelihood methods since the
polynomial parameters are used to estimate the baseline functions. Even
so, this is said to be semi-parametric since this approach does not rely
on any distribution. Unlike the Cox model, the BP based models provide
smooth hazard and survival curve estimates. For Bayesian fits, users
should routinely inspect divergences, split-\\\hat R\\, and effective
sample sizes, and tighten NUTS controls (e.g., increasing `adapt_delta`,
`iter`, and `warmup`) when warnings are present.

## References

Panaro R.V. (2020). spsurv: An R package for semi-parametric survival
analysis. arXiv preprint arXiv:2003.10548.

Demarqui, F. N., & Mayrink, V. D. (2019). A fully likelihood-based
approach to model survival data with crossing survival curves. arXiv
preprint arXiv:1910.02406.

Demarqui, F. N., Mayrink, V. D., & Ghosh, S. K. (2019). An Unified
Semiparametric Approach to Model Lifetime Data with Crossing Survival
Curves. arXiv preprint arXiv:1910.04475.

Osman, M., & Ghosh, S. K. (2012). Nonparametric regression models for
right-censored data using Bernstein polynomials. Computational
Statistics & Data Analysis, 56(3), 559-573.

Lorentz, G. G. (1953). Bernstein polynomials. American Mathematical
Society.

## Author

rvpanaro@gmail.com
