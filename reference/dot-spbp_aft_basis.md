# Internal: AFT basis evaluator aligned with Stan transformed parameters

Internal: AFT basis evaluator aligned with Stan transformed parameters

## Usage

``` r
.spbp_aft_basis(time_aft, tau_a, tau_b, P)
```

## Arguments

- time_aft:

  Numeric vector of AFT residual-scale times y = log(t) - eta.

- tau_a:

  Lower bound used by the fitted model.

- tau_b:

  Upper bound used by the fitted model.

- P:

  Power-basis transformation matrix from
  [`pw.basis()`](https://rvpanaro.github.io/spsurv/reference/pw.basis.md).

## Value

`list(g = ..., G = ...)` where `g` is hazard basis and `G` is
cumulative-hazard basis on the AFT residual scale.
