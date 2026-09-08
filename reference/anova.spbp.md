# Analysis of deviance table for nested spbp models

Compare nested MLE fits. With one model, builds a sequential table by
refitting intercept-only and nested submodels along the formula term
order (same idea as `anova(fit)` for
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)). With
multiple models, compares each consecutive pair (same convention as
`anova(fit0, fit1)` for `survreg`).

## Usage

``` r
# S3 method for class 'spbp'
anova(object, ..., test = "Chisq")
```

## Arguments

- object:

  A fitted `"spbp"` object.

- ...:

  Additional nested fitted models (MLE, same family recommended).

- test:

  Which test to report (only `"Chisq"` is supported).

## Value

An `"anova"` object (also a `data.frame`). Pairwise comparisons use
`Terms`, `Resid. Df`, `-2*LL`, `Test`, `Df`, `Deviance`, and `Pr(>Chi)`
columns (as in
[`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)).
Sequential single-model tables use `Df`, `Deviance`, `Resid. Df`,
`-2*LL`, and `Pr(>Chi)` with row names `NULL` and term labels.
