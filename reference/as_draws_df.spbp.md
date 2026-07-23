# Convert spbp Bayesian fit to posterior draws

Convert spbp Bayesian fit to posterior draws

## Usage

``` r
as_draws_df.spbp(x, variables = c("beta", "gamma"), ...)
```

## Arguments

- x:

  A fitted `"spbp"` object from `approach = "bayes"`.

- variables:

  Character vector of posterior components to include (`"beta"`,
  `"gamma"`, `"log_lik"`).

- ...:

  Not used.

## Value

A
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
object on the back-transformed scale.
