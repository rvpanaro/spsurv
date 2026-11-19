# Bernstein polynomial baseline specification

Helper for specifying a Bernstein polynomial baseline via the `dist` /
`baseline` arguments accepted by
[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md).

## Usage

``` r
bernstein(m = NULL)
```

## Arguments

- m:

  Bernstein polynomial degree (number of basis coefficients).

## Value

A list with components `baseline` and `m`.

## Examples

``` r
bernstein(5)
#> $baseline
#> [1] "bernstein"
#> 
#> $m
#> [1] 5
#> 
```
