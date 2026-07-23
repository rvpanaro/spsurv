# Bernstein basis polynomials calculations

Bernstein basis polynomials calculations

## Usage

``` r
bp.basis(time, degree, tau = max(time))
```

## Arguments

- time:

  a vector of times.

- degree:

  Bernstein polynomial degree

- tau:

  must be greater than times maximum value observed.

## Value

A list containing matrices g and G corresponding BP basis and
corresponding tau value used to compute them.
