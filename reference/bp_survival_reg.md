# Bernstein survival regression specification

Convenience constructor mapping `family` to the appropriate parsnip
model (`proportional_hazards`, `proportional_odds`, or `survival_reg`).

## Usage

``` r
bp_survival_reg(
  family = c("ph", "po", "aft"),
  mode = "censored regression",
  engine = "spsurv"
)
```

## Arguments

- family:

  Model family: `"ph"`, `"po"`, or `"aft"`.

- mode:

  Model mode (default `"censored regression"`).

- engine:

  parsnip engine (default `"spsurv"` for MLE).

## Value

A parsnip model specification.
