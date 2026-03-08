# spbp: The BP Based Semiparametric Survival Analysis Function

spbp: The BP Based Semiparametric Survival Analysis Function

## Usage

``` r
# Default S3 method
spbp(
  formula,
  degree,
  data,
  approach = c("mle", "bayes"),
  model = c("ph", "po", "aft"),
  priors = list(beta = c("normal(0,4)"), gamma = c("lognormal(0,4)"), frailty =
    c("gamma(0.01,0.01)")),
  cores = min(parallel::detectCores() - 1, 4),
  scale = TRUE,
  verbose = FALSE,
  chains = 4,
  ...
)
```

## Arguments

- formula:

  a Surv object with time to event, status and explanatory terms

- degree:

  Bernstein Polynomial degree

- data:

  a data.frame object

- approach:

  Bayesian or Maximum Likelihood estimation methods, default is approach
  = "bayes"

- model:

  Proportional Hazards or Proportional Odds BP based regression, default
  is model = "ph"

- priors:

  prior settings for the Bayesian approach; \`normal\` or \`cauchy\` for
  beta; \`lognormal\` or \`loglogistic\` for gamma (BP coefficients)

- cores:

  number of core threads to use

- scale:

  logical; indicates whether to center and scale the data

- verbose:

  verbose passed to stan

- chains:

  number of chains passed to stan

- ...:

  further arguments passed to or from other methods

## Value

An object of class `spbp`
