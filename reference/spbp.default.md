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
  cores = .spbp_default_cores(),
  scale = TRUE,
  dist = NULL,
  baseline = NULL,
  verbose = FALSE,
  chains = 4,
  ...
)
```

## Arguments

- formula:

  a Surv object with time to event, status and explanatory terms

- degree:

  Bernstein polynomial degree (integer). If omitted and neither `dist`
  nor `baseline` supplies
  [`bernstein`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)`(m)`,
  the default is `ceiling(sqrt(n))` where `n` is the number of rows in
  `data`.

- data:

  a data.frame object

- approach:

  Bayesian or Maximum Likelihood estimation methods; default is `"mle"`

- model:

  Bernstein PH (`"ph"`), PO (`"po"`), or AFT (`"aft"`) model; default is
  `"ph"`

- priors:

  prior settings for the Bayesian approach; \`normal\` or \`cauchy\` for
  beta; \`lognormal\` or \`loglogistic\` for gamma (BP coefficients)

- cores:

  number of core threads to use (Bayes sampling)

- scale:

  logical; indicates whether to center and scale the data

- dist:

  optional baseline specification; use
  [`bernstein`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)`(m)`
  for the Bernstein polynomial degree

- baseline:

  optional alias for `dist`

- verbose:

  passed to Stan

- chains:

  number of MCMC chains (Bayes)

- ...:

  further arguments passed to
  [`rstan::optimizing`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)
  (MLE) or
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  (Bayes), e.g. `iter`, `warmup`, `init`.

## Value

An object of class `spbp`. Component `degree` records the Bernstein
polynomial degree used in the fit (also stored in `call$degree`).

## Details

Right-censored survival data are modeled with a Bernstein-polynomial
baseline and regression on covariates. With `approach = "mle"`,
parameters are estimated by Stan's optimizer and approximate inference
uses the Hessian when available. With `approach = "bayes"`, posterior
samples are drawn with NUTS; use `summary`,
[`tidy.spbp`](https://rvpanaro.github.io/spsurv/reference/tidy.spbp.md),
and
[`glance.spbp`](https://rvpanaro.github.io/spsurv/reference/glance.spbp.md)
for output. Covariates are centered and scaled when `scale = TRUE`
(default).

The returned object includes:

- `coefficients`:

  Regression estimates on the original covariate scale.

- `bp.param`:

  Bernstein baseline coefficients (`gamma`).

- `degree`:

  Polynomial degree used (also in `call$degree`).

- `loglik`:

  MLE: intercept-only and full-model log-likelihoods; Bayes: posterior
  mean pointwise log-likelihoods.

- `call`:

  Matched call with `approach`, `model`, and `degree`.

## See also

[`bpph`](https://rvpanaro.github.io/spsurv/reference/bpph.md),
[`bppo`](https://rvpanaro.github.io/spsurv/reference/bppo.md),
[`bpaft`](https://rvpanaro.github.io/spsurv/reference/bpaft.md),
[`bernstein`](https://rvpanaro.github.io/spsurv/reference/bernstein.md),
[`summary.spbp`](https://rvpanaro.github.io/spsurv/reference/summary.spbp.md)
