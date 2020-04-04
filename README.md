## spsurv <a href='https://rvpanaro.github.io/spsurv'><img src='https://raw.githubusercontent.com/rvpanaro/spsurv/master/inst/figures/logo.png' align="right" height="139" /></a>

---

[![CRAN status](https://www.r-pkg.org/badges/version/spsurv)](https://cran.r-project.org/package=spsurv)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/spsurv)](http://www.r-pkg.org/pkg/spsurv)
[![R build status](https://github.com/rvpanaro/spsurv/workflows/R-CMD-check/badge.svg)](https://github.com/rvpanaro/spsurv/actions)
[![Travis build status](https://travis-ci.org/rvpanaro/spsurv.svg?branch=master)](https://travis-ci.org/rvpanaro/spsurv)
[![Codecov test coverage](https://codecov.io/gh/rvpanaro/spsurv/branch/master/graph/badge.svg)](https://codecov.io/gh/rvpanaro/spsurv?branch=master)

## Overview
<div style="text-align: justify">
spsurv: An R package for semi-parametric survival analysis.

The **spsurv** package was designed to contribute with a flexible set of semi-parametric survival regression options, including proportional hazards (PH), proportional odds (PO), and accelerated failure time (AFT) models for right-censored data.

The proposed package provides:

- Survival classes (PH, PO, AFT) extensions based on a fully likelihood-based approach for either Bayesian or maximum likelihood (ML) estimation procedures,

- smooth estimates for the unknown baseline functions based on the Bernstein polynomial (BP),

-  integration with Stan software aiming more flexibility in terms of user-defined modeling,

- six distinct prior specification options in a Bayesian analysis.

[Stan](https://mc-stan.org) is an open-source platform that has a specific language and many built-in log-probability functions that can be used to define custom likelihood functions and prior
specifications. The program has extensive supporting literature
available online such as reference manuals, forums, articles, and books for users and developers. The Stan currently defaults to No-U-Turn (NUTS) sampling, which consists of a Hamiltonian Monte Carlo (HMC) extension that explores the posterior distribution more efficiently, access to Stan can be established through several modules integration such as [rstan](https://mc-stan.org/users/interfaces/rstan). </div>

## Installation

- Install the CRAN version:
```r
install.packages("spsurv")
```

- Installation using the devtools package:

```r
install.packages("devtools")
devtools::install_github("rvpanaro/spsurv")
```

## Usage 

- Fit a BP based survival regression PH model using:

```r
library("KMsurv")
data("larynx")

library(spsurv)

fit <- bpph(Surv(time, delta) ~ age + factor(stage), model = "ph", data = larynx)
summary(fit)      
```

- Alternatively, one can use the `spbp` function:

```r
fit <- spbp(Surv(time, delta) ~ age + factor(stage), model = "ph",  data = larynx)
summary(fit)      
```
- Access to Bayesian analysis using `approach` argument:

```
## NUTS sampling (Bayesian)
fit2 <- spbp(Surv(time, delta) ~ age + factor(stage),
                     approach = "bayes",  data = larynx,
                     iter = 2000, chains = 1, warmup = 1000)
summary(fit2)
```

See the [reference manual](https://rvpanaro.github.io/spsurv/reference/index.html) for more examples.

### Upcoming features (To do)

- \[ \]  `survivor` class.
- \[ \]  `survivor` ggplot method.
- \[ \]  Deviance, martingale and standard residuals.
- \[ \]  Frailty models coding coverage.

## Troubleshooting

Please report issues at https://github.com/rvpanaro/spsurv/issues or mail it to [renatovp@ime.usp.br](mailto:lunde@adobe.com?subject=[GitHub]%20spsurv%20issue%20found), we will answer as soon as possible.

