# spsurv <a href='https://rvpanaro.github.io/spsurv'><img src='inst/figures/logo.png' align="right" height="139" /></a>
---

[![CRAN status](https://www.r-pkg.org/badges/version/spsurv)](https://cran.r-project.org/package=spsurv)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/spsurv)](http://www.r-pkg.org/pkg/spsurv)
[![R build status](https://github.com/rvpanaro/spsurv/workflows/R-CMD-check/badge.svg)](https://github.com/rvpanaro/spsurv/actions)
[![Travis build status](https://travis-ci.org/rvpanaro/spsurv.svg?branch=master)](https://travis-ci.org/rvpanaro/spsurv)
[![Codecov test coverage](https://codecov.io/gh/rvpanaro/spsurv/branch/master/graph/badge.svg)](https://codecov.io/gh/rvpanaro/spsurv?branch=master)

## Overview

An R package for semi-parametric survival analysis.

The *spsurv* package was designed to contribute with a flexible set of semi-parametric survival regression modelings, including proportional hazards (PH), proportional odds (PO), and accelerated failure time (AFT) models for right-censored data.

## Installation

- Install CRAN version of the package.
```r
install.packages("spsurv")
```

- Install and load the *spsurv* package using the devtools package.

```r
install.packages("devtools")

library(devtools)
install_github("rvpanaro/spsurv")
```

## Usage 

- Check out the main fitter function examples.

```r
library("KMsurv")
data("larynx")

library(spsurv)

## Maximum Likelihood
fit <- spbp(Surv(time, delta)~age+factor(stage),
                    approach = "mle",  data = larynx)
summary(fit)      

## NUTS sampling (Bayesian)
fit2 <- spbp(Surv(time, delta)~age+factor(stage),
                     approach = "bayes",  data = larynx,
                     iter = 2000, chains = 1, warmup = 1000)
summary(fit2)
```

The *spsurv* already provides:
- Integration with Stan software.
- Estimates either in Bayesian or Frequentist (point estimate) inferential approaches.
- Three survival regression classes: PH, PO and AFT.
- Six distinct prior specifications in a Bayesian analysis.
