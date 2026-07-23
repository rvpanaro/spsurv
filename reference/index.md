# Package index

## Details

- [`spsurv-package`](https://rvpanaro.github.io/spsurv/reference/spsurv-package.md)
  [`spsurv`](https://rvpanaro.github.io/spsurv/reference/spsurv-package.md)
  : The 'spsurv' package.

## Fitter functions

Functions for fitting a BP based survival regression model.

- [`bpph()`](https://rvpanaro.github.io/spsurv/reference/bpph.md) :
  Bernstein PH Model
- [`bppo()`](https://rvpanaro.github.io/spsurv/reference/bppo.md) :
  Bernstein PO Model
- [`bpaft()`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) :
  Bernstein AFT Model
- [`spbp()`](https://rvpanaro.github.io/spsurv/reference/spbp.md) :
  spbp: The BP Based Survival Analysis Function
- [`spbp(`*`<default>`*`)`](https://rvpanaro.github.io/spsurv/reference/spbp.default.md)
  : spbp: The BP Based Semiparametric Survival Analysis Function

## Extraction and inference

Coefficients, confidence and credible intervals, residuals, survival
curves.

- [`coef(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/coef.spbp.md)
  : Estimated regression coefficients
- [`confint(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/confint.spbp.md)
  : Confidence intervals for the regression coefficients
- [`credint()`](https://rvpanaro.github.io/spsurv/reference/credint.md)
  : Generic S3 method credint
- [`credint(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/credint.spbp.md)
  : Credible intervals for the regression coefficients
- [`vcov(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/vcov.spbp.md)
  : Covariance of the regression coefficients
- [`residuals(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/residuals.spbp.md)
  : BP based models residuals.
- [`survfit(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)
  [`as.data.frame(`*`<survfitbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)
  : BP-based model survival curves
- [`predict(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/predict.spbp.md)
  : Predicted survival as a tidy data frame
- [`augment(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/augment.spbp.md)
  : Augment data with spbp model information
- [`tidy(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/tidy.spbp.md)
  : Tidy summary for fitted spbp models
- [`glance(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/glance.spbp.md)
  : Glance at a fitted spbp model
- [`model.matrix(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/model.matrix.spbp.md)
  : Model.matrix method for fitted spbp models
- [`as_draws_df.spbp()`](https://rvpanaro.github.io/spsurv/reference/as_draws_df.spbp.md)
  : Convert spbp Bayesian fit to posterior draws
- [`spread_surv_draws.spbp()`](https://rvpanaro.github.io/spsurv/reference/spread_surv_draws.spbp.md)
  : Posterior survival curves in long tidy format

## tidymodels

parsnip model specifications and engines.

- [`proportional_odds()`](https://rvpanaro.github.io/spsurv/reference/proportional_odds.md)
  : Proportional odds regression model specification
- [`bp_survival_reg()`](https://rvpanaro.github.io/spsurv/reference/bp_survival_reg.md)
  : Bernstein survival regression specification

## Model comparison and extraction

Bernstein baseline helper, parameter estimates, and information
criteria.

- [`bernstein()`](https://rvpanaro.github.io/spsurv/reference/bernstein.md)
  : Bernstein polynomial baseline specification
- [`estimates()`](https://rvpanaro.github.io/spsurv/reference/estimates.md)
  : Parameter estimates for fitted spbp models
- [`se()`](https://rvpanaro.github.io/spsurv/reference/se.md) : Standard
  errors for fitted spbp model parameters
- [`logLik(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/logLik.spbp.md)
  : Log-likelihood for fitted spbp models
- [`AIC(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/AIC.spbp.md)
  : Akaike information criterion for fitted spbp models
- [`extractAIC(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/extractAIC.spbp.md)
  : Extract AIC from a fitted spbp model
- [`anova(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/anova.spbp.md)
  : Analysis of deviance table for nested spbp models
- [`rank_models()`](https://rvpanaro.github.io/spsurv/reference/rank_models.md)
  : Rank fitted spbp models by AIC

## Bernstein polynomial basis

Functions for the BP basis.

- [`bp.basis()`](https://rvpanaro.github.io/spsurv/reference/bp.basis.md)
  : Bernstein basis polynomials calculations
- [`pw.basis()`](https://rvpanaro.github.io/spsurv/reference/pw.basis.md)
  : Power basis polynomials calculations

## Plot and print

S3 methods for spbp objects.

- [`plot(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/plot.spbp.md)
  : BP based models plot.
- [`ggresiduals()`](https://rvpanaro.github.io/spsurv/reference/ggresiduals.md)
  : ggplot2 residual diagnostic plots for spbp models
- [`print(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.spbp.md)
  : Bernstein Polynomial Based Regression Object Print
- [`summary(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/summary.spbp.md)
  : Bernstein Polynomial Based Regression Object Summary
- [`print(`*`<summary.bpaft.bayes>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bpaft.bayes.md)
  : Bernstein Polynomial Based Regression Object Summary BPAFT Bayes
- [`print(`*`<summary.bpaft.mle>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bpaft.mle.md)
  : Bernstein Polynomial Based Regression Object Summary BPAFT MLE
- [`print(`*`<summary.bpph.bayes>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bpph.bayes.md)
  : Bernstein Polynomial Based Regression Object Summary BPPH Bayes
- [`print(`*`<summary.bpph.mle>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bpph.mle.md)
  : Bernstein Polynomial Based Regression Object Summary BPPH MLE
- [`print(`*`<summary.bppo.bayes>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bppo.bayes.md)
  : Bernstein Polynomial Based Regression Object Summary BPPO Bayes
- [`print(`*`<summary.bppo.mle>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.bppo.mle.md)
  : Bernstein Polynomial Based Regression Object BPPO MLE
- [`print(`*`<summary.spbp.bayes>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.spbp.bayes.md)
  : Bernstein Polynomial Based Regression Object Summary Bayes
- [`print(`*`<summary.spbp.mle>`*`)`](https://rvpanaro.github.io/spsurv/reference/print.summary.spbp.mle.md)
  : Bernstein Polynomial Based Regression Object Summary MLE
