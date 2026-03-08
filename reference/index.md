# Package index

## Details

- [`spsurv-package`](https://rvpanaro.github.io/spsurv/reference/spsurv-package.md)
  [`spsurv`](https://rvpanaro.github.io/spsurv/reference/spsurv-package.md)
  : The 'spsurv' package.

## Fitter functions

Functions for fitting a BP based survival regression model.

- [`bpph()`](https://rvpanaro.github.io/spsurv/reference/bpph.md) :
  Bernstein Polynomial Based Proportional Hazards Model
- [`bppo()`](https://rvpanaro.github.io/spsurv/reference/bppo.md) :
  Bernstein Polynomial Based Proportional Odds Model
- [`bpaft()`](https://rvpanaro.github.io/spsurv/reference/bpaft.md) :
  Bernstein Polynomial Based Accelerated Failure Time Model
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
  : Confidence intervals for the regression coefficients
- [`vcov(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/vcov.spbp.md)
  : Covariance of the regression coefficients
- [`residuals(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/residuals.spbp.md)
  : BP based models residuals.
- [`survfit(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/survfit.spbp.md)
  : BP-based model survival curves
- [`model.matrix(`*`<spbp>`*`)`](https://rvpanaro.github.io/spsurv/reference/model.matrix.spbp.md)
  : Model.matrix method for fitted spbp models

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
