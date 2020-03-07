setwd("~/Documents/spsurv")

if(!require(roxygen2)) install.packages('roxygen2')
if(!require(rstan)) install.packages('rstan')
if(!require(pkgbuild)) install.packages('pkgbuild')
if(!require(rstantools)) install.packages('rstantools')
if(!require(Rcpp)) install.packages('Rcpp')

try(remove.packages("spsurv", lib="~/R/x86_64-pc-linux-gnu-library/3.6"),
    silent = T)
try(roxygen2::roxygenize(clean = T), silent = T)

pkgbuild::compile_dll()
devtools::document()
devtools::install()
devtools::load_all(".")
devtools::build()
?spsurv::spbp
