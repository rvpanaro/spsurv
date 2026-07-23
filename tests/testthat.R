library(testthat)
library(parsnip)
library(generics)
library(spsurv)
veteran <- survival::veteran

test_check("spsurv")
