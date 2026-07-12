library(testthat)
library(parsnip)
library(generics)
library(spsurv)
data("veteran", package = "survival")

test_check("spsurv")
