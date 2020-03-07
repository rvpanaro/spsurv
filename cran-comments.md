## Test environments

* using R version 3.5.2 (2018-12-20)
* using platform: x86_64-pc-linux-gnu (64-bit)
* using session charset: UTF-8

## R CMD check results

There were no ERRORs. 

There were 2 WARNINGs:

* checking S3 generic/method consistency ... WARNING
  survivor:
    function(spbp, ...)
  survivor.default:
    function(time, arg, newdata, model, approach, ...)

* checking compilation flags used ... WARNING
  Compilation used the following non-portable flag(s):
    ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’
    
There were 3 NOTEs:

* checking installed package size ... NOTE
    installed size is 77.8Mb
    sub-directories of 1Mb or more:
      libs  77.4Mb
      
* checking R code for possible problems ... NOTE
    no visible binding for global variables

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

0 errors ✓ | 2 warnings x | 3 notes x
* checking installed package size ... NOTE
    installed size is 65.0Mb
    sub-directories of 1Mb or more:
      libs  64.7Mb

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
  
  
