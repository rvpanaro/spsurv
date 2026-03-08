# Model.matrix method for fitted spbp models

Model.matrix of a fitted
[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md) model.

## Usage

``` r
# S3 method for class 'spbp'
model.matrix(object, ...)
```

## Arguments

- object:

  an object of class \`spbp\`, see
  [`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md).

- ...:

  arguments passed to parent method.

## Value

The model matrix.

## See also

[`spbp`](https://rvpanaro.github.io/spsurv/reference/spbp.md),
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html)

## Examples

``` r
library("spsurv")
data("veteran", package = "survival")
#> Warning: data set ‘veteran’ not found

fit <- bpph(Surv(time, status) ~ karno + factor(celltype),
  data = veteran
)
#> Priors are ignored because the MLE approach is used.

model.matrix(fit)
#>     karno factor(celltype)smallcell factor(celltype)adeno factor(celltype)large
#> 1      60                         0                     0                     0
#> 2      70                         0                     0                     0
#> 3      60                         0                     0                     0
#> 4      60                         0                     0                     0
#> 5      70                         0                     0                     0
#> 6      20                         0                     0                     0
#> 7      40                         0                     0                     0
#> 8      80                         0                     0                     0
#> 9      50                         0                     0                     0
#> 10     70                         0                     0                     0
#> 11     60                         0                     0                     0
#> 12     40                         0                     0                     0
#> 13     30                         0                     0                     0
#> 14     80                         0                     0                     0
#> 15     70                         0                     0                     0
#> 16     60                         1                     0                     0
#> 17     60                         1                     0                     0
#> 18     40                         1                     0                     0
#> 19     80                         1                     0                     0
#> 20     60                         1                     0                     0
#> 21     40                         1                     0                     0
#> 22     60                         1                     0                     0
#> 23     60                         1                     0                     0
#> 24     30                         1                     0                     0
#> 25     80                         1                     0                     0
#> 26     30                         1                     0                     0
#> 27     50                         1                     0                     0
#> 28     60                         1                     0                     0
#> 29     80                         1                     0                     0
#> 30     40                         1                     0                     0
#> 31     20                         1                     0                     0
#> 32     80                         1                     0                     0
#> 33     30                         1                     0                     0
#> 34     75                         1                     0                     0
#> 35     70                         1                     0                     0
#> 36     60                         1                     0                     0
#> 37     30                         1                     0                     0
#> 38     60                         1                     0                     0
#> 39     80                         1                     0                     0
#> 40     60                         1                     0                     0
#> 41     70                         1                     0                     0
#> 42     50                         1                     0                     0
#> 43     50                         1                     0                     0
#> 44     40                         1                     0                     0
#> 45     40                         1                     0                     0
#> 46     20                         0                     1                     0
#> 47     70                         0                     1                     0
#> 48     40                         0                     1                     0
#> 49     80                         0                     1                     0
#> 50     80                         0                     1                     0
#> 51     50                         0                     1                     0
#> 52     80                         0                     1                     0
#> 53     30                         0                     1                     0
#> 54     80                         0                     1                     0
#> 55     50                         0                     0                     1
#> 56     80                         0                     0                     1
#> 57     50                         0                     0                     1
#> 58     70                         0                     0                     1
#> 59     60                         0                     0                     1
#> 60     40                         0                     0                     1
#> 61     80                         0                     0                     1
#> 62     80                         0                     0                     1
#> 63     70                         0                     0                     1
#> 64     90                         0                     0                     1
#> 65     90                         0                     0                     1
#> 66     80                         0                     0                     1
#> 67     80                         0                     0                     1
#> 68     70                         0                     0                     1
#> 69     60                         0                     0                     1
#> 70     90                         0                     0                     0
#> 71     80                         0                     0                     0
#> 72     80                         0                     0                     0
#> 73     50                         0                     0                     0
#> 74     50                         0                     0                     0
#> 75     70                         0                     0                     0
#> 76     70                         0                     0                     0
#> 77     20                         0                     0                     0
#> 78     60                         0                     0                     0
#> 79     90                         0                     0                     0
#> 80     30                         0                     0                     0
#> 81     20                         0                     0                     0
#> 82     70                         0                     0                     0
#> 83     90                         0                     0                     0
#> 84     80                         0                     0                     0
#> 85     50                         0                     0                     0
#> 86     70                         0                     0                     0
#> 87     60                         0                     0                     0
#> 88     90                         0                     0                     0
#> 89     50                         0                     0                     0
#> 90     30                         1                     0                     0
#> 91     70                         1                     0                     0
#> 92     20                         1                     0                     0
#> 93     30                         1                     0                     0
#> 94     60                         1                     0                     0
#> 95     40                         1                     0                     0
#> 96     30                         1                     0                     0
#> 97     20                         1                     0                     0
#> 98     60                         1                     0                     0
#> 99     70                         1                     0                     0
#> 100    80                         1                     0                     0
#> 101    85                         1                     0                     0
#> 102    70                         1                     0                     0
#> 103    70                         1                     0                     0
#> 104    70                         1                     0                     0
#> 105    50                         1                     0                     0
#> 106    30                         1                     0                     0
#> 107    40                         1                     0                     0
#> 108    40                         0                     1                     0
#> 109    40                         0                     1                     0
#> 110    99                         0                     1                     0
#> 111    80                         0                     1                     0
#> 112    60                         0                     1                     0
#> 113    60                         0                     1                     0
#> 114    60                         0                     1                     0
#> 115    60                         0                     1                     0
#> 116    50                         0                     1                     0
#> 117    70                         0                     1                     0
#> 118    10                         0                     1                     0
#> 119    40                         0                     1                     0
#> 120    70                         0                     1                     0
#> 121    90                         0                     1                     0
#> 122    80                         0                     1                     0
#> 123    50                         0                     1                     0
#> 124    40                         0                     1                     0
#> 125    40                         0                     1                     0
#> 126    60                         0                     0                     1
#> 127    70                         0                     0                     1
#> 128    30                         0                     0                     1
#> 129    60                         0                     0                     1
#> 130    30                         0                     0                     1
#> 131    60                         0                     0                     1
#> 132    80                         0                     0                     1
#> 133    75                         0                     0                     1
#> 134    60                         0                     0                     1
#> 135    70                         0                     0                     1
#> 136    80                         0                     0                     1
#> 137    30                         0                     0                     1
```
