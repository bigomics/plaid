# Normalize column medians of matrix

This function normalizes the column medians of matrix x. It calls
optimized functions from the matrixStats package.

## Usage

``` r
normalize_medians(x, ignore.zero = NULL)
```

## Arguments

- x:

  Input matrix

- ignore.zero:

  Logical indicating whether to ignore zeros to exclude for median
  calculation

## Value

Matrix with normalized column medians.

## Examples

``` r
# Create example matrix
set.seed(123)
x <- matrix(rnorm(100), nrow = 10, ncol = 10)
x[1:3, 1:3] <- 0  # Add some zeros

# Normalize medians
x_norm <- normalize_medians(x)
head(x_norm)
#>            [,1]        [,2]        [,3]        [,4]       [,5]       [,6]
#> [1,] 0.08971477  0.08971477  0.08971477  0.02598805 -0.2995912  0.1731354
#> [2,] 0.08971477  0.08971477  0.08971477 -0.69554765  0.1871985 -0.1087299
#> [3,] 0.08971477  0.08971477  0.08971477  0.49464949 -0.8702805 -0.1230536
#> [4,] 0.16022316  0.20039748 -0.63917646  0.47765732  2.5640718  1.2884191
#> [5,] 0.21900250 -0.46612637 -0.53532450  0.42110491  1.6030778 -0.3059541
#> [6,] 1.80477975  1.87662790 -1.59697854  0.28816408 -0.7279928  1.4362875
#>            [,7]        [,8]        [,9]       [,10]
#> [1,]  0.2910878 -0.01341431 -0.26305223  0.68965428
#> [2,] -0.5908751 -1.83155202  0.11646399  0.24454738
#> [3,] -0.4217591  1.48335538 -0.63947645 -0.06511785
#> [4,] -1.1071271 -0.23158391  0.37556013 -0.93175566
#> [5,] -1.1603429 -0.21039176 -0.48930298  1.05680287
#> [6,]  0.2149770  1.50318822  0.06296555 -0.90410917
```
