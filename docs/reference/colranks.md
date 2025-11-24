# Compute columnwise ranks of matrix

Computes columnwise rank of matrix. Can be sparse. Tries to call
optimized functions from Rfast or matrixStats.

## Usage

``` r
colranks(
  X,
  sparse = NULL,
  signed = FALSE,
  keep.zero = FALSE,
  ties.method = "average"
)
```

## Arguments

- X:

  Input matrix

- sparse:

  Logical indicating to use sparse methods

- signed:

  Logical indicating using signed ranks

- keep.zero:

  Logical indicating whether to keep zero as ranked zero

- ties.method:

  Character Choice of ties.method

## Value

Matrix of columnwise ranks with same dimensions as input.

## Examples

``` r
# Create example matrix
set.seed(123)
X <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(X) <- paste0("Gene", 1:10)
colnames(X) <- paste0("Sample", 1:10)

# Compute column ranks
ranks <- colranks(X)
print(ranks[1:5, 1:5])
#>       Sample1 Sample2 Sample3 Sample4 Sample5
#> Gene1       3       9       3       5       3
#> Gene2       5       5       7       3       6
#> Gene3       9       6       4      10       1
#> Gene4       6       4       5       9      10
#> Gene5       7       2       6       8       9

# Compute signed ranks
signed_ranks <- colranks(X, signed = TRUE)
print(signed_ranks[1:5, 1:5])
#>       Sample1 Sample2 Sample3 Sample4 Sample5
#> Gene1      -6       8      -7       5      -5
#> Gene2      -3       2      -2      -2      -2
#> Gene3       9       3      -6      10      -9
#> Gene4       1       1      -4       9      10
#> Gene5       2      -6      -3       8       8
```
