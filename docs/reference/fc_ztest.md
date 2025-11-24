# Z-test statistical testing of differentially enrichment

This function performs statistical testing for differential enrichment
using plaid

## Usage

``` r
fc_ztest(fc, G, zmat = FALSE, alpha = 0.5)
```

## Arguments

- fc:

  Vector of logFC values

- G:

  Sparse matrix of gene sets. Non-zero entry indicates gene/feature is
  part of gene sets. Features on rows, gene sets on columns.

- zmat:

  Logical indicating to return z-matrix

- alpha:

  Scalar weight for SD estimation. Default 0.5.

## Value

List with element: z_statistic (z-statistic from one-sample z-test),
p_value (p-value from one-sample z-test), and zmat (z-matrix).
