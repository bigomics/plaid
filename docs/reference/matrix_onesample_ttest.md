# Perform one-sample t-test on matrix with gene sets

Perform one-sample t-test on matrix with gene sets

## Usage

``` r
matrix_onesample_ttest(Fm, G)
```

## Arguments

- Fm:

  Vector of feature values (e.g., fold changes)

- G:

  Sparse matrix of gene sets with genes on rows and gene sets on columns

## Value

List containing mean, t-statistic, and p-value matrices.
