# Perform t-test on gene set scores

Perform t-test on gene set scores

## Usage

``` r
gset_ttest(gsetX, y)
```

## Arguments

- gsetX:

  Matrix of gene set scores with gene sets on rows and samples on
  columns

- y:

  Binary vector (0/1) indicating group membership

## Value

Data frame with columns: diff (difference in means), statistic
(t-statistic), pvalue (p-value), and other t-test results.
