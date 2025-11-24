# Statistical testing of differentially enrichment

This function performs statistical testing for differential enrichment
using plaid

## Usage

``` r
fc_ttest(fc, G, sort.by = "pvalue")
```

## Arguments

- fc:

  Vector of logFC values

- G:

  Sparse matrix of gene sets. Non-zero entry indicates gene/feature is
  part of gene sets. Features on rows, gene sets on columns.

- sort.by:

  Column name to sort results by ("pvalue", "gsetFC", or "none")

## Value

Data frame with columns: gsetFC (gene set fold change), pvalue (p-value
from one-sample t-test), and qvalue (FDR-adjusted p-value).
