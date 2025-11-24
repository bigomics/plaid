# Compute geneset expression as the average log-ration of genes in the geneset. Requires log-expression matrix X and (sparse) geneset matrix matG.

Compute geneset expression as the average log-ration of genes in the
geneset. Requires log-expression matrix X and (sparse) geneset matrix
matG.

## Usage

``` r
gset_averageCLR(X, matG, center = TRUE)
```

## Arguments

- X:

  Log-expression matrix with genes on rows and samples on columns

- matG:

  Sparse gene set matrix with genes on rows and gene sets on columns

- center:

  Logical indicating whether to center the results

- use.rank:

  Logical indicating whether to use rank transformation

## Value

Matrix of gene set expression scores with gene sets on rows and samples
on columns.
