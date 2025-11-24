# Calculate gene set rank correlation

Compute rank correlation between a gene rank vector/matrix and gene sets

## Usage

``` r
gset.rankcor(rnk, gset, compute.p = FALSE, use.rank = TRUE)
```

## Arguments

- rnk:

  Numeric vector or matrix of gene ranks, with genes as row names

- gset:

  Numeric matrix of gene sets, with genes as row/column names

- compute.p:

  Logical indicating whether to compute p-values

- use.rank:

  Logical indicating whether to rank transform rnk before correlation

## Value

Named list with components:

- rho - Matrix of correlation coefficients between rnk and gset

- p.value - Matrix of p-values for correlation (if compute.p = TRUE)

- q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)

## Details

This function calculates sparse rank correlation between rnk and each
column of gset using
[`qlcMatrix::corSparse()`](https://rdrr.io/pkg/qlcMatrix/man/corSparse.html).
It handles missing values in rnk by computing column-wise correlations.

P-values are computed from statistical distribution

## Examples

``` r
# Create example rank vector
set.seed(123)
ranks <- rnorm(100)
names(ranks) <- paste0("GENE", 1:100)

# Create example gene sets as sparse matrix
gmt <- list(
  "Pathway1" = paste0("GENE", 1:20),
  "Pathway2" = paste0("GENE", 15:35),
  "Pathway3" = paste0("GENE", 30:50)
)
genesets <- gmt2mat(gmt)

# Calculate rank correlation
result <- gset.rankcor(ranks, genesets, compute.p = TRUE)
print(result$rho)
#>                  rnk
#> Pathway2 -0.06037224
#> Pathway3  0.17550071
#> Pathway1  0.08486979
print(result$p.value)
#>                rnk
#> Pathway2 0.6751902
#> Pathway3 0.2168031
#> Pathway1 0.5551073
```
