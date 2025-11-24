# Fast calculation of singscore

Calculates single-sample enrichment singscore (Fouratan et al., 2018)
using plaid back-end. The computation is 10-100x faster than the
original code.

## Usage

``` r
replaid.sing(X, matG, assay = "logcounts", min.genes = 5, max.genes = 500)
```

## Arguments

- X:

  Gene or protein expression matrix. Generally log transformed. See
  details. Genes on rows, samples on columns. Also accepts
  SummarizedExperiment or SingleCellExperiment objects.

- matG:

  Gene sets sparse matrix. Genes on rows, gene sets on columns. Also
  accepts BiocSet objects or GMT lists.

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample singscore enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing the singscore requires to compute the ranks of the expression
matrix. We have wrapped this in a single convenience function.

We have extensively compared the results of `replaid.sing` and from the
original `singscore` R package and we showed identical result in the
score, logFC and p-values.

## Examples

``` r
# Create example expression matrix
set.seed(123)
X <- matrix(rnorm(500), nrow = 50, ncol = 10)
rownames(X) <- paste0("GENE", 1:50)
colnames(X) <- paste0("Sample", 1:10)

# Create example gene sets
gmt <- list(
  "Pathway1" = paste0("GENE", 1:15),
  "Pathway2" = paste0("GENE", 10:25)
)
matG <- gmt2mat(gmt)

# Compute singscore
scores <- replaid.sing(X, matG)
print(scores[1:2, 1:5])
#>              Sample1     Sample2      Sample3    Sample4    Sample5
#> Pathway2 -0.03789474 -0.06315789  0.044210526 0.03368421 0.01684211
#> Pathway1  0.03000000 -0.03777778 -0.003333333 0.03888889 0.06888889
```
