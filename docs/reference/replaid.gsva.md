# Fast approximation of GSVA

Calculates single-sample enrichment GSVA (Hänzelmann et al., 2013) using
plaid back-end. The computation is 10-100x faster than the original
code.

## Usage

``` r
replaid.gsva(
  X,
  matG,
  tau = 0,
  rowtf = c("z", "ecdf")[1],
  assay = "logcounts",
  min.genes = 5,
  max.genes = 500
)
```

## Arguments

- X:

  Gene or protein expression matrix. Generally log transformed. See
  details. Genes on rows, samples on columns. Also accepts
  SummarizedExperiment or SingleCellExperiment objects.

- matG:

  Gene sets sparse matrix. Genes on rows, gene sets on columns. Also
  accepts BiocSet objects or GMT lists.

- tau:

  Rank weight parameter (see GSVA publication). Default tau=0.

- rowtf:

  Row transformation method ("z" or "ecdf"). Default "z".

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample GSVA enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing the GSVA score requires to compute the CDF of the expression
matrix, ranking and scoring the genesets. We have wrapped this in a
single convenience function.

We have extensively compared the results of `replaid.gsva` and from the
original `GSVA` R package and we showed good concordance of results in
the score, logFC and p-values.

In the original formulation, GSVA uses an emperical CDF to transform
expression of each feature to a (0;1) relative expression value. For
efficiency reasons, this is here approximated by a z-transform
(center+scale) of each row.

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

# Compute GSVA scores
scores <- replaid.gsva(X, matG)
print(scores[1:2, 1:5])
#>            Sample1   Sample2   Sample3   Sample4   Sample5
#> Pathway2 -9.556889 -9.610941 -9.390508 -9.572427 -9.636507
#> Pathway1 -9.402912 -9.556497 -9.459573 -9.568334 -9.464869
```
