# Fast calculation of AUCell

Calculates single-sample enrichment AUCell (Aibar et al., 2017) using
plaid back-end. The computation is 10-100x faster than the original
code.

## Usage

``` r
replaid.aucell(
  X,
  matG,
  aucMaxRank = NULL,
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

- aucMaxRank:

  Rank threshold (see AUCell paper). Default aucMaxRank = 0.05\*nrow(X).

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample AUCell enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing the AUCell score requires to compute the ranks of the
expression matrix and approximating the AUC of a gene set. We have
wrapped this in a single convenience function.

We have extensively compared the results of `replaid.aucell` and from
the original `AUCell` R package and we showed good concordance of
results in the score, logFC and p-values.

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

# Compute AUCell scores
scores <- replaid.aucell(X, matG)
print(scores[1:2, 1:5])
#>            Sample1   Sample2   Sample3   Sample4   Sample5
#> Pathway2 -9.368685 -9.352829 -9.351168 -9.357331 -9.416749
#> Pathway1 -9.386580 -9.390724 -9.370116 -9.392067 -9.356749
```
