# Fast calculation of ssGSEA

Calculates single-sample enrichment singscore (Barbie et al., 2009;
Hänzelmann et al., 2013) using plaid back-end. The computation is
10-100x faster than the original code.

## Usage

``` r
replaid.ssgsea(
  X,
  matG,
  alpha = 0,
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

- alpha:

  Weighting factor for exponential weighting of ranks

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample ssGSEA enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing ssGSEA score requires to compute the ranks of the expression
matrix and weighting of the ranks. We have wrapped this in a single
convenience function.

We have extensively compared the results of `replaid.ssgsea` and from
the original `GSVA` R package and we showed highly similar results in
the score, logFC and p-values. For alpha=0 we obtain exact results, for
alpha\>0 the results are highly similar but not exactly the same.

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

# Compute ssGSEA scores (alpha = 0)
scores <- replaid.ssgsea(X, matG, alpha = 0)
print(scores[1:2, 1:5])
#>            Sample1  Sample2   Sample3   Sample4   Sample5
#> Pathway2 -9.375109 -9.39161 -9.313848 -9.385681 -9.390127
#> Pathway1 -9.307214 -9.36623 -9.361391 -9.380477 -9.338080

# Compute ssGSEA scores with weighting (alpha = 0.25)
scores_weighted <- replaid.ssgsea(X, matG, alpha = 0.25)
print(scores_weighted[1:2, 1:5])
#>            Sample1   Sample2   Sample3   Sample4   Sample5
#> Pathway2 -9.316151 -9.330741 -9.249503 -9.326221 -9.331940
#> Pathway1 -9.249067 -9.308255 -9.298057 -9.317977 -9.280176
```
