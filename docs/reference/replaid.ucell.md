# Fast calculation of UCell

Calculates single-sample enrichment UCell (Andreatta et al., 2021) using
plaid back-end. The computation is 10-100x faster than the original
code.

## Usage

``` r
replaid.ucell(
  X,
  matG,
  rmax = 1500,
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

- rmax:

  Rank threshold (see Ucell paper). Default rmax = 1500.

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample UCell enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing ssGSEA score requires to compute the ranks of the expression
matrix and truncation of the ranks. We have wrapped this in a single
convenience function.

We have extensively compared the results of `replaid.ucell` and from the
original `UCell` R package and we showed near exacct results in the
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

# Compute UCell scores (default rmax = 1500)
scores <- replaid.ucell(X, matG)
print(scores[1:2, 1:5])
#>           Sample1  Sample2  Sample3  Sample4  Sample5
#> Pathway2 1.055790 1.055309 1.058114 1.057417 1.056583
#> Pathway1 1.057866 1.055968 1.056342 1.057404 1.058131

# Compute UCell scores with custom rmax
scores_custom <- replaid.ucell(X, matG, rmax = 1000)
print(scores_custom[1:2, 1:5])
#>           Sample1  Sample2  Sample3  Sample4  Sample5
#> Pathway2 1.083685 1.082964 1.087171 1.086126 1.084875
#> Pathway1 1.086799 1.083953 1.084513 1.086106 1.087196
```
