# Fast calculation of scSE score

Calculates Single-Cell Signature Explorer (Pont et al., 2019) scores
using plaid back-end. The computation is 10-100x faster than the
original code.

## Usage

``` r
replaid.scse(
  X,
  matG,
  removeLog2 = NULL,
  scoreMean = FALSE,
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

- removeLog2:

  Logical for whether to remove the Log2, i.e. will apply power
  transform (base2) on input (default TRUE).

- scoreMean:

  Logical for whether computing sum or mean as score (default FALSE).

- assay:

  Character: assay name for Bioconductor objects. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set. Default 5.

- max.genes:

  Integer: maximum genes per gene set. Default 500.

## Value

Matrix of single-sample scSE enrichment scores. Gene sets on rows,
samples on columns.

## Details

Computing the scSE requires running plaid on the linear (not
logarithmic) score and perform additional normalization by the total UMI
per sample. We have wrapped this in a single convenience function:

To replicate the original "sum-of-UMI" scSE score, set `removeLog2=TRUE`
and `scoreMean=FALSE`. scSE and plaid scores become more similar for
`removeLog2=FALSE` and `scoreMean=TRUE`.

We have extensively compared the results from `replaid.scse` and from
the original scSE (implemented in GO lang) and we showed almost
identical results in the score, logFC and p-values.

## Examples

``` r
# Create example expression matrix (log-transformed)
set.seed(123)
X <- log2(matrix(rpois(500, lambda = 10) + 1, nrow = 50, ncol = 10))
rownames(X) <- paste0("GENE", 1:50)
colnames(X) <- paste0("Sample", 1:10)

# Create example gene sets
gmt <- list(
  "Pathway1" = paste0("GENE", 1:15),
  "Pathway2" = paste0("GENE", 10:25)
)
matG <- gmt2mat(gmt)

# Compute scSE scores (original method)
scores <- replaid.scse(X, matG, removeLog2 = TRUE, scoreMean = FALSE)
#> [replaid.scse] Converting data to linear scale (removing log2)...
print(scores[1:2, 1:5])
#>           Sample1  Sample2  Sample3  Sample4  Sample5
#> Pathway2 28.82012 30.90909 36.29764 32.15613 31.62879
#> Pathway1 32.68859 27.81818 31.39746 31.22677 31.81818

# Compute scSE scores (mean method)
scores_mean <- replaid.scse(X, matG, removeLog2 = TRUE, scoreMean = TRUE)
#> [replaid.scse] Converting data to linear scale (removing log2)...
print(scores_mean[1:2, 1:5])
#>            Sample1   Sample2   Sample3   Sample4   Sample5
#> Pathway2 0.7584241 0.8133971 0.9552011 0.8462140 0.8323365
#> Pathway1 0.9080163 0.7727273 0.8721516 0.8674102 0.8838384
```
