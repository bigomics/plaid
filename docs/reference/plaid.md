# Compute PLAID single-sample enrichment score

Compute single-sample geneset expression as the average log-expression f
genes in the geneset. Requires log-expression matrix X and (sparse)
geneset matrix matG. If you have gene sets as a gmt list, please convert
it first using the function
[`gmt2mat()`](https://bigomics.github.io/plaid/reference/gmt2mat.md).

## Usage

``` r
plaid(
  X,
  matG,
  stats = c("mean", "sum"),
  chunk = NULL,
  normalize = TRUE,
  nsmooth = 3,
  assay = "logcounts",
  min.genes = 5,
  max.genes = 500
)
```

## Arguments

- X:

  Log-transformed expr. matrix. Genes on rows, samples on columns. Also
  accepts SummarizedExperiment or SingleCellExperiment objects.

- matG:

  Gene sets sparse matrix. Genes on rows, gene sets on columns. Also
  accepts BiocSet objects or GMT lists (named list of gene vectors).

- stats:

  Score computation stats: mean or sum of intensity. Default 'mean'.

- chunk:

  Logical: use chunks for large matrices. Default 'NULL' for autodetect.

- normalize:

  Logical: median normalize results or not. Default 'TRUE'.

- nsmooth:

  Smoothing parameter for more stable average when stats="mean". Default
  3.

- assay:

  Character: assay name to extract from
  SummarizedExperiment/SingleCellExperiment. Default "logcounts".

- min.genes:

  Integer: minimum genes per gene set (for BiocSet/GMT input). Default
  5.

- max.genes:

  Integer: maximum genes per gene set (for BiocSet/GMT input). Default
  500.

## Value

Matrix of single-sample enrichment scores. Gene sets on rows, samples on
columns.

## Details

PLAID needs the gene sets as sparse matrix. If you have your collection
of gene sets a a list, we need first to convert the gmt list to matrix
format.

We recommend to run PLAID on the log transformed expression matrix, not
on the counts, as the average in the logarithmic space is more robust
and is in concordance to calculating the geometric mean.

It is not necessary to normalize your expression matrix before running
PLAID because PLAID performs median normalization of the enrichment
scores afterwards.

It is recommended to use sparse matrix as PLAID relies on sparse matrix
computations. But, PLAID is also fast for dense matrices.

PLAID can also be run on the ranked matrix. This corresponds to the
singscore (Fouratan et al., 2018). PLAID can also be run on the
(non-logarithmic) counts which can be used to calculate the scSE score
(Pont et al., 2019).

PLAID is fast and memery efficient because it uses efficient sparse
matrix computation. When input matrix is very large, PLAID performs
'chunked' computation by splitting the matrix in chunks.

Although `X` and `matG` are generally sparse, the result matrix `gsetX`
generally is dense and can thus be very large. Example: computing gene
set scores for 10K gene sets on 1M cells will create a 10K x 1M dense
matrix which requires ~75GB memory.

PLAID now automatically detects and handles Bioconductor objects. If X
is a SummarizedExperiment or SingleCellExperiment, it will extract the
appropriate assay. If matG is a BiocSet object or GMT list, it will be
converted to sparse matrix format automatically.

## Examples

``` r
library(plaid)

# Create example expression matrix
set.seed(123)
X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(X) <- paste0("GENE", 1:100)
colnames(X) <- paste0("Sample", 1:10)

# Create example gene sets
gmt <- list(
  "Pathway1" = paste0("GENE", 1:20),
  "Pathway2" = paste0("GENE", 15:35),
  "Pathway3" = paste0("GENE", 30:50)
)
matG <- gmt2mat(gmt)

# Compute PLAID scores
gsetX <- plaid(X, matG)
print(dim(gsetX))
#> [1]  3 10
print(gsetX[1:3, 1:5])
#>            Sample1   Sample2   Sample3   Sample4   Sample5
#> Pathway2 -12.01495 -11.85949 -12.04643 -11.73753 -11.82438
#> Pathway3 -11.76847 -11.91129 -12.06332 -12.04601 -12.11773
#> Pathway1 -11.82811 -12.03737 -11.90176 -11.99814 -12.07793

# Use sum statistics instead of mean
gsetX_sum <- plaid(X, matG, stats = "sum")

# \donttest{
# Using real data (if available in package)
extdata_path <- system.file("extdata", "pbmc3k-50cells.rda", package = "plaid")
if (file.exists(extdata_path)) {
  load(extdata_path)
  hallmarks <- system.file("extdata", "hallmarks.gmt", package = "plaid")
  gmt <- read.gmt(hallmarks)
  matG <- gmt2mat(gmt)
  gsetX <- plaid(X, matG)
}
# }
```
