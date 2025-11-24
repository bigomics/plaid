# Convert GMT to Binary Matrix

Convert a GMT file (Gene Matrix Transposed) to a binary matrix, where
rows represent genes and columns represent gene sets. The binary matrix
indicates presence or absence of genes in a gene set.

## Usage

``` r
gmt2mat(
  gmt,
  max.genes = -1,
  ntop = -1,
  sparse = TRUE,
  bg = NULL,
  use.multicore = TRUE
)
```

## Arguments

- gmt:

  List representing the GMT file: each element is a character vector
  representing a gene set.

- max.genes:

  Max number of genes to include in the binary matrix. Default = -1 to
  include all genes.

- ntop:

  Number of top genes to consider for each gene set. Default = -1 to
  include all genes.

- sparse:

  Logical: create a sparse matrix. Default `TRUE`. If `FALSE` creates a
  dense matrix.

- bg:

  Character vector of background gene set. Default `NULL` to consider
  all unique genes.

- use.multicore:

  Logical: use parallel processing ('parallel' R package). Default
  `TRUE`.

## Value

A binary matrix representing the presence or absence of genes in each
gene set. Rows correspond to genes, and columns correspond to gene sets.

## Examples

``` r
# Create example GMT data
gmt <- list(
  "Pathway1" = c("GENE1", "GENE2", "GENE3"),
  "Pathway2" = c("GENE2", "GENE4", "GENE5"),
  "Pathway3" = c("GENE1", "GENE5", "GENE6")
)

# Convert to binary matrix
mat <- gmt2mat(gmt)
print(mat)
#> 6 x 3 sparse Matrix of class "dgCMatrix"
#>       Pathway1 Pathway2 Pathway3
#> GENE1        1        .        1
#> GENE2        1        1        .
#> GENE5        .        1        1
#> GENE3        1        .        .
#> GENE4        .        1        .
#> GENE6        .        .        1

# Create dense matrix instead of sparse
mat_dense <- gmt2mat(gmt, sparse = FALSE)
print(mat_dense)
#>       Pathway1 Pathway2 Pathway3
#> GENE1        1        0        1
#> GENE2        1        1        0
#> GENE5        0        1        1
#> GENE3        1        0        0
#> GENE4        0        1        0
#> GENE6        0        0        1
```
