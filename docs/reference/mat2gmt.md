# Convert Binary Matrix to GMT

Convert binary matrix to a GMT (Gene Matrix Transposed) list. The binary
matrix indicates presence or absence of genes in each gene set. Rows
represent genes and columns represent gene sets.

## Usage

``` r
mat2gmt(mat)
```

## Arguments

- mat:

  Matrix with non-zero entries representing genes in each gene set. Rows
  represent genes and columns represent gene sets.

## Value

A list of vector representing each gene set. Each list element
correspond to a gene set and is a vector of genes

## Examples

``` r
# Create example binary matrix
mat <- matrix(0, nrow = 6, ncol = 3)
rownames(mat) <- paste0("GENE", 1:6)
colnames(mat) <- paste0("Pathway", 1:3)
mat[1:3, 1] <- 1  # Pathway1: GENE1, GENE2, GENE3
mat[c(2,4,5), 2] <- 1  # Pathway2: GENE2, GENE4, GENE5
mat[c(1,5,6), 3] <- 1  # Pathway3: GENE1, GENE5, GENE6

# Convert to GMT list
gmt <- mat2gmt(mat)
print(gmt)
#> $Pathway1
#> [1] "GENE1" "GENE2" "GENE3"
#> 
#> $Pathway2
#> [1] "GENE2" "GENE4" "GENE5"
#> 
#> $Pathway3
#> [1] "GENE1" "GENE5" "GENE6"
#> 
```
