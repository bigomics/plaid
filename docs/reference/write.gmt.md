# Write GMT File

Write gene sets to GMT file (Gene Matrix Transposed). The GMT format is
commonly used to store gene sets or gene annotations.

## Usage

``` r
write.gmt(gmt, file, source = NA)
```

## Arguments

- gmt:

  A list of gene sets in GMT format: each gene set is represented as a
  vector of gene names.

- file:

  The file path to write the GMT file.

- source:

  A character vector specifying the source of each gene set. If not
  provided, the names of the gene sets are used as the source.

## Value

Does not return anything.

## Examples

``` r
# Create example GMT data
gmt <- list(
  "Pathway1" = c("GENE1", "GENE2", "GENE3"),
  "Pathway2" = c("GENE2", "GENE4", "GENE5"),
  "Pathway3" = c("GENE1", "GENE5", "GENE6")
)

# \donttest{
# Write to GMT file (creates file in temp directory)
temp_file <- tempfile(fileext = ".gmt")
write.gmt(gmt, temp_file)

# Write with custom source information
temp_file2 <- tempfile(fileext = ".gmt")
write.gmt(gmt, temp_file2, source = c("DB1", "DB2", "DB3"))

# Clean up
unlink(c(temp_file, temp_file2))
# }
```
