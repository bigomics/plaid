# Reimplementation of dualGSEA (Bull et al., 2024) but defaults with replaid backend. For the preranked test we still use fgsea. Should be much faster than original using fgsea + GSVA::ssGSEA.

Reimplementation of dualGSEA (Bull et al., 2024) but defaults with
replaid backend. For the preranked test we still use fgsea. Should be
much faster than original using fgsea + GSVA::ssGSEA.

## Usage

``` r
dualGSEA(
  X,
  y,
  gmt,
  matG = NULL,
  fc.test = "fgsea",
  gsea.method = "replaid.ssgsea"
)
```

## Arguments

- X:

  Expression matrix with genes on rows and samples on columns

- y:

  Binary vector (0/1) indicating group membership

- gmt:

  List of gene sets in GMT format

- matG:

  Optional sparse matrix of gene sets (will be computed from gmt if
  NULL)

- fc.test:

  Method for fold change testing ("fgsea", "ztest", "rankcor", "cor")

- gsea.method:

  Method for single-sample enrichment ("replaid.ssgsea", "replaid.gsva",
  "ssgsea", "gsva")

## Value

Data frame with results from dual testing including fold changes,
p-values, and combined statistical measures.

## Examples

``` r
# Create example expression matrix
set.seed(123)
X <- matrix(rnorm(1000), nrow = 100, ncol = 20)
rownames(X) <- paste0("GENE", 1:100)
colnames(X) <- paste0("Sample", 1:20)

# Create binary group vector
y <- rep(c(0, 1), each = 10)

# Create example gene sets
gmt <- list(
  "Pathway1" = paste0("GENE", 1:20),
  "Pathway2" = paste0("GENE", 15:35),
  "Pathway3" = paste0("GENE", 30:50)
)

# Perform dualGSEA with correlation test (fast method)
results_cor <- dualGSEA(X, y, gmt, fc.test = "cor", gsea.method = "replaid.gsva")
#> fc.test using cor
#> computing matG
#> single-sample testing using replaid.gsva
print(head(results_cor))
#>          gsetFC size p.fc p.ss p.dual q.dual
#> Pathway2      0   21    1    1      1      1
#> Pathway3      0   21    1    1      1      1
#> Pathway1      0   20    1    1      1      1

# \donttest{
# Perform dualGSEA with fgsea (requires fgsea package)
if (requireNamespace("fgsea", quietly = TRUE)) {
  results <- dualGSEA(X, y, gmt, fc.test = "fgsea", gsea.method = "replaid.ssgsea")
  print(head(results))
}
#> fc.test using fgsea
#> Warning: For some of the pathways the P-values were likely overestimated. For such pathways log2err is set to NA.
#> computing matG
#> single-sample testing using replaid.ssgsea
#>          gsetFC size        p.fc p.ss p.dual q.dual
#> Pathway2      0   21 0.005633866    1      1      1
#> Pathway3      0   21 0.005633866    1      1      1
#> Pathway1      0   20 0.005391312    1      1      1
# }
```
