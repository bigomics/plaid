# Reimplementation of dualGSEA (Bull et al., 2024) but defaults with replaid backend. For the preranked test we still use fgsea. Should be much faster than original using fgsea + GSVA::ssGSEA.

Reimplementation of dualGSEA (Bull et al., 2024) but defaults with
replaid backend. For the preranked test we still use fgsea. Should be
much faster than original using fgsea + GSVA::ssGSEA.

## Usage

``` r
dualGSEA(
  X,
  y,
  G,
  gmt = NULL,
  gsetX = NULL,
  fc.method = c("fgsea", "rankcor", "ztest", "ttest", "cor")[2],
  ss.method = c("plaid", "replaid.ssgsea", "replaid.gsva", "ssgsea", "gsva")[1],
  metap.method = c("stouffer", "fisher", "maxp")[1],
  pv1 = NULL,
  pv2 = NULL,
  sort.by = "p.dual"
)
```

## Arguments

- X:

  Expression matrix with genes on rows and samples on columns

- y:

  Binary vector (0/1) indicating group membership

- G:

  Sparse matrix of gene sets. Non-zero entry indicates gene/feature is
  part of gene sets. Features on rows, gene sets on columns.

- gmt:

  List of gene sets in GMT format

- fc.method:

  Method for fold change testing ("fgsea", "ztest", "ttest", "rankcor",
  "cor")

- ss.method:

  Method for single-sample enrichment ("plaid", "replaid.ssgsea",
  "replaid.gsva", "ssgsea", "gsva")

- metap.method:

  Method for combining p-values ("stouffer", "fisher" or "maxp").
  Default "stouffer".

- pv1:

  Pre-computed p-values from fold change test. If NULL, will be computed
  based on fc.test.

- pv2:

  Pre-computed p-values from single sample test. If NULL, will be
  computed using gset_ttest.

- sort.by:

  Column name to sort results by ("p.dual", "gsetFC", "p.fc", "p.ss").
  Default "p.dual".

## Value

Data frame with columns: gsetFC (gene set fold change), size (gene set
size), p.fc (p-value from fold change test), p.ss (p-value from single
sample test), p.dual (combined p-value), and q.dual (FDR-adjusted
combined p-value).

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
results_cor <- dualGSEA(X, y, gmt, fc.method = "cor", ss.method = "replaid.gsva")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': 'list' object cannot be coerced to type 'double'
print(head(results_cor))
#> Error: object 'results_cor' not found

# \donttest{
# Perform dualGSEA with fgsea (requires fgsea package)
if (requireNamespace("fgsea", quietly = TRUE)) {
  results <- dualGSEA(X, y, gmt, fc.method = "fgsea", ss.method = "replaid.ssgsea")
  print(head(results))
}
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'which': 'list' object cannot be coerced to type 'double'
# }
```
