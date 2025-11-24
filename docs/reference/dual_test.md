# Statistical testing of differentially enrichment

This function performs statistical testing for differential enrichment
using plaid

## Usage

``` r
dual_test(
  X,
  y,
  G,
  gsetX = NULL,
  fc.test = "cor",
  pv1 = NULL,
  pv2 = NULL,
  metap.method = "stouffer",
  sort.by = "p.dual"
)
```

## Arguments

- X:

  Matrix of log expression value

- y:

  Vector of 0s and 1s indicating group

- G:

  Sparse matrix of gene sets. Non-zero entry indicates gene/feature is
  part of gene sets. Features on rows, gene sets on columns.

- gsetX:

  Gene set score matrix which is output of
  [`plaid()`](https://bigomics.github.io/plaid/reference/plaid.md). Can
  be NULL in that case it will be recomputed from X and G (default
  required).

- fc.test:

  Method for fold change testing ("ztest", "rankcor", "cor"). Default
  "cor".

- pv1:

  Pre-computed p-values from fold change test. If NULL, will be computed
  based on fc.test.

- pv2:

  Pre-computed p-values from single sample test. If NULL, will be
  computed using gset_ttest.

- metap.method:

  Method for combining p-values ("stouffer" or "fisher"). Default
  "stouffer".

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
G <- gmt2mat(gmt)

# Perform dual test
results <- dual_test(X, y, G)
#> computing gsetX using plaid. please precompute for efficiency.
print(head(results))
#>          gsetFC size p.fc p.ss p.dual q.dual
#> Pathway2      0   21    1    1      1      1
#> Pathway3      0   21    1    1      1      1
#> Pathway1      0   20    1    1      1      1

# Perform dual test with correlation test
results_cor <- dual_test(X, y, G, fc.test = "cor")
#> computing gsetX using plaid. please precompute for efficiency.
print(head(results_cor))
#>          gsetFC size p.fc p.ss p.dual q.dual
#> Pathway2      0   21    1    1      1      1
#> Pathway3      0   21    1    1      1      1
#> Pathway1      0   20    1    1      1      1
```
