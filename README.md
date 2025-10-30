[![codecov](https://codecov.io/github/bigomics/plaid/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/plaid)

# PLAID

Plaid (Pathway Level Average Intensity Detection) is an ultra-fast
method to compute single-sample enrichment scores for gene expression
or proteomics data. For each sample, plaid computes the gene set score
as the average intensity of the genes/proteins in the gene set. The
output is a gene set score matrix suitable for further analyses.

Plaid is freely available on GitHub. It's a main gene sets scoring
algorithm in OmicsPlayground, our Bioinformatics platform at BigOmics
Analytics. In OmicsPlayground, you can perform Plaid without coding
needs.

## Installation

You can install the Plaid R package with the following steps:
1. Download Plaid from https://github.com/bigomics/plaid or use "git
   clone" in the command line;
2. Enter the directory where Plaid has been downloaded;
3. In your terminal, type: "R CMD INSTALL plaid" to install Plaid.

You can also install plaid from R using devtools with the following
command:

```r
devtools::install_github('bigomics/plaid')
```

Ultimately, it is also possible to install `plaid` on a `conda` environment. For convenience we placed an `environment.yml` file on the repository that will take care of configuring the environment. In order to make use of it, follow this bash commands:

```
conda env create -f environment.yml
conda activate plaid-env
Rscript -e 'remotes::install_github("bigomics/plaid")'
```

## Usage example

We provide a basic example on how to use Plaid. This example uses a
subset of the pbmc3k dataset from Seurat. For the gene sets, as
example, we included the hallmarks genesets from MSigDB.

However, we invite you to use your own bigger datasets and download
bigger gene set collections as this shows the speed advantage of
plaid. Subsequently, we show how the single-sample scores can be used
for differential enrichment testing.

```r
library(plaid)
library(Matrix)
load(system.file("extdata", "pbmc3k-50cells.rda", package = "plaid"),verbose=TRUE)
dim(X)

hallmarks <- system.file("extdata", "hallmarks.gmt", package = "plaid")
gmt <- read.gmt(hallmarks)
matG <- gmt2mat(gmt)
dim(matG)

## run plaid
gsetX <- plaid(X, matG)
dim(gsetX)

## simulate other scores
s1 <- replaid.sing(X, matG)
s2 <- replaid.ssgsea(X, matG, alpha=0)
s3 <- replaid.scse(X, matG, removeLog2=FALSE, scoreMean=TRUE)
S <- cbind(plaid=gsetX[,1], sing=s1[,1], ssgsea=s2[,1], scSE=s3[,1])
pairs(S)

## differential enrichment testing
table(celltype)
y <- 1*(celltype == "B")
res <- dual_test(X, y, matG, gsetX)
head(res)
```

## Using Plaid with Bioconductor objects

Plaid automatically detects and handles Bioconductor classes such as
`SummarizedExperiment`, `SingleCellExperiment`, and `BiocSet` for gene
sets. Simply pass them to any plaid function - no special wrappers needed!

```r
library(plaid)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Create a SingleCellExperiment object
counts <- matrix(rpois(1000, lambda=10), nrow=100, ncol=10)
rownames(counts) <- paste0("GENE", 1:100)
colnames(counts) <- paste0("Cell", 1:10)
sce <- SingleCellExperiment(assays=list(counts=counts))

# Add log-transformed data
logcounts(sce) <- log2(counts(sce) + 1)

# Define gene sets (as GMT list or BiocSet)
gmt <- list(
  "Pathway1" = paste0("GENE", 1:20),
  "Pathway2" = paste0("GENE", 15:35),
  "Pathway3" = paste0("GENE", 30:50)
)

# Just use the regular plaid functions!
# They auto-detect SingleCellExperiment and extract the right assay
scores <- plaid(sce, gmt, assay="logcounts")
dim(scores)

# All other methods work the same way
sing_scores <- replaid.sing(sce, gmt)
ssgsea_scores <- replaid.ssgsea(sce, gmt)
ucell_scores <- replaid.ucell(sce, gmt)
gsva_scores <- replaid.gsva(sce, gmt)
```

**What happens automatically:**
- Detects SummarizedExperiment/SingleCellExperiment and extracts the appropriate assay
- Handles BiocSet objects for gene sets
- Converts GMT lists to sparse matrix format
- Filters gene sets by size (default: 5-500 genes)
- Works with regular matrices too (100% backward compatible)

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch
