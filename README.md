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

## Usage example

We provide a basic example on how to use Plaid. This example uses the
pbmc3k dataset from Seurat which is a dataset of 2,700 PBMC single cells. 
For the gene sets, as example, we included the hallmarks genesets from MSigDB.

However, we invite you to use your own bigger datasets and download
bigger gene set collections as this shows the speed advantage of
plaid. Subsequently, we show how the single-sample scores can be used
for differential enrichment testing.

```r
library("plaid")
library(Seurat)
library(SeuratData)
data("pbmc3k.final")
pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)
X <- pbmc3k.final[['RNA']]@data
dim(X)

hallmarks <- system.file("extdata", "hallmarks.gmt", package = "plaid")
gmt <- read.gmt(hallmarks)
matG <- gmt2mat(gmt)
dim(matG)

## run plaid
gsetX <- plaid(X, matG)
dim(gsetX)

## differential enrichment testing
celltype <- pbmc3k.final$seurat_annotations
y <- (celltype == "B")
res <- plaid.test(X, y, matG, gsetX=gsetX)
head(res)

## simulate other scores
s1 <- replaid.sing(X, matG)
s2 <- replaid.ssgsea(X, matG, alpha=0)
s3 <- replaid.scse(X, matG)
S <- cbind(plaid=gsetX[,1], sing=s1[,1], ssgsea=s2[,1], scSE=s3[,1])
pairs(S)
```

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch
