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

## Usage example

We provide a basic example on how to use Plaid to correct batch
effects (BEs) in a biological dataset.  We use the GSE10846 dataset
(Lenz et al., 2008), which includes array gene expression profiling
data from diffuse large B-cell lymphoma (DLBCL) samples from patients
pre-treated with the pharmacological regimens CHOP and
Rituximab-CHOP. Dataset includes two biological types of DLBCL: ABC
and GCB. Because treatment was performed prior to expression profiling
and samples were split in the two treatment groups, "treatment"
represents a batch variable. We show that BEs in the uncorrected data
appear evident with samples clustering by treatment. Following Plaid
batch correction, the samples cluster by DLBCL type, reflecting their
biological heterogeneity.

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

## Batch correction with Plaid
gsetX <- normalize_medians(gsetX)

celltype <- pbmc3k.final$seurat_annotations
y <- (celltype == "B")
res <- plaid.test(gsetX, y, matG, gsetX=gsetX)
```

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics:

- Antonino Zito, PhD:  antonino.zito@bigomics.ch
- Ivo Kwee, PhD: ivo.kwee@bigomics.ch
