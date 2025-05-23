# PLAID: Pathway Level Average Intensity Detection

Plaid (Pathway Level Average Intensity Detection) is an ultra-fast
method to compute single-sample enrichment scores for gene expression
or proteomics data. For each sample, plaid computes the gene set score
as the average intensity of the genes/proteins in the gene set. The
output is a gene set score matrix suitable for further analyses.

Plaid is freely available on GitHub. It's a main batch correction
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

``` r
# Load Plaid 
library("plaid")

## X: raw data matrix, with features in rows and samples in columns.
## Meta: matrix or dataframe with the metadata associated with X. 
## We need to ensure that the samples in X and Meta are aligned.
X <- read.table("./data/GSE10846.Expression.txt", sep="\t")

gmt <- read.gmt("hallmark.gmt")
matG <- gmt2mat(gmt)

## run plaid
gsetX <- plaid(X, matG)

## Batch correction with Plaid
gsetX <- normalize_medians(gsetX)

y <- sample(1:2, ncol(X))
res <- rowttest(gsetX, y)
```

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics:

- Antonino Zito, PhD:  antonino.zito@bigomics.ch
- Ivo Kwee, PhD: ivo.kwee@bigomics.ch
