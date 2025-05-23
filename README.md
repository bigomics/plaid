# NPM
Latent Batch Effect Correction of Omics data by Nearest-Pair Matching

Batch effects (BEs) are a predominant source of noise in omics data and often mask real biological signals. BEs remain common in existing datasets. Current methods for BE correction mostly rely on specific assumptions or complex models, and may not detect and adjust BEs adequately, impacting downstream analysis and discovery power. To address these challenges we developed NPM, a nearest-neighbor matching-based method that adjusts BEs satisfactorily and outperforms current methods in a wide range of datasets.

NPM (Nearest-Pair Matching) relies on distance-based matching to deterministically search for nearest neighbors with opposite labels, so-called “nearest-pair”, among samples. NPM requires knowledge of the phenotypes but not of the batch assignment. Differently to many other algorithms, NPM does not rely on specific models or underlying distribution. It does not require special experimental designs, randomized controlled experiments, control genes or batch information. NPM is based on the simple rationale that samples sharing a biological state (e.g., phenotype, condition) should empirically pair based on distance in biological profiles, such as transcriptomics profiles.

NPM is freely available on GitHub. It's a main batch correction algorithm in OmicsPlayground, our Bioinformatics platform at BigOmics Analytics. In OmicsPlayground, you can perform NPM without coding needs.

## Installation
You can install the NPM R package with the following steps:
1. Download NPmatch from https://github.com/bigomics/NPM or use "git clone" in the command line;
2. Enter the directory where NPM has been downloaded;
3. In your terminal, type: "R CMD INSTALL NPM" to install NPM.

## Usage example
We provide a basic example on how to use NPmatch to correct batch effects (BEs) in a biological dataset.
We use the GSE10846 dataset (Lenz et al., 2008), which includes array gene expression profiling data from diffuse large B-cell lymphoma (DLBCL) samples from patients pre-treated with the pharmacological regimens CHOP and Rituximab-CHOP. Dataset includes two biological types of DLBCL: ABC and GCB. Because treatment was performed prior to expression profiling and samples were split in the two treatment groups, "treatment" represents a batch variable. We show that BEs in the uncorrected data appear evident with samples clustering by treatment. Following NPmatch batch correction, the samples cluster by DLBCL type, reflecting their biological heterogeneity.

``` r
## Load NPmatch and limma
library("NPmatch")
library("limma")

## X: raw data matrix, with features in rows and samples in columns.
## Meta: matrix or dataframe with the metadata associated with X. 
## We need to ensure that the samples in X and Meta are aligned.
X <- read.table("./data/GSE10846.Expression.txt", sep="\t")
Meta <- read.table("./data/GSE10846.Metadata.txt", sep="\t")
dim(X); class(X)
dim(Meta); class(Meta)
table(rownames(Meta) == colnames(X))

## To correct BEs, NPmatch requires a vector of phenotype labels per sample.
## To assess  BE correction, we will also need a vector of batch labels (see below).
## "pheno": phenotype labels.
## "batch": batch labels.
pheno <- Meta[,"dlbcl.type"]
batch <- Meta[,"Chemotherapy"]

## Intra-sample normalization of the raw data.
## We use the normalize.log2CPM.R function provided
nX <- normalize.log2CPM(X)

## Inter-sample normalization by quantile normalization
nX <- limma::normalizeQuantiles(nX)

## Batch correction with NPmatch
cX <- NPmatch(X=nX, y=pheno, dist.method="cor", sdtop=5000)
table(rownames(Meta) == colnames(cX))

## Check BEs in the raw and batch-corrected data by UMAP or t-SNE
LL <- list(X, cX)
names(LL) <- c("Uncorrected", "Batch-corrected")
Var <- c("Batch", "Pheno")

x11(width = 10, height = 10)
par(mfrow = c(2,2))
i=1
for(i in 1:length(LL)) {
     
      nb <- max(1, min(30, round(ncol(LL[[i]]) / 5)))
      # pos <- Rtsne::Rtsne(t(LL[[i]]), perplexity=nb)$Y
      pos <- uwot::tumap(t(LL[[i]]), n_neighbors = max(2, nb)) 
      
      pos <- data.frame(Dim1=pos[,1], Dim2=pos[,2], Pheno=pheno, Batch=batch)
      table(rownames(pos) == colnames(cX))
      pos[,1:2] <- apply(pos[,1:2], 2, function(x) as.numeric(x))
      pos$Col.Pheno <- as.numeric(factor(pos$Pheno))
      pos$Col.Batch <- as.numeric(factor(pos$Batch))
        
      v=1
      for(v in 1:length(Var)) {
            Col <- pos[,paste0("Col.",Var[v])]
            plot(pos$Dim1,
                 pos$Dim2,
                 col = Col,
                 xlab = "Dim1", 
                 ylab = "Dim2",
                 pch = 18, 
                 cex = 0.8, 
                 cex.lab = 1.3,
                 cex.axis = 1.3,
                 las = 1, 
                 tcl = -0.1,
                 mgp = c(1.5,0.5,0))
            
            mtext(names(LL)[i], 
                  font = 2,
                  adj = 0.5, 
                  cex = 1)
    
            legend("bottomleft",
                   unique(pos[,Var[v]]),
                   cex = 1,
                   bty = "n",
                   fill = unique(Col),
                   col = unique(Col))
    
            grid(lwd = 1.2)
       }
}
```

### Support
For support feel free to reach our Bioinformatics Data Science Team at BigOmics Analytics:

Antonino Zito, PhD:  antonino.zito@bigomics.ch

Ivo Kwee, PhD: ivo.kwee@bigomics.ch
