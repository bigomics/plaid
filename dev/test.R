##library("plaid")
library(devtools)
load_all()

library(Seurat)
library(SeuratData)
data("pbmc3k.final")
pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)
X <- pbmc3k.final[['RNA']]@data
dim(X)




gmt <- read.gmt(system.file("extdata", "hallmarks.gmt", package = "plaid"))
gmt <- read.gmt("~/Playground/data/gmt/c5.all.v6.1.symbols.gmt")
matG <- gmt2mat(gmt)
dim(matG)

## run plaid
gsetX <- plaid(X, matG)
dim(gsetX)

## good practice to normalize median
gsetX <- normalize_medians(gsetX)

## do test
celltype <- pbmc3k.final$seurat_annotations
y <- (celltype == "B")
res <- plaid.test(X, y, G=matG, gsetX=gsetX)
head(res)

pairs(res[,-1])
pairs(-log10(res[,-1]))

head(res[order(res[,2]),1])
head(res[order(res[,3]),1])
head(res[order(res[,4]),1])
head(res[order(res[,5]),1])

fc <- rowMeans(X[,y==1]) - rowMeans(X[,y==0])
res.fgsea <- fgsea::fgsea( gmt, fc)
head(res.fgsea)
nes <- res.fgsea$NES
p.nes <- res.fgsea$pval
names(nes)=names(p.nes)=res.fgsea$pathway

