library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("functions.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
DimPlot(pbmc3k.final, reduction = "umap",
        group.by = "seurat_annotations",label = T) + NoLegend()
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations
dim(pbmc3k.final)

matG <- Matrix::t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))

sel <- grep("HALLMARK",colnames(matG))
sel <- grep("GO_BP",colnames(matG))
#sel <- head(sel,1000)
matG <- matG[,sel]
dim(matG)

X <- pbmc3k.final[['RNA']]@data
X <- X[rowSums(X)>0, ]
X <- X[,1:25]
##X <- logCPM(playbase::COUNTS)

gg <- intersect(rownames(X),rownames(matG))
X <- X[gg,]
matG <- matG[gg,]
matG <- matG[,Matrix::colSums(matG!=0)>10]
dim(matG)
dim(X)

gmt <- mat2gmt(matG)
length(gmt)

rX  <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="min")) / nrow(X)
cX <- as.matrix(X - rowMeans(X, na.rm=TRUE))
rcX <- Matrix::t(matrixStats::colRanks(cX, ties.method="min")) / nrow(X)
xlist <- list(X=X, rX=rX, cX=cX, rcX=rcX)
#xlist <- list(X=X)

dim(X)
length(gmt)
runx <- lapply(xlist, function(x) run.methods(x, gmt))
length(runx)
runx[[1]]$timings

head(pbmc3k.final)
celltype <- pbmc3k.final$seurat_annotations
head(celltype,25)
y <- 1*(celltype[colnames(X)] == "B")
table(y)

samples <- playbase::SAMPLES
y <- 1*(samples$activated=="act")
table(y)

median_normalize <- function(x) t(t(x) - apply(x,2,median)) + median(x)

pdf("compare-pbmc3k-boxplots.pdf", w=12, h=8)
k=1
for(k in 1:length(runx)) {
  rr <- runx[[k]]$results
  par(mfrow=c(4,4), mar=c(3,2,3,1))
  i=1
  for(i in 1:length(rr)) {
    boxplot( rr[[i]], main = names(rr)[i])
  }
  par(mfrow=c(4,4), mar=c(3,2,3,1))
  i=1
  for(i in 1:length(rr)) {
    rx <- rr[[i]]
    rx <- rx - Matrix::rowMeans(rx)
    boxplot( rx, main = names(rr)[i])
  }
}
dev.off()

pdf("compare-pbmc3k-boxplots-normalized.pdf", w=12, h=8)
k=1
for(k in 1:length(runx)) {
  rr <- runx[[k]]$results
  par(mfrow=c(4,4), mar=c(3,2,3,1))
  i=1
  for(i in 1:length(rr)) {
    rx <- median_normalize(rr[[i]])
    boxplot( rx, main = names(rr)[i])
  }
  par(mfrow=c(4,4), mar=c(3,2,3,1))
  i=1
  for(i in 1:length(rr)) {
    rx <- median_normalize(rr[[i]])
    rx <- rx - Matrix::rowMeans(rx)
    boxplot( rx, main = names(rr)[i])
  }
}
dev.off()



##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------
