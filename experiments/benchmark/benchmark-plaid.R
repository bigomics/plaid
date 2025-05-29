library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("R/plaid.R")
source("R/functions.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations

X <- pbmc3k.final[['RNA']]@data
dim(X)

matG <- Matrix::t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))
dim(matG)

gg <- intersect(rownames(X),rownames(matG))
X <- X[gg,]
matG <- matG[gg,]
dim(matG)

dim(X)
XL <- do.call(cbind, rep(list(X), 10))
XL <- do.call(cbind, rep(list(XL), 10))
XL <- do.call(cbind, rep(list(XL), 4))
dim(XL)

gmt <- mat2gmt(matG)
gmt <- head(gmt,50000)
length(gmt)
matG <- matG[,names(gmt)]

## 12010x1055200 takes 570 seconds...
peakRAM(rX <- colranks(XL[,], keep.zero=TRUE))
dim(rX)

## 50k gsets takes 50.914 seconds...
peakRAM(gmt2mat(gmt))

##----------------------------------------------------
## Computation vs. number of cells
##----------------------------------------------------
timings <- c()
nlist <- c(1000,100000,200000,400000,600000,800000,1000000)
i=1
n=100
for(n in nlist) {
  X1 <- XL[,1:n]
  G1 <- matG[,1:1000]
  message("----------------------")
  message("num.samples = ", ncol(X1))
  message("num.gsets = ", ncol(G1))
  tt <- peakRAM::peakRAM( plaid(X1, matG=G1, chunk=NULL) )
  tt$Total_RAM_Used_MiB <- NULL
  tt <- cbind(tt, nsets=ncol(G1), nrow=nrow(X1), ncol=ncol(X1))
  timings <- rbind(timings, tt)
  write.csv(timings, file="benchmark-plaid.csv")
}
timings

##----------------------------------------------------
## Computation vs. geneset collection size
##----------------------------------------------------

glist <- c(100,1000,10000,20000,30000,40000,50000)
i=1
g=100
for(g in glist) {
  X1 <- XL[,1:1000]
  G1 <- matG[,1:g]
  message("----------------------")
  message("num.samples = ", ncol(X1))
  message("num.gsets = ", ncol(G1))
  tt <- peakRAM::peakRAM( plaid(X1, matG=G1, chunk=NULL) )
  tt$Total_RAM_Used_MiB <- NULL
  tt <- cbind(tt, nsets=ncol(G1), nrow=nrow(X1), ncol=ncol(X1))
  timings <- rbind(timings, tt)
  write.csv(timings, file="benchmark-plaid.csv")
}
timings


if(0) {

  timings <- read.csv(file="benchmark-plaid.csv", row.names=1)
  timings
  
  pdf("benchmark-plaid.pdf",w=12,h=7.5,pointsize=12)

  par(mfrow=c(2,2), mar=c(4,7,3,1), mgp=c(2.4,1,0))
  sel <- 1:7
  xx <- timings$ncol[sel]
  tt <- timings$Elapsed_Time_sec[sel]
  plot(xx, tt, type="b", cex=2, pch=19,
    xlab="number of samples/cells", ylab="Runtime (seconds)")
  title("Runtime vs. number of cells (1000 genesets)",line=1, cex.main=1.4)
  
  sel <- 8:14
  xx <- timings$nsets[sel]
  tt <- timings$Elapsed_Time_sec[sel]
  plot(xx, tt, type="b", cex=2, pch=19,
    xlab="number of gene sets", ylab="Runtime (seconds)")
  title("Runtime vs. number of gene sets (1000 cells)",line=1, cex.main=1.4)


  sel <- 1:7
  xx <- timings$ncol[sel]
  tt <- timings$Peak_RAM_Used_MiB[sel]
  plot(xx, tt, type="b", cex=2, pch=19,
    xlab="number of samples/cells", ylab="Peak memory (Mb)")
  title("Peak memory vs. number of cells (1000 genesets)",line=1, cex.main=1.4)
  
  sel <- 8:14
  xx <- timings$nsets[sel]
  tt <- timings$Peak_RAM_Used_MiB[sel]
  plot(xx, tt, type="b", cex=2, pch=19,
    xlab="number of gene sets", ylab="Peak memory (Mb)")
  title("Peak memory vs. number of gene sets (1000 cells)",line=1, cex.main=1.4)
  
  dev.off()

  
  
}
