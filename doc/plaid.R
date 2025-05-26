## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# devtools::install_github('BigOmics/plaid')

## -----------------------------------------------------------------------------
library("plaid")
load(system.file("extdata", "pbmc3k-50cells.rda", package = "plaid"),verbose=TRUE)
if(FALSE) {
  library(Seurat)
  library(SeuratData)
  data("pbmc3k.final")
  pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)
  X <- pbmc3k.final[['RNA']]@data
  celltype <- pbmc3k.final$seurat_annotations
}
dim(X)

## -----------------------------------------------------------------------------
hallmarks <- system.file("extdata", "hallmarks.gmt", package = "plaid")
gmt <- read.gmt(hallmarks)
matG <- gmt2mat(gmt)
dim(matG)

## -----------------------------------------------------------------------------
gsetX <- plaid(X, matG, normalize=TRUE)
dim(gsetX)

## -----------------------------------------------------------------------------
y <- 1*(celltype == "B")
res <- plaid.test(X, y, matG, gsetX=gsetX, tests=c("one","lm"))

## -----------------------------------------------------------------------------
res <- res[order(res[,"p.meta"]),] 
head(res)

## -----------------------------------------------------------------------------
fc <- res[,"gsetFC"]
pv <- res[,"p.meta"]
plot( fc, -log10(pv), xlab="logFC", ylab="-log10p", pch=19)
abline(h=0, v=0, lty=2)
text( fc[1:5], -log10(pv[1:5]), rownames(res)[1:5],pos=2)

## -----------------------------------------------------------------------------
sing <- replaid.sing(X, matG)

## -----------------------------------------------------------------------------
scse <- replaid.scse(X, matG, removeLog2=TRUE, scoreMean=FALSE)

