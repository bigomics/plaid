library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("functions.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations

X <- pbmc3k.final[['RNA']]@data
dim(X)

matG <- t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))
dim(matG)

gg <- intersect(rownames(X),rownames(matG))
X <- X[gg,]
sel <- grep("GO_BP",colnames(matG))
matG <- matG[gg,sel]
dim(matG)
gmt <- mat2gmt(matG)

G1 <- plaid(X, matG)
G2 <- cor_sparse_matrix(matG, X)
G3 <- cosine_similarity(matG, X)

gg <- cbind( plaid=G1[,1], cor=G2[,1], cosine=G3[,1])
pairs(gg)




## ------------------------------------------------
## --------- TIMEOUT ------------------------------
## ------------------------------------------------

library(R.utils)
dim(X)
length(gmt)
#stop execution after one second
tt <- peakRAM(
  R.utils::withTimeout({
    gsva1 <- gset.gsva(X[,1:5], gmt[1:50], method="gsva")
  }, timeout=1, onTimeout="warning"),
  R.utils::withTimeout({
    gsva2 <- gset.gsva(X[,1:5], gmt[1:50], method="gsva")
  }, timeout=2, onTimeout="warning")
)
tt
dim(gsva)

system.time(local({
  setTimeLimit(elapsed = 1, transient = TRUE)
  gsva <- gset.gsva(X[,1:5], gmt[1:50], method="gsva")
}))

require(R.utils)
for(i in 1:5) {
  tryCatch(
    expr = {
      withTimeout({Sys.sleep(i); cat(i, "\n")}, 
        timeout = 3.1)
    }, 
    TimeoutException = function(ex) cat("Timeout. Skipping.\n")
  )
}
  
