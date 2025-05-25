

library(Seurat)
library(SeuratData)
SeuratData::InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)
X <- pbmc3k.final[['RNA']]@data
y <- pbmc3k.final$seurat_annotations
table(y)
sel <- c(head(which(y=="B"),25), head(which(y=="Naive CD4 T"),25))
X <- X[,sel]
X <- X[rowSums(X)>0,]
celltype <- c(rep("_Bcell",25),rep("_Tcell",25))
save(X, celltype, file="../inst/extdata/pbmc3k-50cells.rda")
