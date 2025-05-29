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
y <- 1*(celltype == "B")
res <- plaid.test(X, y, matG, gsetX=gsetX)
head(res)

## plot top score on UMAP
pbmc3k.final[["enrichment"]] <- CreateAssayObject(data = gsetX)
DefaultAssay(pbmc3k.final) <- "enrichment"
p1 <- DimPlot(pbmc3k.final, reduction = "umap", group.by="seurat_annotations")
top.gs <- gsub("_","-",rownames(res)[1])
p2 <- FeaturePlot(pbmc3k.final, reduction = "umap", features = top.gs)
p1 + p2

