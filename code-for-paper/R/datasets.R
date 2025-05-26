
get_dataset <- function(ds, n=20) {

  pgx.file <- file.path("~/Playground/pgx",paste0(ds,".pgx"))
  has.pgx <- file.exists(pgx.file)
  has.pgx

  if(has.pgx) {
    pgx <- pgx.load(pgx.file)
    X <- pgx$X
    topX <- head(X[order(-apply(X,1,sd)),],400)
    hc <- hclust(dist(t(scale(topX))))
    y <- cutree(hc,2)
    out <- list(X=X, y=y, name=ds)
  }

  if(ds=="pbmc3k") {
    #InstallData("pbmc3k")
    require(Seurat)
    require(SeuratData)
    data("pbmc3k.final")
    pbmc3k.final <- Seurat::UpdateSeuratObject(pbmc3k.final)
    DimPlot(pbmc3k.final, reduction = "umap",
            group.by = "seurat_annotations",label = T) + NoLegend()
    Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations
    dim(pbmc3k.final)
    
    X <- pbmc3k.final[['RNA']]@data
    X <- X[rowSums(X)>0, ]
    celltype <- pbmc3k.final$seurat_annotations
    y <- 1*(celltype == "B")
    out <- list(X=X, y=y, name="pbmc3k")
  }
  
  if(ds=="geiger") {
    X <- logCPM(playbase::COUNTS)
    dim(X)
    samples <- playbase::SAMPLES
    y <- 1*(samples$activated=="act")
    out <- list(X=X, y=y, name="geiger")
  }

  if(ds=="testis50") {
    X <- read.csv("~/Playground/projects/SingleCellSignatureScorer/data/50first_cells_in_testis.tsv", row.names=1, sep='\t')
    X <- t(as.matrix(X))
    topX <- head(X[order(-apply(X,1,sd)),],400)
    topX <- t(scale(t(topX)))
    y <- cutree(hclust(dist(t(topX))),2)
    out <- list(X=X, y=y, name="testis50")
  }

  ## subsetiting
  if( ncol(out$X) > n) {
    sel <- unlist(tapply(1:length(out$y),out$y,head,n))
    out$y <- out$y[sel]
    out$X <- out$X[,sel]
  }

  return(out)
}

