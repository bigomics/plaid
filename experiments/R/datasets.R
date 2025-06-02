
get_dataset <- function(ds, n=20) {

  pgx.file <- file.path("~/Playground/pgx",paste0(ds,".pgx"))
  has.pgx <- file.exists(pgx.file)
  has.pgx
  out <- NULL
  
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
  
  if(require(playbase) && ds=="geiger") {
    X <- logCPM(playbase::COUNTS)
    dim(X)
    samples <- playbase::SAMPLES
    y <- 1*(samples$activated=="act")
    out <- list(X=X, y=y, name="geiger")
  }

  if(ds=="testis50") {
    fn <- "~/Projects/SingleCellSignatureScorer/data/50first_cells_in_testis.tsv"
    if(!file.exists(fn)) stop("missing data file")
    X <- read.csv(fn, row.names=1, sep='\t')
    X <- t(as.matrix(X))
    topX <- head(X[order(-apply(X,1,sd)),],400)
    topX <- t(scale(t(topX)))
    y <- cutree(hclust(dist(t(topX))),2)
    out <- list(X=X, y=y, name="testis50")
  }

  if(is.null(out)) {
    stop("Error: could not find dataset '", ds, "'")
  }

  ## subsetiting
  if(ncol(out$X) > n) {
    sel <- unlist(tapply(1:length(out$y),out$y,head,n))
    out$y <- out$y[sel]
    out$X <- out$X[,sel]
  }

  return(out)
}


#' @export
logCPM <- function(counts, total = 1e6, prior = 1, log = TRUE) {
  ## Transform to logCPM (log count-per-million) if total counts is
  ## larger than 1e6, otherwise scale to previous avarage total count.
  ##
  ##
  if (is.null(total)) {
    total0 <- mean(Matrix::colSums(counts, na.rm = TRUE)) ## previous sum
    total <- ifelse(total0 < 1e6, total0, 1e6)
    message("[logCPM] setting column sums to = ", round(total, 2))
  }
  if (any(class(counts) == "dgCMatrix")) {
    ## fast/sparse calculate CPM
    cpm <- counts
    cpm[is.na(cpm)] <- 0 ## OK??
    cpm@x <- total * cpm@x / rep.int(Matrix::colSums(cpm), diff(cpm@p)) ## fast divide by columns sum
    if (log) cpm@x <- log2(prior + cpm@x)
    return(cpm)
  } else {
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    ## cpm <- t(t(counts) / totcounts * total)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * total
    if (log) cpm <- log2(prior + cpm)
    return(cpm)
  }
}

