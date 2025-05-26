library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)
##library(plaid)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../R/functions.R")
source("../R/datasets.R")

matG <- Matrix::t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))
sel <- grep("HALLMARK",colnames(matG))
#sel <- grep("GO_BP",colnames(matG))
#sel <- head(sel,1000)
full.matG <- matG[,sel]
dim(full.matG)

DATASETS <- c("pbmc3k","GSE72056-scmelanoma")
ds=DATASETS[1]
ds=DATASETS[2]

for(ds in DATASETS) {

  cat("*****",ds,"*****\n")
  dataset <- get_dataset(ds, n=Inf)
  str(dataset)
  cat("name = ",dataset$name,"\n")

  X <- dataset$X
  y <- dataset$y
  
  gg <- intersect(rownames(X),rownames(full.matG))
  X <- X[gg,]
  matG <- full.matG[gg,]
  matG <- matG[,Matrix::colSums(matG!=0)>10]
  dim(matG)
  dim(X)
  
  gmt <- mat2gmt(matG)
  length(gmt)
  
  ## RUN!!
  run <- run.methods(X, gmt)
  names(run)
  run$timings

  dim(X)
  sv <- irlba::irlba(X, nv=40)$v
  pos <- Rtsne::Rtsne(sv)$Y
  #pos <- uwot::umap(sv)  
  dim(pos)  

  xx <- run$results
  xx <- xx[!sapply(xx,is.null)]
  xx <- xx[order(names(xx))]
  names(xx)
  xx <- lapply(xx, function(x) x[rownames(xx[[1]]),])
  xx <- lapply(xx, plaid::normalize_medians)  ## important

  ## get best geneset
  res.list <- lapply(xx, function(x) gx.limmaF(x, y, fdr=1, lfc=0, sort.by='none'))
  P <- sapply(res.list, function(x) x$P.Value)
  rownames(P) <- rownames(res.list[[1]])
  jj <- order(rowMeans(log10(P)))
  topgs <- rownames(P)[jj][1]
  topgs
  
  pdf(paste0("compare-umap-",ds,".pdf"), w=10, h=10)
  par(mfrow=c(4,4), mar=c(2,2,2,2))
  m=names(xx)[1]
  plot.new(); text(0.5,0.5,topgs)
  for(m in names(xx)) {
    x1 <- xx[[m]][topgs,]
    pal <- colorRampPalette(c("blue", "grey80","red"))(64)
    col1 <- pal[1 + 63*(x1-min(x1))/(max(x1)-min(x1))]
    plot( pos, cex=0.2, col=col1, pch=20, main=m)
  }
  dev.off()
  
}


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------
