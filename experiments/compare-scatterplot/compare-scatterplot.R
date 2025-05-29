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
sel <- grep("GO_BP",colnames(matG))
#sel <- head(sel,1000)
full.matG <- matG[,sel]
dim(full.matG)

DATASETS <- c("testis50","pbmc3k","geiger","GSE10846-dlbcl",
              "tcga-brca_pub","GSE72056-scmelanoma",
              "GSE102908-ibet")

ds=DATASETS[1]
ds=DATASETS[2]
ds=DATASETS[6]

for(ds in DATASETS) {

  cat("*****",ds,"*****\n")
  dataset <- get_dataset(ds)
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
  
  pdf(paste0("compare-scatterplot-",ds,".pdf"), w=10, h=10)

  xx <- run$results
  xx <- xx[!sapply(xx,is.null)]
  xx <- xx[order(names(xx))]
  names(xx)
  if(1) {
    sel <- grep("replaid",names(xx),invert=TRUE,value=TRUE)
    sel <- setdiff(sel, c("plaid.c","plaid.rc"))
    xx <- xx[sel]
  }
  xx <- lapply(xx, function(x) x[rownames(xx[[1]]),])
  xx <- lapply(xx, plaid::normalize_medians)  ## important
  lm <- lapply(xx, function(x)
    gx.limma(x, y, fdr=1, lfc=0, sort.by='none'))
  lm <- lapply(lm, function(m) m[rownames(lm[[1]]),])
  FF <- sapply(lm, function(m) m$logFC)
  PP <- sapply(lm, function(m) m$P.Value)               
  rownames(FF)=rownames(PP)=rownames(lm[[1]])
  FF[is.na(FF)] <- 0
  PP[is.na(PP)] <- 1
  XX <- do.call(cbind, lapply(xx, function(x) x[,1]))

  rxx <- lapply( xx, function(x) x - rowMeans(x))
  relXX <- do.call(cbind, lapply(rxx, function(x) x[,1]))
  
  gsize <- sapply(gmt,length)
  cex1 <- 0.05*sqrt(gsize[rownames(FF)])
  avgfc <- rowMeans(plaid::colranks(FF,signed=TRUE))
  topfc <- tail(names(sort(avgfc)),50)
  is.topfc <- rownames(XX) %in% topfc

  avgx <- rowMeans(plaid::colranks(XX,signed=FALSE))
  topx <- tail(names(sort(avgx)),50)
  is.topx <- rownames(XX) %in% topx
  col1 <- c("grey50","red2")[1 + 1*is.topfc]
  col1[is.topx] <- "blue"
  col1[is.topfc] <- "red2"  
  
  pairs(XX, main=paste0("geneset score (",ds,") - sample1"),
        cex=cex1, col=col1, gap=0.5)  
  pairs(relXX, main=paste0("relative geneset score (",ds,") - sample1"),
        cex=cex1, col=col1, gap=0.5)  
  pairs(FF, main=paste0("geneset logFC (",ds,")"), cex=cex1, col=col1, gap=0.5)
  pairs(PP, main=paste0("geneset p-value (",ds,")"), pch='.',
        cex=2*cex1, col=col1, gap=0.5)
  pairs(-log(PP), main=paste0("geneset -log(p) (",ds,")"), pch='.' ,
        cex=8*cex1, col=col1, gap=0.5)  
  
  ## focussed FC plot
  sel <- sort(c("cor","rankcor","scse.mean","plaid","plaid.r","gsva","ssgsea","sing"))
  sel <- intersect(sel,colnames(FF))
  F2 <- FF[,sel]
  pairs(F2, cex=1*cex1, col=col1)
  
  dev.off()

  rho.list <- list(
    XX = cor(XX),
    relXX = cor(relXX),
    FF = cor(FF),
    PP = cor(PP),
    logPP = cor(-log(PP))
  )
  save(rho.list, file=paste0("compare-correlation-",ds,".rda"))
}


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------
