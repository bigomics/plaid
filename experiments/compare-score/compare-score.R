library(plaid)

BiocManager::install(c("Seurat","SeuratData"))
InstallData("pbmc3k")

##source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../R/functions.R")
source("../R/datasets.R")
source("../../R/plaid.R")

matG <- Matrix::t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))
sel <- grep("HALLMARK",colnames(matG))
sel <- grep("GO_BP",colnames(matG))
#sel <- head(sel,1000)
full.matG <- matG[,sel]

DATASETS <- c("testis50","pbmc3k","geiger","GSE10846-dlbcl",
              "tcga-brca_pub","GSE72056-scmelanoma","GSE102908-ibet")
ds=DATASETS[2]
ds=DATASETS[1]

all.methods <- c("scse","scse.mean","sing","ssgsea")
method="sing"
method="ssgsea"

pdf("compare-score-n20-sd01-REL.pdf", h=15, w=12, pointsize=12)
par(mfrow=c(7,6), mar=c(3.4,4,2.4,1), mgp=c(2.1,0.8,0))
for(ds in DATASETS) {
  
  cat("*****",ds,"*****\n")
  dataset <- try(get_dataset(ds, n=Inf))
  if(inherits(dataset,"try-error") || is.null(dataset)) next

  str(dataset)
  cat("name = ",dataset$name,"\n")
  
  X <- dataset$X
  y <- dataset$y
  dim(X)
  X <- X[,1:min(20,ncol(X))]
  X <- X[mat.rowsds(X)>0.1,]  # important!
  dim(X)
  
  gg <- intersect(rownames(X),rownames(full.matG))
  X <- X[gg,]
  matG <- full.matG[gg,]
  matG <- matG[,Matrix::colSums(matG!=0)>10]
  dim(matG)
  dim(X)

#  sdx <- matrixStats::rowSds(X)
#  mx <- Matrix::rowMeans(X)  
#  avg.sdx <- (t(matG!=0) %*% sdx / colSums(matG!=0))[,1]
#  avg.mx <- (t(matG!=0) %*% mx / colSums(matG!=0))[,1]  
    
  gmt <- mat2gmt(matG)
  length(gmt)
  
  ## run ranked plaid and sing
  time0=time1=NULL
  rel=TRUE
  
  S1 <- gset.singscore(X, gmt, return.score = TRUE)
  S2 <- replaid.sing(X, matG) 
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))
  plot( S1[jj,1], S2[jj,1], xlab="singscore", ylab="replaid.sing")
  title("singscore")
  legend("topleft", legend=toupper(ds), cex=1, bty="n")
  
  S1 <- run.ssgsea(X, gmt, alpha=0) 
  S2 <- replaid.ssgsea(X, matG, alpha=0)
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))  
  plot( S1[jj,1], S2[jj,1], xlab="ssgsea", ylab="replaid.ssgsea")
  title("ssGSEA (alpha=0)")
  
  ##S1 <- gset.gsva(X, gmt) 
  gsvapar <- GSVA::gsvaParam(X, gmt, tau=0, maxDiff=TRUE)
  S1 <- GSVA::gsva(gsvapar, verbose = FALSE)
  S2 <- replaid.gsva(X, matG, tau=0)
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))
  plot( S1[jj,1], S2[jj,1], xlab="GSVA", ylab="replaid.gsva")
  title("GSVA (tau=0)")

  #plot( S1[jj,1], S2[jj,1], cex=0.1+5*avg.sdx[jj]**2 )
  #plot( S1[jj,1], S2[jj,1], cex=0.1+2*avg.mx[jj]**2 )  

  S1 <- t(UCell::ScoreSignatures_UCell(X, gmt)) 
  rownames(S1) <- sub("_UCell$","",rownames(S1))
  S2 <- replaid.ucell(X, matG)
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))  
  plot( S1[jj,1], S2[jj,1], xlab="UCell", ylab="replaid.ucell")
  title("UCell")

  aucMaxRank = ceiling(0.05 * nrow(X))
  S1 <- AUCell::getAUC(AUCell::AUCell_run(X, gmt, aucMaxRank=aucMaxRank))  
  S2 <- replaid.aucell(X, matG, aucMaxRank=aucMaxRank )
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))    
  plot( S1[jj,1], S2[jj,1], xlab="AUCell", ylab="replaid.aucell")
  title("AUCell")  

  prepare.SCSE(X, gmt, path="../scse")
  S1 <- run.SCSE(X, gmt, removeLog2=TRUE, scoreMean=FALSE, path="../scse")
  S2 <- replaid.scse(X, matG, removeLog2=TRUE, scoreMean=FALSE)
  if(rel) {
    S1 <- S1-rowMeans(S1)
    S2 <- S2-rowMeans(S2)    
  }
  jj <- intersect(rownames(S1),rownames(S2))
  plot( S1[jj,1], S2[jj,1], xlab="scSE", ylab="replaid.scse")
  title("scSE")
  
}
dev.off()


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------
