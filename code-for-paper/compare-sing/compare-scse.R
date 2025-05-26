library(playbase)
library(plaid)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../R/functions.R")
source("../R/datasets.R")

matG <- Matrix::t(playdata::GSETxGENE)
colnames(matG) <- gsub("[/.]","_",colnames(matG))
sel <- grep("HALLMARK",colnames(matG))
sel <- grep("GO_BP",colnames(matG))
#sel <- head(sel,1000)
full.matG <- matG[,sel]

DATASETS <- c("testis50","pbmc3k","geiger","GSE10846-dlbcl",
              "tcga-brca_pub","GSE72056-scmelanoma","GSE102908-ibet")

ds=DATASETS[1]

pdf("compare-plaid-vs-scseSUM.pdf", h=16, w=9, pointsize=14)
par(mfrow=c(7,4), mar=c(3.4,4,2.4,1), mgp=c(2.1,0.8,0))

for(ds in DATASETS) {

  cat("*****",ds,"*****\n")
  dataset <- get_dataset(ds)
  str(dataset)
  cat("name = ",dataset$name,"\n")

  X <- dataset$X
  y <- dataset$y
  cat("ncolX = ",ncol(X),"\n")
  
  gg <- intersect(rownames(X),rownames(full.matG))
  X <- X[gg,]
  matG <- full.matG[gg,]
  matG <- matG[,Matrix::colSums(matG!=0)>10]
  dim(matG)
  dim(X)
  
  gmt <- mat2gmt(matG)
  length(gmt)
  
  ## run ranked plaid and sing
  prepare.SCSE(X, gmt, path="../ext")
#  peakRAM::peakRAM(G1 <- run.SCSE(X, gmt, removeLog2=FALSE, scoreMean=TRUE, path="../ext"))
#  peakRAM::peakRAM(G2 <- replaid.scse(X, matG, removeLog2=FALSE, scoreMean=TRUE))
  peakRAM::peakRAM(G1 <- run.SCSE(X, gmt, removeLog2=1, scoreMean=FALSE, path="../ext"))
  peakRAM::peakRAM(G2 <- replaid.scse(X, matG, removeLog2=1, scoreMean=FALSE))
  gg <- intersect(rownames(G1),rownames(G2))
  G1 <- G1[gg,]
  G2 <- G2[gg,]
  
  m1 <- gx.limma(G1, y, lfc=0, fdr=1, sort.by="none")
  m2 <- gx.limma(G2, y, lfc=0, fdr=1, sort.by="none")
  gg <- intersect(rownames(m1),rownames(m2))
  m1 <- m1[gg,]
  m2 <- m2[gg,]
  
  plot( G1[,1], G2[,1], xlab="scSE score", ylab="replaid.scse")
  title("score", cex.main=1.2)
  legend("topleft", legend=toupper(ds), cex=0.8, bty="n")
  
  plot( m1$logFC, m2$logFC, xlab="scSE score", ylab="replaid.scse" )
  title("logFC", cex.main=1.2)
  plot( m1$P.Value, m2$P.Value, xlab="scSE score", ylab="replaid.scse" )
  title("p-value", cex.main=1.2)
  plot( -log10(m1$P.Value), -log10(m2$P.Value), xlab="scSE score", ylab="replaid.scse" )
  title("-log10p", cex.main=1.2)  
}
dev.off()


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------
