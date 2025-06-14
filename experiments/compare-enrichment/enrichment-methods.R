library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)
library(plaid)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../R/functions.R")
source("../R/datasets.R")

##----------------------------------------------------
## Compare with direct methods
##----------------------------------------------------

DATASETS <- c("testis50","pbmc3k","geiger","GSE10846-dlbcl",
              "tcga-brca_pub","GSE72056-scmelanoma","GSE102908-ibet")
ds=DATASETS[2]

run.enrichment_methods <- function(X, y, matG, gmt) {
  y <- 1*(y==y[1])
  de <- gx.limma(as.matrix(X), y, fdr=1, lfc=0)
  sig.up <- rownames(de)[de$logFC > 0.2 & de$P.Value < 0.05]
  sig.dn <- rownames(de)[de$logFC < -0.2 & de$P.Value < 0.05]  
  fc <- Matrix::rowMeans(X[,y==1]) - Matrix::rowMeans(X[,y==0])
  gsetX <- plaid(X, matG)
  res <- list()
  timings <- peakRAM(
    res[['fgsea']] <- fgsea::fgsea(gmt, fc),
    res[['blitz']] <- run.blitzgsea( fc, gmt ),
    res[['fisher']] <- gset.fisher2(sig.up, sig.dn, gmt, fdr=1, min.genes=0, max.genes=Inf),
    res[['cor']] <- gset.rankcor(fc, matG, compute.p=TRUE, use.rank=FALSE),
    res[['rankcor']] <- gset.rankcor(fc, matG, compute.p=TRUE),
    ##res[['avgFC']] <- run.avgFC(fc, matG),
    res[['plaid']] <- plaid.test(X, y, matG, gsetX=gsetX),
    res[['gsva']] <- gsva.limma(X, y, gmt, method="gsva"),
    res[['ssgsea']] <- gsva.limma(X, y, gmt, method="ssgsea"),    
    res[['sing']] <- sing.limma(X, y, gmt) 
  )
  timings  
  res[['fgsea']] <- as.data.frame(res[['fgsea']])
  rownames(res[['fgsea']]) <- res[['fgsea']]$pathway
  res <- res[!sapply(res, is.null)]
  rr <- sapply(res,rownames)
  pp <- Reduce(intersect, rr[!sapply(rr,is.null)])
  score <- cbind(
    fgsea = res[['fgsea']][pp,"NES"],
    blitz.fgsea = res[['blitz']][pp,"nes"],    
    fisher = log(res[['fisher']][pp,"odd.ratio"]) * sign(res[['fisher']][pp,"sign"]),
    cor = res[['cor']]$rho[pp,1],
    rankcor = res[['rankcor']]$rho[pp,1],
    ##avgFC = res[['avgFC']][pp,"logFC"],
    plaid.test = res[['plaid']][pp,"gsetFC"],
    gsva = res[['gsva']][pp,"logFC"],
    ssgsea = res[['ssgsea']][pp,"logFC"],    
    sing = res[['sing']][pp,"logFC"]    
  )
  pvalue <- cbind(
    fgsea = res[['fgsea']][pp,"pval"],
    blitz.fgsea = res[['blitz']][pp,"pval"],        
    fisher = res[['fisher']][pp,"p.value"],
    cor = res[['cor']]$p.value[pp,1],
    rankcor = res[['rankcor']]$p.value[pp,1],
    ##avgFC = res[['avgFC']][pp,"p.value"],
    plaid.one = res[['plaid']][pp,"p.one"],
    plaid.two = res[['plaid']][pp,"p.two"],
    plaid.lm = res[['plaid']][pp,"p.lm"],    
    plaid.meta = res[['plaid']][pp,"p.meta"],    
    gsva = res[['gsva']][pp,"P.Value"],
    ssgsea = res[['ssgsea']][pp,"P.Value"],
    sing = res[['sing']][pp,"P.Value"]    
  )  
  timings$Function_Call <- colnames(score)
  timings <- cbind(timings, nsets=ncol(matG), nrow=nrow(X), ncol=ncol(X))
  list(score=score, pvalue=pvalue, timings = timings)
}

for(ds in DATASETS) {
  
  cat("*****",ds,"*****\n")
  dataset <- get_dataset(ds)
  str(dataset)
  cat("name = ",dataset$name,"\n")
  
  X <- dataset$X
  y <- dataset$y
  cat("ncolX = ",ncol(X),"\n")
  cat("len.y = ",length(y),"\n")
    
  matG <- Matrix::t(playdata::GSETxGENE)
  sel <- grep("GO_BP",colnames(matG))
  matG <- matG[,sel]
  dim(matG)
  ##colnames(matG) <- gsub("[_ /]","-",colnames(matG))

  gg <- intersect(rownames(X),rownames(matG))
  X <- X[gg,]
  matG <- matG[gg,]
  
  gmt <- mat2gmt(matG)
  length(gmt)
  matG <- matG[,names(gmt)]
  dim(matG)
  
  enr <- run.enrichment_methods(X, y, matG, gmt)
  str(enr)
  enr$timings

  P <- enr$pvalue
  apply(P,2,function(x) head(names(sort(x))))
  F <- enr$score
  R <- apply(F, 2, rank)
  rnkP <- apply(P, 2, rank)  
  
  pdf(paste0("enrichment-methods-",ds,".pdf"),h=10,w=10)
  pairs( F, pch='.', cex=2, main="enrichment score")
  pairs( R, pch='.', cex=1, main="enrichment score (ranks)")
  pairs( P, pch='.', cex=1, main="enrichment p-value")  
  pairs( rnkP, pch='.', cex=1, main="enrichment p-value (ranked)")
  pairs( -log(P), pch='.', cex=3, main="enrichment -log(p)")
  dev.off()
  
}


##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------


