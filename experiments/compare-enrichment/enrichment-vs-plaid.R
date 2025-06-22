library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../../R/plaid.R")
source("../R/functions.R")
source("../R/datasets.R")

##----------------------------------------------------
## Compare with direct methods
##----------------------------------------------------

load(system.file("extdata", "pbmc3k-50cells.rda", package = "plaid"),verbose=TRUE)

X <- X[rowSums(X)>0, ]
dim(X)

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

y <- 1*(celltype == "B")
table(y)

run.enrichment <- function(X, y, matG, gmt) {

  ##X <- normalize_medians(X, ignore.zero=NULL)   
  rX <- apply(X,2,signedRank)
  gsetX <- plaid(X, matG)
  
  de <- gx.limma(as.matrix(X), y, fdr=1, lfc=0)
  sig.up <- rownames(de)[de$logFC > 0.2 & de$P.Value < 0.05]
  sig.dn <- rownames(de)[de$logFC < -0.2 & de$P.Value < 0.05]  
  fc <- Matrix::rowMeans(X[,y==1]) - Matrix::rowMeans(X[,y==0])
  timings <- peakRAM(
    res.fgsea <- fgsea::fgsea(gmt, fc),
    res.fisher <- gset.fisher2(sig.up, sig.dn, gmt, fdr=1,
                               min.genes=0, max.genes=9999),
    res.sing <- sing.limma(X, y, gmt),     

    res.one <- plaid.test(X, y, matG, gsetX, tests="one"),
    res.two <- plaid.test(X, y, matG, gsetX, tests="two"),
    res.lm  <- plaid.test(X, y, matG, gsetX, tests="lm"),
    res.test <- plaid.test(X, y, matG, gsetX ),
    res2.one <- plaid.test(rX, y, matG, gsetX, tests="one"),
    res2.two <- plaid.test(rX, y, matG, gsetX, tests="two"),
    res2.lm  <- plaid.test(rX, y, matG, gsetX, tests="lm"),
    res2.test <- plaid.test(rX, y, matG, gsetX)
  )
  timings  
  res.fgsea <- as.data.frame(res.fgsea)
  rownames(res.fgsea) <- res.fgsea$pathway
  res <- list(
    fgsea = res.fgsea,
    fisher = res.fisher,
    sing.limma = res.sing,    
    plaid.one = res.one,
    plaid.two = res.two,
    plaid.lm = res.lm,
    plaid.test = res.test,
    rplaid.one = res2.one,
    rplaid.two = res2.two,
    rplaid.lm = res2.lm,
    rplaid.test = res2.test
  )
  res <- res[sapply(res, nrow)>0]
  pp <- Reduce(intersect, lapply(res,rownames))
  score <- cbind(
    fgsea = res.fgsea[pp,"NES"],
    fisher = res.fisher[pp,"sign"],
    sing.limma = res.sing[pp,"logFC"],    
    plaid.one = res.one[pp,"gsetFC"],
    plaid.two = res.two[pp,"gsetFC"],
    plaid.lm = res.lm[pp,"gsetFC"],
    plaid.test = res.test[pp,"gsetFC"],
    rplaid.one = res2.one[pp,"gsetFC"],
    rplaid.two = res2.two[pp,"gsetFC"],
    rplaid.lm = res2.lm[pp,"gsetFC"],
    rplaid.test = res2.test[pp,"gsetFC"]
  )  
  pvalue <- cbind(
    fgsea = res.fgsea[pp,"pval"],
    fisher = res.fisher[pp,"p.value"],
    sing.limma = res.sing[pp,"P.Value"],    
    plaid.one = res.one[pp,"p.one"],
    plaid.two = res.two[pp,"p.two"],
    plaid.lm = res.lm[pp,"p.lm"],
    plaid.test = res.test[pp,"p.meta"],
    rplaid.one = res2.one[pp,"p.one"],
    rplaid.two = res2.two[pp,"p.two"],
    rplaid.lm = res2.lm[pp,"p.lm"],
    rplaid.test = res2.test[pp,"p.meta"]
  )  
  timings$Function_Call <- colnames(score)
  timings <- cbind(timings, nsets=ncol(matG), nrow=nrow(X), ncol=ncol(X))
  list(score=score, pvalue=pvalue, timings = timings)
}

dim(matG)
enr <- run.enrichment(X, y, matG, gmt)
str(enr)
enr$timings

P <- enr$pvalue
apply(P,2,function(x) head(names(sort(x))))
F <- enr$score
R <- apply(F, 2, rank)
rnkP <- apply(P, 2, rank)

pdf("enrichment-plaid.pdf",h=10,w=10)
pairs( F, pch='.', cex=1, main="enrichment score")
pairs( R, pch='.', cex=1, main="enrichment score (ranks)")
pairs( P, pch='.', cex=1, main="enrichment p-value")
pairs( -log(P), pch='.', cex=1, main="enrichment p-value")
pairs( rnkP, pch='.', cex=1, main="enrichment p-value (ranks)")
dev.off()


##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------


P1 <- P[,4:7]

gs.size <- colSums(matG!=0)
cex1 <- sqrt(4+gs.size[rownames(P1)])

pairs(P1, cex=0.2*cex1)
pairs(-log(P1), cex=0.2*cex1)

