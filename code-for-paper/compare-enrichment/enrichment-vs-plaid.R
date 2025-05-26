library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../plaid/R/plaid.R")
source("R/functions.R")
source("R/datasets.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)

head(pbmc3k.final)

DimPlot(pbmc3k.final, reduction = "umap",
  group.by = "seurat_annotations",
  label = TRUE, label.size=8) + NoLegend()
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations
dim(pbmc3k.final)

##----------------------------------------------------
## Compare with direct methods
##----------------------------------------------------

X <- pbmc3k.final[['RNA']]@data
X <- X[rowSums(X)>0, ]
dim(X)

matG <- Matrix::t(playdata::GSETxGENE)
sel <- grep("GO_BP",colnames(matG))
matG <- matG[,sel]
dim(matG)
##colnames(matG) <- gsub("[_ /]","-",colnames(matG))

gg <- intersect(rownames(X),rownames(matG))
X <- X[gg,1:25]
matG <- matG[gg,]

gmt <- mat2gmt(matG)
length(gmt)
matG <- matG[,names(gmt)]
dim(matG)

y <- 1*(pbmc3k.final$seurat_annotations == "B")
y <- y[1:ncol(X)]
table(y)

gsetX <- plaid(X, matG)
res <- plaid.test(X, y, matG, gsetX=gsetX, tests=c("one","two","lm"))
head(res)
head(res[,"p.meta"],6)

pairs( res[,-1] )
pairs( -log(res[,-1]) )

run.enrichment <- function(X, y, matG, gmt) {
  ##X <- normalize_medians(X, ignore.zero=NULL)   
  de <- gx.limma(as.matrix(X), y, fdr=1, lfc=0)
  sig.up <- rownames(de)[de$logFC > 0.2 & de$P.Value < 0.05]
  sig.dn <- rownames(de)[de$logFC < -0.2 & de$P.Value < 0.05]  
  fc <- Matrix::rowMeans(X[,y==1]) - Matrix::rowMeans(X[,y==0])
  timings <- peakRAM(
    res.fgsea <- fgsea::fgsea(gmt, fc),
    res.fisher <- gset.fisher2(sig.up, sig.dn, gmt, fdr=1,
                               min.genes=0, max.genes=9999),
    res.limma <- plaid.limma(X, y, matG),
    res.tt <- plaid.ttest(X, y, matG),
    res.one <- plaid.test(X, y, matG, normalize=TRUE, tests="one"),
    res.two <- plaid.test(X, y, matG, normalize=TRUE, tests="two"),
    res.lm <- plaid.test(X, y, matG, normalize=TRUE, tests="lm"),
    res.test <- plaid.test(X, y, matG, normalize=TRUE),
    res.sing <- sing.limma(X, y, gmt)     
  )
  timings  
  res.fgsea <- as.data.frame(res.fgsea)
  rownames(res.fgsea) <- res.fgsea$pathway
  res <- list(
    fgsea = res.fgsea,
    fisher = res.fisher,
    plaid.limma = res.limma,
    plaid.ttest = res.tt,
    plaid.one = res.one,
    plaid.two = res.two,
    plaid.lm = res.lm,
    plaid.test = res.test,
    sing.limma = res.sing    
  )
  res <- res[sapply(res, nrow)>0]
  pp <- Reduce(intersect, lapply(res,rownames))
  score <- cbind(
    fgsea = res.fgsea[pp,"NES"],
    fisher = res.fisher[pp,"sign"],
    plaid.limma = res.limma[pp,"logFC"],
    plaid.ttest = res.tt[pp,"logFC"],
    plaid.one = res.one[pp,"gsetFC"],
    plaid.two = res.two[pp,"gsetFC"],
    plaid.lm = res.lm[pp,"gsetFC"],
    plaid.test = res.test[pp,"gsetFC"],
    sing.limma = res.sing[pp,"logFC"]    
  )  
  pvalue <- cbind(
    fgsea = res.fgsea[pp,"pval"],
    fisher = res.fisher[pp,"p.value"],
    plaid.limma = res.limma[pp,"P.Value"],
    plaid.ttest = res.tt[pp,"P.Value"],
    plaid.one = res.one[pp,"p.one"],
    plaid.two = res.two[pp,"p.two"],
    plaid.lm = res.lm[pp,"p.lm"],
    plaid.test = res.test[pp,"p.meta"],
    sing.limma = res.sing[pp,"P.Value"]    
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

pdf("enrichment-plaid.pdf",h=10,w=10)
pairs( F, pch='.', cex=1, main="enrichment score")
pairs( R, pch='.', cex=1, main="enrichment score (ranks)")
pairs( P, pch='.', cex=1, main="enrichment p-value")
pairs( -log(P), pch='.', cex=1, main="enrichment -log(p)")
dev.off()


##---------------------------------------------------------
##---------------------------------------------------------
##---------------------------------------------------------


