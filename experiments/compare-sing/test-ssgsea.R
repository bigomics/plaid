library(playbase)
library(plaid)

#source("~/Playground/playbase/dev/include.R", chdir=TRUE)
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
ds=DATASETS[2]
ds=DATASETS[4]

method="ssgsea"
alpha=0.25

pdf("compare-vs-datasets-SSGEA.pdf", h=16, w=11, pointsize=12)
par(mfrow=c(7,5), mar=c(3.4,4,2.4,1), mgp=c(2.1,0.8,0))
for(alpha in c(0,0.25,1)) {
  cat("============ ",alpha," =============\n")
  
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
    
    method1 = "ssGSEA"
    method2 = "replaid.ssgsea"

    ## run ranked plaid and sing
    time0 <- peakRAM::peakRAM( G1 <- run.ssgsea(X, gmt, alpha=alpha) )
    ##time0 <- peakRAM::peakRAM( G1 <- gao.ssgsea(X, gmt, alpha=alpha) )    
    time1 <- peakRAM::peakRAM( G2 <- plaid::replaid.ssgsea(X, matG, alpha=alpha) )
        
    ss <- intersect(rownames(G1),rownames(G2))
    G1 <- G1[ss,]
    G2 <- G2[ss,]

    #GG <- cbind(GSVA=G1[,1], Gao=G2[,1], replaid=G3[,1])
    #pairs(GG)
    relG1 <- G1 - rowMeans(G1)
    relG2 <- G2 - rowMeans(G2)
  
    G1 <- plaid::normalize_medians(G1)
    G2 <- plaid::normalize_medians(G2)
    m1 <- gx.limma(G1, y, lfc=0, fdr=1, sort.by="none")
    m2 <- gx.limma(G2, y, lfc=0, fdr=1, sort.by="none")
    gg <- intersect(rownames(m1),rownames(m2))
    m1 <- m1[gg,]
    m2 <- m2[gg,]
    
    rr <- c(
      score = cor(G1[,1], G2[,1]),
      rel.score = cor(relG1[,1], relG2[,1]),
      fc = cor(m1$logFC, m2$logFC),
      pv = cor(m1$P.Value, m2$P.Value),
      mlog.pv = cor(-log10(m1$P.Value), -log10(m2$P.Value))
    )
    rr
    rr <- round(rr,3)
    rr.str <- ifelse(rr==1, "(r>0.999)", paste0("(r=",rr,")"))

    par(mar=c(3.4,4,2.4,1))
    plot( G1[,1], G2[,1], xlab=paste(method1,"score"), ylab=method2)
    title(paste0("score  ",rr.str[1]), cex.main=1.2)
    legend("topleft", legend=toupper(ds), cex=0.8, bty="n")
    
    plot( relG1[,1], relG2[,1], xlab=paste("relative",method1,"score"), ylab=paste("relative",method2))
    title(paste0("relative score  ",rr.str[2]), cex.main=1.2)
    
    plot( m1$logFC, m2$logFC, xlab=paste(method1,"score"), ylab=method2)
    title(paste0("logFC  ",rr.str[3]), cex.main=1.2)
    
#    plot( m1$P.Value, m2$P.Value, xlab=paste(method1,"score"), ylab=method2)
#    title(paste0("p-value  ",rr.str[4]), cex.main=1.2)
    
    plot( -log10(m1$P.Value), -log10(m2$P.Value), xlab=paste(method1,"score"),
      ylab=method2)
    title(paste0("-log10p  ",rr.str[5]), cex.main=1.2)

    ## write timings
    t0 <- round(time0[1,"Elapsed_Time_sec"],3)
    m0 <- round(time0[1,"Peak_RAM_Used_MiB"],1)
    t1 <- round(time1[1,"Elapsed_Time_sec"],3)
    m1 <- round(time1[1,"Peak_RAM_Used_MiB"],1)
    plot.new();par(mar=c(1,0,1,1)*2)
    text(0,0.9,paste0(method1,": ",t0,"sec"),adj=0, cex=1)
    text(0,0.8,paste0(method1,": ",m0,"MiB"),adj=0, cex=1)
    text(0,0.6,paste0(method2,": ",t1,"sec"),adj=0, cex=1)
    text(0,0.5,paste0(method2,": ",m1,"MiB"),adj=0, cex=1)    

    text(0,0.3,paste0("nsets: ",ncol(matG)),adj=0, cex=1)
    text(0,0.2,paste0("nsamples: ",ncol(X)),adj=0, cex=1)
    text(0,0.1,paste0("ngenes: ",nrow(X)),adj=0, cex=1)        
  }
}
dev.off()


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------


## run ranked plaid and sing
alpha=0.25
G1 <- run.ssgsea(X, gmt, alpha=alpha)
G2 <- gao.ssgsea(X, gmt, alpha=alpha)    
G3 <- plaid::replaid.ssgsea(X, matG, alpha=alpha) 

ss <- intersect(rownames(G1),rownames(G2))
ss <- intersect(ss,rownames(G3))    
G1 <- G1[ss,]
G2 <- G2[ss,]
G3 <- G3[ss,]    

GG <- cbind(GSVA=G1[,1], Gao=G2[,1], replaid=G3[,1])
pairs(GG)


##-----------------------------------------------------------------
##-----------------------------------------------------------------
##-----------------------------------------------------------------

ds = "GSE10846-dlbcl",
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
G1 <- run.ssgsea(X, gmt, alpha=0)
G2 <- run.ssgsea(X, gmt, alpha=0.25)
G3 <- run.ssgsea(X, gmt, alpha=0.50) 

R1 <- plaid::replaid.ssgsea(X, matG, alpha=0)
R2 <- plaid::replaid.ssgsea(X, matG, alpha=0.25)
R3 <- plaid::replaid.ssgsea(X, matG, alpha=4*0.25) 

ss <- intersect(rownames(G1),rownames(R1))
res.list <- list(ssgsea.a0=G1, ssgsea.a025=G2, ssgsea.a050=G3,
                 replaid.a0=R1, replaid.a025=R2, replaid.a050=R3)
res.list <- lapply(res.list, function(x) x[ss,])
res.list <- lapply(res.list, plaid::normalize_medians)

relX.list <- lapply(res.list, function(x) x - rowMeans(x))

lm.list <- lapply(res.list, function(x) gx.limma(x, y, lfc=0, fdr=1,
                                                 sort.by="none"))
lm.list <- lapply(lm.list, function(m) m[rownames(lm.list[[1]]),])

GG <- sapply(res.list, function(x) x[,1])
RR <- sapply(relX.list, function(x) x[,1])
FF <- sapply(lm.list, function(x) x$logFC)
PP <- sapply(lm.list, function(x) x$P.Value)

pairs(GG, main="score")
pairs(RR, main="relative score")
pairs(FF, main="logFC")
pairs(-log10(PP), main="-log10p")

