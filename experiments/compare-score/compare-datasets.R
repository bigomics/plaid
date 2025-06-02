library(plaid)

##source("~/Playground/playbase/dev/include.R", chdir=TRUE)
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

all.methods <- c("scse","scse.mean","sing","ssgsea","gsva","aucell","ucell")
method="sing"
method="ssgsea"

pdf("compare-vs-datasets.pdf", h=16, w=11, pointsize=12)
for(method in all.methods) {
  par(mfrow=c(7,5), mar=c(3.4,4,2.4,1), mgp=c(2.1,0.8,0))
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
    time0=time1=NULL
    if(method == "scse.mean") {
      prepare.SCSE(X, gmt, path="../scse")
      time0 <- peakRAM::peakRAM(G1 <- run.SCSE(X, gmt, removeLog2=FALSE, scoreMean=TRUE, path="../scse"))
      time1 <- peakRAM::peakRAM(G2 <- replaid.scse(X, matG, removeLog2=FALSE, scoreMean=TRUE))
      method1 = "scse.mean"
      method2 = "replaid.scse"
    }
    if(method == "scse") {
      prepare.SCSE(X, gmt, path="../scse")
      time0 <- peakRAM::peakRAM(G1 <- run.SCSE(X, gmt, removeLog2=1, scoreMean=FALSE, path="../scse"))
      time1 <- peakRAM::peakRAM(G2 <- replaid.scse(X, matG, removeLog2=1, scoreMean=FALSE))
      method1 = "scse"
      method2 = "replaid.scse"
    }
    if(method == "sing") {
      time0 <- peakRAM::peakRAM( G1 <- gset.singscore(X, gmt, return.score = TRUE))
      time1 <- peakRAM::peakRAM( G2 <- replaid.sing(X, matG) )
      method1 = "sing"
      method2 = "replaid.sing"
    }
    if(method == "ssgsea") {
      ##G1 <- gset.gsva(X, gmt, method="ssgsea")
      time0 <- peakRAM::peakRAM( G1 <- run.ssgsea(X, gmt, alpha=0) )
      time1 <- peakRAM::peakRAM( rX <- colranks(X, keep.zero=TRUE, ties.method="average") )
      G2 <- plaid::plaid(rX, matG)
      method1 = "ssGSEA"
      method2 = "replaid.ssgsea"
    }
    if(method == "gsva") {
      ##G1 <- gset.gsva(X, gmt, method="ssgsea")
      time0 <- peakRAM::peakRAM( G1 <- run.gsva(X, gmt, tau=0) )
      time1 <- peakRAM::peakRAM( G2 <- replaid.gsva(X, matG))
      method1 = "GSVA"
      method2 = "replaid.gsva"
    }
    if(method == "ucell") {
      time0 <- peakRAM::peakRAM(
        G1 <- t(UCell::ScoreSignatures_UCell(X, gmt))  ## needs logx
      )
      rownames(G1) <- sub("_UCell$","",rownames(G1))
      time1 <- peakRAM::peakRAM( G2 <- replaid.ucell(X, matG))
      method1 = "UCell"
      method2 = "replaid.ucell"
    }
    if(method == "aucell") {
      time0 <- peakRAM::peakRAM(G1 <- AUCell::getAUC(AUCell::AUCell_run(X, gmt)))
      time1 <- peakRAM::peakRAM(G2 <- replaid.aucell(X, matG) )
      method1 = "AUCell"
      method2 = "replaid.aucell"
    }
    
    gg <- intersect(rownames(G1),rownames(G2))
    G1 <- G1[gg,]
    G2 <- G2[gg,]
    ##plot(G1[,1], G2[,1])
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
