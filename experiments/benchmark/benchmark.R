library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

#source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("../R/functions.R")
source("../R/datasets.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations

##----------------------------------------------------
## Computation vs. geneset collection size
##----------------------------------------------------

DATASETS <- c("pbmc3k","tcga-brca_pub")
DATASETS <- c("tcga-brca_pub")
ds=DATASETS[1]
ds=DATASETS[2]

for(ds in DATASETS) {

  cat("*****",ds,"*****\n")
  dataset <- get_dataset(ds, n=Inf)
  str(dataset)
  cat("name = ",dataset$name,"\n")

  X <- dataset$X
  y <- dataset$y
  dim(X)

  matG <- Matrix::t(playdata::GSETxGENE)
  colnames(matG) <- gsub("[/.]","_",colnames(matG))
  dim(matG)

  gg <- intersect(rownames(X),rownames(matG))
  X <- X[gg,]
  matG <- matG[gg,]
  dim(matG)
  gmt <- mat2gmt(matG)

  dim(X)
  XL <- do.call(cbind, rep(list(X), 10))
  XL <- do.call(cbind, rep(list(XL), 4))
  dim(XL)

  xlist <- list( X[,1:100], XL[,1:1000], XL[,1:10000])
  #xlist <- list( X[,1:10], X[,1:100], X[,1:1000])
  #xlist <- list( XL[,1:1000])
  sapply(xlist, ncol)
  remove(XL)

  gmt <- mat2gmt(matG)
  gmt.list <- list(
    hallmark = gmt[grep("HALLMARK",names(gmt))],
    gobp = gmt[grep("GO_BP",names(gmt))],
    gmt = gmt
  )
  sapply(gmt.list, length)

  fn = paste0("benchmark-",ds,"@p14.csv")
  fn
  
  timings <- c()
  ##timings <- read.csv(file=fn,row.names=1)
  m=1;k=1
  for(k in 1:length(gmt.list)) {
    for(m in 1:length(xlist)) {
      gmt1 <- gmt.list[[k]]
      X1 <- xlist[[m]]
      colnames(X1) <- paste(1:ncol(X1),"_",colnames(X1))
      names(gmt1) <- paste0("geneset_",1:length(gmt1))
      message("----------------------")
      message("num.gsets = ", length(gmt1))
      message("num.samples = ", ncol(X1))

      methods=c("plaid","plaid.r","plaid.c","plaid.rc")
      methods=c("cor","rankcor","sing","gsva","ssgsea","ucell","aucell",
                "scse", "scse.mean",
                "replaid.scse","replaid.sing",
                "plaid","plaid.r","plaid.c","plaid.rc")
      tt <- run.timings(X1, gmt1, timeout=(60*60), methods=methods)
      tt <- cbind(tt, nsets=length(gmt1), nrow=nrow(X1), ncol=ncol(X1), dataset=ds)
      timings <- rbind(timings, tt)
      write.csv(timings, file=fn)
    }
  }
  timings
}


if(0) {

  timings <- read.csv(file="benchmark-v2@tokyo.csv", row.names=1)  
  timings <- read.csv(file="benchmark-pbmc3k@p14.csv", row.names=1)
  timings <- read.csv(file="benchmark-tcga-brca_pub@p14.csv", row.names=1)  
  head(timings)
  timings

  sel.methods <- c("sing","gsva","ssgsea", ## "cor",
    "scse", "aucell", "ucell",  #"replaid.scse","replaid.sing",
    "plaid")
  timings <- timings[timings$Function_Call %in% sel.methods,]
  timings
  
  pdf("benchmark-brca@p14.pdf",w=11,h=5.5,pointsize=13)
  nsets <- unique(timings$nsets)
  nsamples <- unique(timings$ncol)
  nc <- length(nsamples)
  nr <- 3
  k=nsets[1];m=nsamples[3]
  par(mfrow=c(nr,2*nc), mar=c(3.8,5,2,2.4), mgp=c(2.1,0.8,0),
      oma=c(0,1,0,1.4))
  for(k in nsets) {
    for(m in nsamples) {
      sel <- which(timings$nsets==k & timings$ncol==m)
      if(length(sel)==0) next
      
      tt <- timings$Elapsed_Time_sec[sel]
      is.timeout <- (timings$Timeout[sel] & tt>3000) | tt>3600
      tt[is.timeout] <- 3600 + tt0[is.timeout]*0.1
      names(tt) <- timings$Function_Call[sel]
      barplot(sort(tt), xlab="Time_sec", las=1, #log=xlog,
              cex.names=1, horiz=TRUE, col="tan2")
      title(paste0("S=",k,"; N=",m),line=0.3, cex.main=1)      
      
      tt1 <- ifelse(tt>10,round(tt,1),round(tt,3))
      lab.tt <- ifelse(is.timeout,paste0(">",tt1),tt1)
      lab.tt <- ifelse(is.timeout,">3600",tt1)        
      text(tt,1.2*(rank(tt)-0.5),labels=lab.tt,pos=4,offset=0.3,
           cex=0.85,xpd=TRUE)        
      
      mb <- timings$Peak_RAM_Used_MiB[sel]
      mb <- ifelse(mb>10,round(mb,1),round(mb,2))
      names(mb) <- timings$Function_Call[sel]      
      barplot(sort(mb), xlab="Peak_RAM_MiB", las=1, #log=xlog,
              cex.names=1, horiz=TRUE, col="dodgerblue3")      
      title(paste0("S=",k,"; N=",m),line=0.3, cex.main=1)      
      lab.mb <- ifelse(is.timeout,paste0(mb,"*"),mb)        
      text(mb,1.2*(rank(mb)-0.5),labels=lab.mb,pos=4,offset=0.3,
           cex=0.85,xpd=TRUE)
      
    }
  }
  dev.off()


  
}


barplot(c(1,10,100),horiz=TRUE,log='x')
