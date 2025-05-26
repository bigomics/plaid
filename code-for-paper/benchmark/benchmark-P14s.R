library(playbase)
library(Seurat)
library(SeuratData)
library(peakRAM)

source("~/Playground/playbase/dev/include.R", chdir=TRUE)
source("functions.R")

#InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final <- UpdateSeuratObject(pbmc3k.final)
Idents(pbmc3k.final) <- pbmc3k.final$seurat_annotations

##----------------------------------------------------
## Computation vs. geneset collection size
##----------------------------------------------------

X <- pbmc3k.final[['RNA']]@data
dim(X)

matG <- playdata::GSETxGENE
rownames(matG) <- gsub("[/.]","_",rownames(matG))
dim(matG)

gg <- intersect(rownames(X),colnames(matG))
X <- X[gg,]
matG <- matG[,gg]
dim(matG)

dim(X)
XL <- do.call(cbind, rep(list(X), 10))
XL <- do.call(cbind, rep(list(XL), 5))
dim(XL)
#xlist <- list( X[,1:10], X[,1:100], X[,1:1000], XL[,1:10000])
#xlist <- list( X[,1:10], X[,1:100], X[,1:1000])
xlist <- list( XL[,1:10000])
sapply(xlist, ncol)
remove(XL)

gmt <- mat2gmt(t(matG))
gmt.list <- list(
  #hallmark = gmt[grep("HALLMARK",names(gmt))],
  #gobp = gmt[head(grep("GO_BP",names(gmt)),1000)],
  gmt = head(gmt,10000)
##  gmt = gmt
)
sapply(gmt.list, length)

timings <- c()
timings <- read.csv(file="benchmark-P14s.csv",row.names=1)

m=2;k=1
for(k in 1:length(gmt.list)) {
  for(m in 1:length(xlist)) {
    gmt1 <- gmt.list[[k]]
    X1 <- xlist[[m]]
    colnames(X1) <- paste(1:ncol(X1),"_",colnames(X1))
    names(gmt1) <- paste0("geneset_",1:length(gmt1))
    message("----------------------")
    message("num.gsets = ", length(gmt1))
    message("num.samples = ", ncol(X1))
    res1 <- run.methods(X1, gmt1)
    tt <- res1$timings
    tt$Total_RAM_Used_MiB <- NULL
    tt <- cbind(tt, nsets=length(gmt1), nrow=nrow(X1), ncol=ncol(X1))
    timings <- rbind(timings, tt)
    write.csv(timings, file="benchmark-P14s.csv")
  }
}
timings


if(0) {

  timings <- read.csv(file="benchmark-P14s.csv", row.names=1)
  timings

  sel.methods <- c("sing","gsva","ssgsea", ## "cor",
    "scse.mean", "aucell", "ucell",
    "replaid.scse","replaid.sing","plaid")
  timings <- timings[timings$Function_Call %in% sel.methods,]

    
  pdf("benchmark-P14s.pdf",w=15,h=5.8,pointsize=13)
  nsets <- unique(timings$nsets)
  nsamples <- unique(timings$ncol)
  nc <- length(nsamples)
  nr <- 3
  par(mfrow=c(nr,2*nc), mar=c(4,5,2,2.5), mgp=c(2.2,1,0),
      oma=c(0,2,0,1))
  for(k in nsets) {
    for(m in nsamples) {
      sel <- which(timings$nsets==k & timings$ncol==m)
      ##sel <- sel[rev(order(timings$Function_Call[sel]))]
      if(length(sel)==0) next
      tt <- timings$Elapsed_Time_sec[sel]
      ##tt <- log10(tt)
      names(tt) <- timings$Function_Call[sel]
      barplot(sort(tt), xlab="Time_sec", las=1,
        cex.names=1, horiz=TRUE, col="tan2")
      title(paste0("S=",k,"; N=",m),line=0.4, cex.main=1)      
      mb <- timings$Peak_RAM_Used_MiB[sel]
      names(mb) <- timings$Function_Call[sel]      
      barplot(sort(mb), xlab="Peak_RAM_MiB", las=1,
        cex.names=1, horiz=TRUE, col="dodgerblue3")      
      title(paste0("S=",k,"; N=",m),line=0.4, cex.main=1)      
    }
  }
  dev.off()
  
}
