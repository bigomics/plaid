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
X <- X[gg,1:10]
matG <- matG[,gg]
dim(matG)

gmt <- mat2gmt(t(matG))
gmt.list <- list(
  hallmark = gmt[grep("HALLMARK",names(gmt))],
  gobp = gmt[grep("GO_BP",names(gmt))],
  gmt = gmt
)
sapply(gmt.list, length)

##timings <- c()
timings <- read.csv(file="timings-methods2-laptop.csv", row.names=1)
m=1;k=1
for(gmt1 in gmt.list) {
  for(m in c(10,100)) {
    X1 <- do.call(cbind, rep(list(X),m))
    colnames(X1) <- paste(1:ncol(X1),"_",colnames(X1))
    names(gmt1) <- paste0("geneset_",1:length(gmt1))
    res1 <- run.methods(X1, gmt1)
    tt <- res1$timings
    tt <- cbind(tt, nsets=length(gmt1), nrow=nrow(X1), ncol=ncol(X1))
    timings <- rbind(timings, tt)
    write.csv(timings, file="timings-methods2-laptop.csv")
  }
}
timings


if(0) {

  timings <- read.csv(file="timings-methods2-laptop.csv", row.names=1)
  timings
  timings <- timings[timings$Function_Call!="cor",]
  timings$Function_Call <- sub("sing","singscore",timings$Function_Call)
  timings$Function_Call <- sub("ssplaid","plaid",timings$Function_Call)  

  grp <- apply(timings[,5:7],1,paste,collapse="x")
  tapply(timings$Elapsed_Time_sec, grp, function(x) sort(x/min(x))[-1])
  tapply(timings$Peak_RAM_Used_MiB, grp, function(x) sort(x/min(x))[-1])  
  
  
  pdf("timings-methods2-laptop.pdf",w=12,h=5.5,pointsize=12)
  nsets <- unique(timings$nsets)
  nsamples <- unique(timings$ncol)
  par(mfrow=c(3,6), mar=c(4,6.5,2,1), mgp=c(2.3,0.9,0))
  for(k in nsets) {
    for(m in nsamples) {
      sel <- which(timings$nsets==k & timings$ncol==m)
      ##sel <- sel[rev(order(timings$Function_Call[sel]))]
      if(length(sel)==0) next
      tt <- timings$Elapsed_Time_sec[sel]
      mb <- timings$Peak_RAM_Used_MiB[sel]
      names(tt) <- timings$Function_Call[sel]
      names(mb) <- timings$Function_Call[sel]      
      barplot(sort(tt), xlab="Elapsed_Time_sec", las=1,
              cex.names=1.3, horiz=TRUE, col="tan2")
      title(paste0(k," x ",m),line=0.5)      
      barplot(sort(mb), xlab="Peak_RAM_Used_MiB", las=1,
              cex.names=1.3, horiz=TRUE, col="dodgerblue3")      
      title(paste0(k," x ",m),line=0.5)      
    }
  }
  dev.off()

  
  
}
