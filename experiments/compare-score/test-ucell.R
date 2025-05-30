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

dataset <- get_dataset("testis50")
dataset <- get_dataset("pbmc3k")
X <- dataset$X
y <- dataset$y
dim(X)

gg <- intersect(rownames(X),rownames(full.matG))
X <- X[gg,]
matG <- full.matG[gg,]
matG <- matG[,Matrix::colSums(matG!=0)>10]
dim(matG)
dim(X)
gmt <- mat2gmt(matG)

##-------------------------------------------------------
##-------------------- UCELL ----------------------------
##-------------------------------------------------------

S1 <- t(UCell::ScoreSignatures_UCell(X, gmt))  ## needs logx
rownames(S1) <- sub("_UCell$","",rownames(S1))

rX <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="average"))
rX <- max(rX) - rX
rmax = 1500
rX <- pmin(rX, rmax+1)
S2 <- plaid(rX, matG)
S2 <- 1 - S2 / rmax + (colSums(matG!=0)+1)/(2*rmax)

##S2 <- replaid.ucell(X, matG, rmax=1500)
S2 <- S2[rownames(S1),]
plot(S1[,1], S2[,1], xlab="UCell", ylab="replaid")

##-------------------------------------------------------
##-------------------- AUCELL ---------------------------
##-------------------------------------------------------

aucMaxRank = ceiling(0.05 * nrow(X))
S1 <- AUCell::getAUC(AUCell::AUCell_run(X, gmt, aucMaxRank=aucMaxRank))

par(mfrow=c(2,2))
rX <- colranks(X, ties.method="average")
wX <- 1.08*pmax((rX - (max(rX) - aucMaxRank)) / aucMaxRank, 0)
S2 <- plaid(wX, matG, stats="mean")
##S2 <- plaid(X, matG)
plot(S1[,1], S2[,1], xlab="AUCell", ylab="plaid")


N=10
S2 <- matrix(0, ncol(matG), ncol(aX))
i=1
for(i in 1:N) {
  aX <- 1 * (rX > (max(rX) - i/N*aucMaxRank))
  dS <- plaid(aX, matG, stats="mean")
  S2 <- S2 + dS/N
}
plot(S1[,1], S2[,1], xlab="AUCell", ylab="plaid")

S3 <- replaid.aucell(X, matG)
plot(S1[,1], S3[,1], xlab="AUCell", ylab="replaid")

