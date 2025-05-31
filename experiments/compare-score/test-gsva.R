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
dataset <- get_dataset("geiger")
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
##-------------------- GSVA -----------------------------
##-------------------------------------------------------
source("../R/functions.R")
source("../../R/plaid.R")

par(mfrow=c(2,2))

gsvapar <- GSVA::gsvaParam(X, gmt, tau=0, kcdf='Gaussian', maxDiff=TRUE)
S1 <- GSVA::gsva(gsvapar, verbose = FALSE)

zX <- (X - rowMeans(X)) / matrixStats::rowSds(X)
rX <- colranks(zX, signed=TRUE, ties.method="average")
rX <- rX / max(abs(rX))
S2 <- plaid(rX, matG)

S2 <- replaid.gsva(X, matG, tau=0)

##S2 <- replaid.ucell(X, matG, rmax=1500)
S2 <- S2[rownames(S1),]
plot(S1[,1], S2[,1], xlab="gsva", ylab="replaid")

par(mfrow=c(2,2))
ss <- c()
for(tau in c(0, 0.5, 1)) {
  gsvapar <- GSVA::gsvaParam(X, gmt, tau=tau, kcdf='Gaussian', maxDiff=TRUE)
  S1 <- GSVA::gsva(gsvapar, verbose = FALSE)  
  S2 <- replaid.gsva(X, matG, tau=tau)
  ss <- cbind(ss, S1[,1], S2[,1])
}
pairs(ss)


