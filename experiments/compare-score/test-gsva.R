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
dataset <- get_dataset("pbmc3k",n=40)
dataset <- get_dataset("geiger")
X <- dataset$X
y <- dataset$y
table(y)
dim(X)
X <- X[mat.rowsds(X)>0.1,]

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

par(mfrow=c(3,3),mar=c(4,4,2,2))

gsvapar <- GSVA::gsvaParam(X, gmt, tau=0, kcdf='Gaussian', maxDiff=TRUE)
S1 <- GSVA::gsva(gsvapar, verbose = FALSE)
gsvapar <- GSVA::gsvaParam(X, gmt, tau=1, kcdf='Gaussian', maxDiff=TRUE)
S2 <- GSVA::gsva(gsvapar, verbose = FALSE)

S3 <- replaid.gsva(X, matG, tau=0, rowtf='z')
S4 <- replaid.gsva(X, matG, tau=0, rowtf='ecdf')
S5 <- replaid.ssgsea(X, matG, alpha=0)
S6 <- gset.gsva(X, gmt, method="ssgsea")
S7 <- gset.gsva(X, gmt, method="zscore")

S2 <- S2[rownames(S1),]
S3 <- S3[rownames(S1),]
S4 <- S4[rownames(S1),]
S5 <- S5[rownames(S1),]
S6 <- S6[rownames(S1),]
S7 <- S7[rownames(S1),]

## try out
sdx <- mat.rowsds(X)
zX <- X
zX <- (X - Matrix::rowMeans(X)) 
zX <- (X - Matrix::rowMeans(X)) / sdx

zX <- t(apply(X,1,function(x) ecdf(x)(x)))
Sx <- plaid(zX, matG)

rX <- colranks(zX, signed=TRUE, ties.method="average")
rX <- rX / max(abs(rX))
Sx <- plaid(rX, matG)

ss <- cbind(
  'gsva.tau0'=S1[,1],
  'gsva.tau1'=S2[,1],  
  'replaid.z'=S3[,1],  
  'replaid.ecdf'=S4[,1],
  'ssgsea'=S5[,1],
  'gsva:ssgsea'=S6[,1],
  'sgva:zscore'=S7[,1],
  'test'=Sx[,1]
)
pairs(ss, cex=0.4)

ff <- sapply(list(S1,S2,S3,S4),function(x) (x-rowMeans(x))[,1])
pairs(ff)



par(mfrow=c(2,2))
ss <- c()
for(tau in c(0, 0.5, 1)) {
  gsvapar <- GSVA::gsvaParam(X, gmt, tau=tau, kcdf='Gaussian', maxDiff=TRUE)
  S1 <- GSVA::gsva(gsvapar, verbose = FALSE)  
  S2 <- replaid.gsva(X, matG, tau=tau)
  ss <- cbind(ss, S1[,1], S2[,1])
}
pairs(ss)


