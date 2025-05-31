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

sdx <- mat.rowsds(X)
zX <- X
zX <- (X - Matrix::rowMeans(X)) 
zX <- (X - Matrix::rowMeans(X)) / mat.rowsds(X)
zX <- t(apply(X,1,function(x) ecdf(x)(x)))
zX <- zX * sdx

rX <- colranks(zX, signed=TRUE, ties.method="average")
rX <- rX / max(abs(rX))
S2 <- plaid(rX, matG)

S2 <- replaid.gsva(X, matG, tau=0, rowtf='z')
S2 <- S2[rownames(S1),]
plot(S1[,1], S2[,1], xlab="gsva", ylab="replaid")

S3 <- replaid.gsva(X, matG, tau=0, rowtf='ecdf')
S4 <- replaid.ssgsea(X, matG, alpha=0)
S3 <- S3[rownames(S2),]
S4 <- S4[rownames(S2),]
ss <- cbind('gsva'=S1[,1], 'replaid.z'=S2[,1],
            'replaid.ecdf'=S3[,1], 'ssgsea'=S4[,1])
pairs(ss)

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


