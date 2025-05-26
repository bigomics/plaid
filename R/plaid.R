##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##


#' Compute geneset expression as the average log-ration of genes in
#' the geneset. Requires log-expression matrix X and (sparse) geneset
#' matrix matG.
#' 
#' @param use.cov Logical for whether using a scaled model matrix of the full pairs
#' 
#' @examples# Load the example dataset provided and run:
#' #gmt <- read.gmt("hallmark.gmt")
#' #matG <- gmt2mat(gmt)
#' #gsetX <- plaid(X, matG)
#'
#' @export
plaid <- function(X, matG, stats="mean", chunk=NULL, normalize=TRUE) {
  if (NCOL(X) == 1) X <- cbind(X)
  gg <- intersect(rownames(X), rownames(matG))
  if (length(gg) == 0) {
    message("[plaid] ERROR. no overlapping features")
    return(NULL)
  }
  X <- X[gg, , drop = FALSE]
  matG <- matG[gg, , drop = FALSE]
  G <- 1 * (matG != 0)
  if(stats == "mean") {
    sumG <- 1e-8 + Matrix::colSums(G, na.rm = TRUE)
    G <- Matrix::colScale(G, 1 / sumG)
  }
  ## This single line calcules the plaid score
  gsetX <- chunked_crossprod(G, X, chunk=NULL)
  gsetX <- as.matrix(gsetX)
  if(normalize) {
    gsetX <- normalize_medians(gsetX, ignore.zero=NULL)
  }
  return(gsetX)
}

#' Compute crossprod (t(x) %*% y) for very large y by computing in
#' chunks.
#'
chunked_crossprod <- function(x, y, chunk=NULL) {
  if(is.null(chunk) || chunk < 0) {
    ## if y is large, we need to chunk computation
    Int_max <- .Machine$integer.max
    chunk <- round(0.8 * Int_max / ncol(x))
  }
  if(ncol(y) < chunk) {
    return(Matrix::crossprod(x, y))
  }
  message("[chunked_crossprod] chunked compute: chunk = ", chunk)
  k <- ceiling(ncol(y) / chunk)
  gsetX <- matrix(NA, nrow=ncol(x), ncol=ncol(y),
    dimnames=list(colnames(x),colnames(y)))
  i=1
  for(i in 1:k) {
    jj <- c(((i-1)*chunk+1):min(ncol(y),(i*chunk)))
    xy <- Matrix::crossprod(x, y[,jj])
    gsetX[,jj] <- as.matrix(xy)
  }
  gsetX
}

plaid.limma <- function(X, y, matG) {
  gsetX <- plaid(X, matG)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)   
  res <- playbase::gx.limma(gsetX, y, fdr=1, lfc=0, sort.by='none')
  res
}

#' @export
plaid.ttest <- function(X, y, matG) {
  gsetX <- plaid(X, matG)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)
  res <- Rfast::ttests( Matrix::t(gsetX), ina=y+1)
  m1 <- rowMeans(gsetX[,y==1],na.rm=TRUE)
  m0 <- rowMeans(gsetX[,y==0],na.rm=TRUE)
  padj <- stats::p.adjust(res[,"pvalue"], method="fdr")
  fc <- m1 - m0
  avg <- rowMeans(gsetX)
  df <- data.frame(logFC=fc, AveExpr=avg, t=res[,"stat"],
    P.Value = res[,"pvalue"], adj.P.Val = padj,
    AveExpr.0 = m0, AveExpr.1 = m1)
  df
}

#'
#'
#' @export
plaid.test <- function(X, y, G, gsetX=NULL, normalize=FALSE,
                       tests = c("one","lm") ) {
  if(!all(unique(y) %in% c(0,1))) stop("elements of y must be 0 or 1")
  if(is.list(G)) {
    message("[plaid.test] converting gmt to sparse matrix...")
    G <- gmt2mat(G)
  } else {
    message("[plaid.test] sparse matrix provided")
  }
  p1=p2=p3=NULL
  df1=df2=df3=NULL
  gg <- intersect(rownames(G),rownames(X))
  X <- X[gg,]
  G <- G[gg,]
  
  m1 <- Matrix::rowMeans(X[,y==1,drop=FALSE])
  m0 <- Matrix::rowMeans(X[,y==0,drop=FALSE])
  fc <- m1 - m0
  
  if("one" %in% tests) {
    mt1 <- matrix_ttest(fc, G, method="one")
    p1 <- mt1$pvalue[,1]
    df1 <- mt1$diff[,1]
  }
  if("two" %in% tests) {
    mt2 <- matrix_ttest(fc, G, method="two")
    p2 <- mt2$pvalue[,1]
    df2 <- mt2$diff[,1]        
  }
  if("lm" %in% tests) {
    if(is.null(gsetX)) {      
      message("[plaid.test] computing plaid scores...")
      gsetX <- plaid(X, G)
    }
    if(normalize) gsetX <- normalize_medians(gsetX, ignore.zero=NULL) 
    message("[plaid.test] computing t-tests...")
    res.lm  <- Rfast::ttests( Matrix::t(gsetX), ina=y+1)
    p3 <- res.lm[,"pvalue"]
    df3 <- rowMeans(gsetX[,y==1]) - rowMeans(gsetX[,y==0])
    names(p3) <- rownames(gsetX)
    names(df3) <- rownames(gsetX)
  }
  
  P <- list("one"=p1,"two"=p2,"lm"=p3)
  P <- P[!sapply(P,is.null)]
  F <- list("one"=df1,"two"=df2,"lm"=df3)
  F <- F[!sapply(F,is.null)]
  
  gg <- Reduce(intersect, lapply(P, names))
  P <- lapply(P, function(x) x[gg])
  F <- lapply(F, function(x) x[gg])
  P <- do.call( cbind, P)
  F <- do.call( cbind, F)
  P[is.na(P)] <- 1
  ##F[is.na(F)] <- 0
  gsetFC <- rowMeans(F)
  
  if(NCOL(P)>1) {
    pmeta <- apply(P, 1, function(x) metap::sumz(x)$p)
  } else {
    pmeta <- P[,1]
  }
  colnames(P) <- paste0("p.",colnames(P))
  qmeta <- stats::p.adjust(pmeta, method="fdr")
  res <- cbind(
    gsetFC = gsetFC,
    pvalues = P,
    p.meta = pmeta,
    q.meta = qmeta    
  )
  res <- res[order(res[,"p.meta"]),]
  res
}


#'
#'
#' @export
replaid.scse <- function(X, matG, removeLog2=NULL, scoreMean=FALSE) {
  if(is.null(removeLog2)) {
    removeLog2 <- min(X,na.rm=TRUE)==0 && max(X,na.rm=TRUE) < 20
  }
  if(removeLog2)  {
    message("[replaid.scse] removing log2")
    if(inherits(X,"dgCMatrix")) {
      X@x <- 2**X@x
    } else {
      nz <- Matrix::which(X>0)
      X[nz] <- 2**X[nz]  ## undo only non-zeros as in scSE code
    }
  }
  if(scoreMean) {
    ## modified scSE with Mean-statistics    
    sX <- plaid(X, matG, stats="mean", normalize=FALSE)
    sumx <- Matrix::colMeans(abs(X)) + 1e-8    
    sX <- sX %*% Matrix::Diagonal(x = 1/sumx)
  } else {
    ## original scSE with Sum-statistics
    sX <- plaid(X, matG, stats="sum", normalize=FALSE)
    sumx <- Matrix::colSums(abs(X)) + 1e-8
    sX <- sX %*% Matrix::Diagonal(x = 1/sumx) * 100      
  }
  colnames(sX) <- colnames(X)
  as.matrix(sX)
}

#'
#'
#' @export
replaid.sing <- function(X, matG) {
  X <- methods::as(X, "CsparseMatrix")
  ## the ties.method=min is important for exact replication
  rX <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="min"))  
  rX <- rX / nrow(X) - 0.5
  plaid(rX, matG=matG, normalize=FALSE)
}

##----------------------------------------------------------------
##-------------------- UTILITIES ---------------------------------
##----------------------------------------------------------------

matrix_ttest <- function(F, G, method="two") {
  t_matrix <- matrix(NA, ncol(G), NCOL(F))
  p_matrix <- matrix(NA, ncol(G), NCOL(F))
  f_matrix <- matrix(NA, ncol(G), NCOL(F))  
  dimnames(t_matrix) <- list(colnames(G),colnames(F))
  dimnames(p_matrix) <- list(colnames(G),colnames(F))
  dimnames(f_matrix) <- list(colnames(G),colnames(F))
  i=1
  F2 <- F
  nf <- NCOL(F)
  if(NCOL(F)==1) F2 <- cbind(F,0) ## Rfast::ttest needs matrix
  for(i in 1:ncol(G)) {
    grp <- 2 - 1*(G[,i]!=0)
    if(method=="two") {
      tt <- try(Rfast::ttests(F2, ina=grp),silent=TRUE)
    } else {
      tt <- try(Rfast::ttest(F2[grp==1,], m=rep(0,ncol(F2))),silent=TRUE)
    }
    if(!"try-error" %in% class(tt)) {
      diff <- colMeans(F2[grp==1,,drop=FALSE]) - colMeans(F2[grp==2,,drop=FALSE] )
      t_matrix[i,] <- tt[1:nf,"stat"]
      p_matrix[i,] <- tt[1:nf,"pvalue"]
      f_matrix[i,] <- diff[1:nf]
    }
  }
  list(diff = f_matrix, stats = t_matrix, pvalue=p_matrix)
}

#' @export
normalize_medians <- function(x, ignore.zero=NULL) {
  if(is.null(ignore.zero)) {
    ignore.zero <- (min(x,na.rm=TRUE)==0)
  }
  x <- as.matrix(x)
  if(ignore.zero) {
    zx <- x
    zx[Matrix::which(x==0)] <- NA
    #medx <- Rfast::colMedians(zx, na.rm=TRUE)
    medx <- matrixStats::colMedians(zx, na.rm=TRUE)    
    medx[is.na(medx)] <- 0
  } else {
    #medx <- Rfast::colMedians(x, na.rm=TRUE)
    medx <- matrixStats::colMedians(x, na.rm=TRUE)    
  }
  sweep(x, 2, medx, '-') + mean(medx, na.rm=TRUE)
}

#' @export
colranks <- function(X, sparse=NULL, signed=FALSE, keep.zero=FALSE,
                     ties.method = "average") {
  if(is.null(sparse)) {
    sparse <- inherits(X,"CsparseMatrix")
  }
  if(sparse) {
    X <- methods::as(X, "CsparseMatrix")
    if(keep.zero) {
      rX <- sparse_colranks(X, signed=signed, ties.method=ties.method)
    } else {
      if(signed) {
        sign.X <- sign(X)
        abs.rX <- Matrix::t(sparseMatrixStats::colRanks(abs(X), ties.method=ties.method))
        rX <- abs.rX * sign.X
      } else {
        rX <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method=ties.method))
      }
    }
  } else {
    if(signed) {
      sign.X <- sign(X)
      abs.rX <- Matrix::t(matrixStats::colRanks(as.matrix(abs(X)), ties.method=ties.method))
      rX <- sign.X * abs.rX
    } else {
      rX <- Matrix::t(matrixStats::colRanks(as.matrix(X), ties.method=ties.method))
    }
  }
  rX
}

sparse_colranks <- function(X, signed=FALSE, ties.method="average") {
  ## https://stackoverflow.com/questions/41772943
  X <- methods::as(X, "CsparseMatrix")
  n <- diff(X@p)  ## number of non-zeros per column
  lst <- split(X@x, rep.int(1:ncol(X), n))  ## columns to list
  ## column-wise ranking and result collapsing
  if(signed) {
    lst.sign <- lapply(lst, sign)
    lst.rnk  <- lapply(lst, function(x) rank(abs(x),ties.method=ties.method))
    rnk <- unlist(mapply('*', lst.sign, lst.rnk, SIMPLIFY=FALSE))  
  } else {
    rnk <- unlist(lapply(lst, rank, ties.method=ties.method))  
  }
  rX <- X  ## copy sparse matrix
  rX@x <- rnk  ## replace non-zero elements with rank
  rX
}

##-------------------------------------------------------------
##------------------ end of file ------------------------------
##-------------------------------------------------------------
