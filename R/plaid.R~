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
#' gmt <- read.gmt("hallmark.gmt")
#' matG <- gmt2mat(gmt)
#' gsetX <- plaid(X, matG)
#'
#' @export
plaid <- function(X, matG, stats="mean", chunk=NULL, as_matrix=TRUE) {
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
  #gsetX <- Matrix::crossprod(G, X) 
  gsetX <- chunked_crossprod(G, X, chunk=NULL)
  if(as_matrix) gsetX <- as.matrix(gsetX)
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
  gsetX <- plaid(X, matG, as_matrix=TRUE)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)   
  gx.limma(gsetX, y, fdr=1, lfc=0, sort.by='none')
}

#' @export
plaid.ttest <- function(X, y, matG) {
  gsetX <- plaid(X, matG, as_matrix=TRUE)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)
  Rfast::ttests(gsetX, ina=y)
}

#'
#'
#' @export
plaid.test <- function(X, y, gmt, gsetX=NULL, normalize=FALSE,
                       tests = c("one","two","lm") ) {
  if(!all(unique(y) %in% c(0,1))) stop("elements of y must be 0 or 1")
  if(class(gmt)=="list") gmt <- gmt2mat(gmt)
  p1=p2=p3=NULL
  diff=NULL
  if(any( c("one","two") %in% tests)) {
    m1 <- Matrix::rowMeans(X[,y==1,drop=FALSE])
    m0 <- Matrix::rowMeans(X[,y==0,drop=FALSE])
    fc <- m1 - m0
  }
  if("one" %in% tests) {
    mt1 <- matrix_ttest(fc, gmt, method="one")
    p1 <- mt1$pvalue[,1]
    diff <- mt1$diff[,1]
  }
  if("two" %in% tests) {
    mt2 <- matrix_ttest(fc, gmt, method="two")
    p2 <- mt2$pvalue[,1]
    diff <- mt2$diff[,1]        
  }
  if("lm" %in% tests) {
    if(is.null(gsetX)) {
      gsetX <- plaid(X, gmt, as_matrix=TRUE)
    }
    if(normalize) gsetX <- normalize_medians(gsetX, ignore.zero=NULL) 
    lm <- gx.limma(gsetX, y, fdr=1, lfc=0, sort.by='none')
#    lm <- try(Rfast::ttests(gsetX, ina=y),silent=TRUE)    
    p3 <- lm[,"P.Value"]
    diff <- lm[,"logFC"]
    names(p3) <- rownames(lm)
    names(diff) <- rownames(lm)
  }
  pvalues <- list("one"=p1,"two"=p2,"lm"=p3)
  pvalues <- pvalues[!sapply(pvalues,is.null)]
  gg <- Reduce(intersect, lapply(pvalues, names))
  diff <- diff[gg]
  pvalues <- lapply(pvalues, function(x) x[gg])
  pvalues <- do.call( cbind, pvalues)
  pvalues[is.na(pvalues)] <- 1
  if(NCOL(pvalues)>1) {
    pmeta <- apply(pvalues, 1, function(x) metap::sumz(x)$p)
  } else {
    pmeta <- pvalues[,1]
  }
  colnames(pvalues) <- paste0("p.",colnames(pvalues))
  res <- cbind(
    gsetFC = diff,
    pvalues,
    p.meta = pmeta
  )
  res <- res[order(res[,"p.meta"]),]
  res
}

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
    sX <- plaid(X, matG, stats="mean", as_matrix=FALSE)
    sumx <- Matrix::colMeans(abs(X)) + 1e-8    
    sX <- sX %*% Matrix::Diagonal(x = 1/sumx)
  } else {
    ## original scSE with Sum-statistics
    sX <- plaid(X, matG, stats="sum", as_matrix=FALSE)
    sumx <- Matrix::colSums(abs(X)) + 1e-8
    sX <- sX %*% Matrix::Diagonal(x = 1/sumx) * 100      
  }
  colnames(sX) <- colnames(X)
  as.matrix(sX)
}

replaid.sing <- function(X, matG) {
  X <- as(X, "CsparseMatrix")
  ## the ties.method=min is important for exact replication
  rX <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="min"))  
  rX <- rX / nrow(X) - 0.5
  plaid(rX, matG=matG)
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

normalize_medians <- function(x, ignore.zero=NULL) {
  if(is.null(ignore.zero)) {
    ignore.zero <- (min(x,na.rm=TRUE)==0)
  }
  if(ignore.zero) {
    medx <- apply(x,2,function(y) median(y[y!=0],na.rm=TRUE))
    med0 <- median(x[x!=0],na.rm=TRUE)
  } else {
    medx <- apply(x,2,median,na.rm=TRUE)
    med0 <- median(x,na.rm=TRUE)
  }  
  t(t(x) - medx) + med0
}

normalize_means <- function(x, ignore.zero=NULL) {
  if(is.null(ignore.zero)) {
    ignore.zero <- (min(x,na.rm=TRUE)==0)
  }
  if(ignore.zero) {
    medx <- apply(x,2,function(y) mean(y[y!=0],na.rm=TRUE))
    med0 <- mean(x[x!=0],na.rm=TRUE)
  } else {
    medx <- apply(x,2,mean,na.rm=TRUE)
    med0 <- mean(x,na.rm=TRUE)
  }  
  t(t(x) - medx) + med0
}

colranks <- function(X, sparse=NULL, signed=FALSE, keep.zero=FALSE,
                     ties.method = "average") {
  if(is.null(sparse)) {
    sparse <- inherits(X,"CsparseMatrix")
  }
  if(sparse) {
    X <- as(X, "CsparseMatrix")
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

#' Test if rows or columns of sparse matrix are duplicated
sparse_colranks <- function(X, signed=FALSE, ties.method="average") {
  ## https://stackoverflow.com/questions/41772943
  X <- as(X, "CsparseMatrix")
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

#' Convert GMT file to binary matrix
#'
#' This function converts a GMT file (Gene Matrix Transposed) to a binary matrix,
#' where rows represent genes and columns represent gene sets. The binary matrix
#' indicates the presence or absence of genes in each gene set.
#'
#' @param gmt A list representing the GMT file, where each element is a character
#'            vector representing a gene set.
#' @param max.genes The maximum number of genes to include in the binary matrix.
#'                  Defaults to -1, which includes all genes.
#' @param ntop The number of top genes to consider for each gene set. Defaults to -1,
#'             which includes all genes.
#' @param sparse A logical value indicating whether to create a sparse matrix
#'               (from the 'Matrix' package) or a dense matrix. Defaults to `TRUE`.
#' @param bg A character vector representing the background set of genes. Defaults to `NULL`,
#'           which considers all unique genes from the gene sets.
#' @param use.multicore A logical value indicating whether to use parallel processing
#'                     (via the 'parallel' package) for faster execution. Defaults to `TRUE`.
#'
#' @export
#'
#' @return A binary matrix representing the presence or absence of genes in each gene set.
#'         Rows correspond to genes, and columns correspond to gene sets.
#'
gmt2mat <- function(gmt, max.genes = -1, ntop = -1, sparse = TRUE,
                    bg = NULL, use.multicore = TRUE) {
  gmt <- gmt[order(-sapply(gmt, length))]
  gmt <- gmt[!duplicated(names(gmt))]
  if (ntop > 0) {
    gmt <- lapply(gmt, utils::head, n = ntop)
  }
  if (is.null(names(gmt))) names(gmt) <- paste("gmt.", 1:length(gmt), sep = "")
  if (is.null(bg)) {
    bg <- names(sort(table(unlist(gmt)), decreasing = TRUE))
  }
  if (max.genes < 0) max.genes <- length(bg)
  gg <- bg
  gg <- Matrix::head(bg, n = max.genes)
  gmt <- lapply(gmt, function(s) intersect(gg, s))
  kk <- unique(names(gmt))
  if (sparse) {
    D <- Matrix::Matrix(0, nrow = length(gg), ncol = length(kk), sparse = TRUE)
  } else {
    D <- matrix(0, nrow = length(gg), ncol = length(kk))
  }

  rownames(D) <- gg
  colnames(D) <- kk
  j <- 1
  if (use.multicore) {
    idx <- lapply(gmt, function(s) match(s, gg))
    idx[sapply(idx, length) == 0] <- 0
    idx <- sapply(1:length(idx), function(i) rbind(idx[[i]], i))
    idx <- matrix(unlist(idx[]), byrow = TRUE, ncol = 2)
    idx <- idx[!is.na(idx[, 1]), ]
    idx <- idx[idx[, 1] > 0, ]
    D[idx] <- 1
  } else {
    for (j in 1:ncol(D)) {
      k0 <- match(kk[j], names(gmt))
      ii0 <- which(gg %in% gmt[[k0]])
      if (length(ii0) > 0) D[ii0, j] <- +1
    }
  }
  D <- D[order(-Matrix::rowSums(D != 0, na.rm = TRUE)), ]
  D
}

#' Convert binary matrix back to GMT list
#'
#' This function converts binary matrix to a GMT (Gene Matrix
#' Transposed) list, The binary matrix indicates the presence or
#' absence of genes in each gene set, where rows represent genes and
#' columns represent gene sets.
#'
#' @param mat A matrix with non-zero entries representing genes in
#'   each gene set where rows represent genes and columns represent
#'   gene sets.
#'
#' @export
#'
#' @return A list of vector representing each gene set. Each list
#'   element correspond to a gene set and is a vector of genes
#'
mat2gmt <- function(mat) {
  idx <- Matrix::which(mat != 0, arr.ind = TRUE)
  gmt <- tapply(rownames(idx), idx[, 2], list)
  names(gmt) <- colnames(mat)[as.integer(names(gmt))]
  gmt
}

##-------------------------------------------------------------
##------------------ end of file ------------------------------
##-------------------------------------------------------------
