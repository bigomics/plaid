##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Compute PLAID single-sample enrichment score 
#'
#' @description Compute single-sample geneset expression as the
#'   average log-expression f genes in the geneset. Requires log-expression
#'   matrix X and (sparse) geneset matrix matG. If you have gene sets
#'   as a gmt list, please convert it first using the function `gmt2mat()`.
#'
#' @details PLAID needs the gene sets as sparse matrix. If you have
#'   your collection of gene sets a a list, we need first to convert
#'   the gmt list to matrix format.
#' 
#' @details We recommend to run PLAID on the log transformed expression matrix,
#' not on the counts, as the average in the logarithmic space is more
#' robust and is in concordance to calculating the geometric mean.
#'
#' @details It is not necessary to normalize your expression matrix before
#' running PLAID because PLAID performs median normalization of the
#' enrichment scores afterwards.
#'
#' @details It is recommended to use sparse matrix as PLAID relies on
#' sparse matrix computations. But, PLAID is also fast for dense matrices.
#'
#' @details PLAID can also be run on the ranked matrix. This corresponds to
#' the singscore (Fouratan et al., 2018). PLAID can also be run on
#' the (non-logarithmic) counts which can be used to calculate the
#' scSE score (Pont et al., 2019).
#'
#' @details PLAID is fast and memery efficient because it uses efficient
#' sparse matrix computation. When input matrix is very large, PLAID
#' performs 'chunked' computation by splitting the matrix in chunks.
#'
#' @details Although `X` and `matG` are generally sparse, the result
#' matrix `gsetX` generally is dense and can thus be very large.
#' Example: computing gene set scores for 10K gene sets on 1M cells
#' will create a 10K x 1M dense matrix which requires ~75GB memory.
#' 
#' @param X Log-transformed expr. matrix. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#' @param stats Score computation stats: mean or sum of intensity. Default 'mean'.
#' @param chunk Logical: use chunks for large matrices. Default 'NULL' for autodetect.
#' @param normalize Logical: median normalize results or not. Default 'TRUE'.
#'
#' @return Matrix of single-sample enrichment scores.
#' Gene sets on rows, samples on columns.
#' 
#' @examples
#' library(plaid)
#' load(system.file("extdata", "pbmc3k-50cells.rda", package = "PLAID"))
#' hallmarks <- system.file("extdata", "hallmarks.gmt", package = "PLAID")
#' gmt <- read.gmt(hallmarks)
#' matG <- gmt2mat(gmt)
#' gsetX <- plaid(X, matG)
#'
#' @export
plaid <- function(X, matG, stats=c("mean","sum"), chunk=NULL, normalize=TRUE) {

  stats <- stats[1]
  if (NCOL(X) == 1) X <- cbind(X)

  gg <- intersect(rownames(X), rownames(matG))
  if (length(gg) == 0) {
    message("[plaid] ERROR. No overlapping features.")
    return(NULL)
  }

  X <- X[gg, , drop = FALSE]
  matG <- matG[gg, , drop = FALSE]
  G <- 1 * (matG != 0)
  if(stats == "mean") {
    sumG <- 1e-8 + Matrix::colSums(G, na.rm = TRUE)
    G <- Matrix::colScale(G, 1 / sumG)
  }

  ## Calculates PLAID score
  gsetX <- chunked_crossprod(G, X, chunk=NULL)
  gsetX <- as.matrix(gsetX)
  
  if(normalize) gsetX <- normalize_medians(gsetX)
 
  return(gsetX)

}

#' Chunked computation of cross product
#'
#' Compute crossprod (t(x) %*% y) for very large y by computing in
#' chunks.
#'
#' @param x Matrix First matrix for multiplication. Can be sparse.
#' @param y Matrix Second matrix for multiplication. Can be sparse.
#' @param chunk Integer Chunk size (max number of columns) for computation.
#'
#' @return Matrix. Result of matrix cross product.
#' 
chunked_crossprod <- function(x, y, chunk=NULL) {
  if(is.null(chunk) || chunk < 0) {
    ## if y is large, we need to chunk computation
    Int_max <- .Machine$integer.max
    chunk <- round(0.8 * Int_max / ncol(x))
  }

  if(ncol(y) < chunk) return(Matrix::crossprod(x, y))
 
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

  return(gsetX)

}

#' Fast calculation of scSE score
#'
#' @description Calculates Single-Cell Signature Explorer (Pont et
#'   al., 2019) scores using plaid back-end. The computation is
#'   10-100x faster than the original code.
#'
#' @details Computing the scSE requires running plaid on the linear
#'   (not logarithmic) score and perform additional normalization by
#'   the total UMI per sample. We have wrapped this in a single
#'   convenience function:
#'
#' To replicate the original "sum-of-UMI" scSE score, set `removeLog2=TRUE`
#' and `scoreMean=FALSE`. scSE and plaid scores become more similar for
#' `removeLog2=FALSE` and `scoreMean=TRUE`.
#'
#' We have extensively compared the results from `replaid.scse` and
#' from the original scSE (implemented in GO lang) and we showed
#' almost identical results in the score, logFC and p-values.
#' 
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns.
#' @param removeLog2 Logical for whether to remove the Log2, i.e. will
#'   apply power transform (base2) on input (default TRUE).
#' @param scoreMean Logical for whether computing sum or mean as score
#'   (default FALSE).
#'
#' @export
replaid.scse <- function(X,
                         matG,
                         removeLog2 = NULL,
                         scoreMean = FALSE) {

  if(is.null(removeLog2))
    removeLog2 <- min(X, na.rm = TRUE)==0 && max(X, na.rm = TRUE) < 20
  
  if(removeLog2)  {
    message("[replaid.scse] Converting data to linear scale (removing log2)...")
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
  sX <- as.matrix(sX)

  return(sX)

}


#' Fast calculation of singscore
#'
#' @description Calculates single-sample enrichment singscore
#'   (Fouratan et al., 2018) using plaid back-end. The computation is
#'   10-100x faster than the original code.
#'
#' @details Computing the singscore requires to compute the ranks of
#'   the expression matrix. We have wrapped this in a single
#'   convenience function.
#'
#' We have extensively compared the results of `replaid.sing` and from
#' the original `singscore` R package and we showed identical result
#' in the score, logFC and p-values.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns.
#' 
#' @export
replaid.sing <- function(X, matG) {
  ## the ties.method=min is important for exact replication
  rX <- colranks(X, ties.method = "min")
  rX <- rX / nrow(X) - 0.5
  gsetX <- plaid(rX, matG = matG, normalize = FALSE)
  return(gsetX)
}

#' Fast calculation of ssGSEA
#'
#' @description Calculates single-sample enrichment singscore (Barbie
#'   et al., 2009; Hänzelmann et al., 2013) using plaid back-end. The
#'   computation is 10-100x faster than the original code.
#'
#' @details Computing ssGSEA score requires to compute the ranks of
#'   the expression matrix and weighting of the ranks. We have wrapped
#'   this in a single convenience function.
#'
#' We have extensively compared the results of `replaid.ssgsea` and
#' from the original `GSVA` R package and we showed highly similar
#' results in the score, logFC and p-values. For alpha=0 we obtain
#' exact results, for alpha>0 the results are highly similar but not
#' exactly the same.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns.
#' @param alpha Weighting factor for exponential weighting of ranks
#' 
#' @export
replaid.ssgsea <- function(X, matG, alpha = 0) {
  rX <- colranks(X, keep.zero = TRUE, ties.method = "average")
  if(alpha != 0) {
    ## This is not exactly like original formula. Not sure how to
    ## efficiently implement original rank weighting
    rX <- rX^(1 + alpha)
  }
  rX <- rX / max(rX) - 0.5
  dimnames(rX) <- dimnames(X)
  gsetX <- plaid(rX, matG, stats = "mean", normalize = TRUE)
  return(gsetX)
}

#' Fast calculation of UCell
#'
#' @description Calculates single-sample enrichment UCell (Andreatta
#'   et al., 2021) using plaid back-end. The computation is
#'   10-100x faster than the original code.
#'
#' @details Computing ssGSEA score requires to compute the ranks of
#'   the expression matrix and truncation of the ranks. We have wrapped
#'   this in a single convenience function.
#'
#' We have extensively compared the results of `replaid.ucell` and
#' from the original `UCell` R package and we showed near exacct
#' results in the score, logFC and p-values. 
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#' @param rmax Rank threshold (see Ucell paper). Default rmax = 1500.
#' @export
replaid.ucell <- function(X, matG, rmax = 1500) {
  rX <- colranks(X, ties.method = "average")
  rX <- pmin( max(rX) - rX, rmax+1 )
  S <- plaid(rX, matG)
  S <- 1 - S / rmax + (Matrix::colSums(matG!=0)+1)/(2*rmax)
  return(S)
}

#' Fast calculation of AUCell
#'
#' @description Calculates single-sample enrichment AUCell (Aibar
#'   et al., 2017) using plaid back-end. The computation is
#'   10-100x faster than the original code.
#'
#' @details Computing the AUCell score requires to compute the ranks
#'   of the expression matrix and approximating the AUC of a gene
#'   set. We have wrapped this in a single convenience function.
#'
#' We have extensively compared the results of `replaid.aucell` and
#' from the original `UCell` R package and we showed good concordance
#' of results in the score, logFC and p-values.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#' @param aucMaxRank Rank threshold (see AUCell paper). Default aucMaxRank = 0.05*nrow(X).
#' 
#' @export
replaid.aucell <- function(X, matG, aucMaxRank = ceiling(0.05*nrow(X))) {
  rX <- colranks(X, ties.method = "average")
  ww <- 1.08*pmax((rX - (max(rX) - aucMaxRank)) / aucMaxRank, 0)
  gsetX <- plaid(ww, matG, stats = "mean")
  return(gsetX)
}

#' Fast approximation of GSVA
#'
#' @description Calculates single-sample enrichment GSVA (Hänzelmann
#'   et al., 2013) using plaid back-end. The computation is
#'   10-100x faster than the original code.
#'
#' @details Computing the GSVA score requires to compute the CDF of
#'   the expression matrix, ranking and scoring the genesets. We have
#'   wrapped this in a single convenience function.
#'
#' We have extensively compared the results of `replaid.gsva` and
#' from the original `GSVA` R package and we showed good concordance
#' of results in the score, logFC and p-values.
#'
#' In the original formulation, GSVA uses an emperical CDF to
#' transform expression of each feature to a [0;1] relative expression
#' value. For efficiency reasons, this is here approximated by a
#' z-transform (center+scale) of each row.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns.
#' @param tau Rank weight parameter (see GSVA publication). Default
#'   tau=0.
#' 
#' @export
replaid.gsva <- function(X, matG, tau = 0, rowtf = c("z", "ecdf")[1]) {
  rowtf <- rowtf[1]

  if(rowtf == "z") {
    ## Faster approximation of relative activation
    zX <- (X - Matrix::rowMeans(X)) / (1e-8 + mat.rowsds(X))
  } else if(rowtf=='ecdf') {
    ## this implements original ECDF idea
    zX <- t(apply(X,1,function(x) ecdf(x)(x))) 
  } else {
    stop("Error: unknown row transform",rowtf)
  }

  rX <- colranks(zX, signed = TRUE, ties.method = "average")
  rX <- rX / max(abs(rX))
  if(tau > 0) {
    ## Note: This is not exactly like original formula. Not sure how
    ## to efficiently implement original rank weighting
    rX <- sign(rX) * abs(rX)^(1 + tau)
  }
  dimnames(rX) <- dimnames(X)
  gsetX <- plaid(rX, matG)

  return(gsetX)

}

mat.rowsds <- function(X) {
  if(inherits(X,"CsparseMatrix"))
    return(sparseMatrixStats::rowSds(X))
  sdx <- matrixStats::rowSds(X)
  return(sdx)
}

##----------------------------------------------------------------
##-------------------- STAT TEST ---------------------------------
##----------------------------------------------------------------

#' Statistical testing of differentially enrichment
#'
#' This function performs statistical testing for differential
#' enrichment using plaid
#'
#' @param X Matrix of log expression value
#' @param y Vector of 0s and 1s indicating group
#' @param G Sparse matrix of gene sets. Non-zero entry indicates
#'   gene/feature is part of gene sets. Features on rows, gene sets on
#'   columns.
#' @param gsetX Gene set score matrix which is output of
#'   `plaid()`. Can be NULL in that case it will be recomputed from X
#'   and G (default required).
#' @param tests Character array indicating which tests to perform.
#' 
#' @export
plaid.test <- function(X, y, G, gsetX, tests = c("one","two","lm"),
                       metap.method="fisher", sort.by='p.meta') {
  if(!all(unique(y) %in% c(0,1))) stop("elements of y must be 0 or 1")
  if(is.list(G)) {
    message("[plaid.test] converting gmt to sparse matrix...")
    G <- gmt2mat(G)
  } else {
    ## message("[plaid.test] sparse matrix provided")
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
    message("[plaid.test] computing one-sample t-tests on logFC")    
    mt1 <- matrix_onesample_ttest(fc, G)
    p1 <- mt1$p[,1]
    df1 <- mt1$mean[,1]
  }
  if("two" %in% tests) {
    message("[plaid.test] computing two-sample t-tests on logFC")
    mt2 <- matrix_twosample_ttest(fc, G)
    p2 <- mt2$p[,1]
    df2 <- mt2$diff[,1]        
  }
  if("lm" %in% tests) {
    if(is.null(gsetX)) {      
      message("[plaid.test] computing plaid scores...")
      gsetX <- plaid(X, G)
    }
    message("[plaid.test] computing gsetX t-tests")
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
  
  for(i in 1:length(P)) {
    P1 <- P[[i]]
    P1[is.na(P1)] <- 1
    P1 <- pmin(pmax(P1, 1e-99), 1-1e-99)
    P[[i]] <- P1
  }

  gg <- Reduce(intersect, lapply(P, names))
  P <- lapply(P, function(x) x[gg])
  F <- lapply(F, function(x) x[gg])
  F <- do.call( cbind, F)
  ##F[is.na(F)] <- 0
  gsetFC <- rowMeans(F)
  
  if(length(P)>1) {
    message("[plaid.test] computing meta-p...")
    pmeta <- matrix_combine_p(P, method=metap.method) 
  } else {
    pmeta <- P[[1]]
  }
  P <- do.call(cbind, P)
  colnames(P) <- paste0("p.",colnames(P))
  qmeta <- stats::p.adjust(pmeta, method="fdr")
  res <- cbind(
    gsetFC = gsetFC,
    pvalues = P,
    p.meta = pmeta,
    q.meta = qmeta    
  )
  if(sort.by %in% colnames(res)) {
    res <- res[order(res[,sort.by]),]
  }
  res
}

matrix_onesample_ttest <- function(F, G) {  
  sumG <- Matrix::colSums(G!=0)
  sum_sq  <- Matrix::crossprod(G!=0, F^2) 
  meanx <- Matrix::crossprod(G!=0, F) / (1e-8 + sumG)
  sdx   <-  sqrt( (sum_sq - meanx^2 * sumG) / (sumG - 1))
  f_stats <- meanx
  t_stats <- meanx / (1e-8 + sdx) * sqrt(sumG)
  p_stats <- apply( abs(t_stats), 2, function(tv)
    2*pt(tv,df=pmax(sumG-1,1),lower.tail=FALSE))
  list(mean = as.matrix(f_stats), t = as.matrix(t_stats), p = p_stats)  
}

matrix_twosample_ttest <- function(F, G) {
  if(is.vector(F)) F <- cbind(F)
  if(nrow(F)!=nrow(G)) stop("dimension mismatch")
  ## see e.g. https://people.umass.edu/bwdillon/.../TwoSampleT-Test.html
  sum1 <- Matrix::colSums(G!=0)
  # sum0 <- Matrix::colSums(G==0)  
  sum0 <- nrow(G) - sum1

  F2 <- F^2
  sum.F2 <- Matrix::colSums(F2)
  ssq1 <- Matrix::crossprod(G!=0, F2)     
  #ssq0 <- Matrix::crossprod(G==0, F2)
  ssq0 <- sweep(-ssq1, 2, sum.F2, '+') # faster

  sum.F <- Matrix::colSums(F)
  mean1 <- Matrix::crossprod(G!=0, F) 
  #mean0 <- Matrix::crossprod(G==0, F) 
  mean0 <- sweep(-mean1, 2, sum.F, '+') 
  mean1 <- mean1 / (1e-8 + sum1)
  mean0 <- mean0 / (1e-8 + sum0)    
  
  var0 <-  (ssq0 - mean0^2 * sum0) / (sum0 - 1)
  var1 <-  (ssq1 - mean1^2 * sum1) / (sum1 - 1)  
  varsum <- ( var0 / sum0 + var1 / sum1 )
  dof <- varsum^2 / ( var0 / sum0 * (sum0-1) + var1 / sum1 * (sum1 - 1) )
  ## NEED CHECKING!!!!
  f_stats <- mean1 - mean0
  t_stats <- f_stats / sqrt(varsum)
  p_stats <- sapply( 1:NCOL(F), function(i)
    2 * pt( abs(t_stats[,i]), df = pmax(dof[,i],1), lower.tail=FALSE))
  res <- list(diff = as.matrix(f_stats), t = as.matrix(t_stats), p = p_stats)
  res
}

matrix_combine_p <- function(plist, method='fisher') {
  if(method %in% c("fisher","sumlog")) {
    chisq <- (-2) * Reduce('+', lapply(plist,log))
    df <- 2 * length(plist)
    pv <- pchisq(chisq, df, lower.tail=FALSE)
  } else if(method %in% c("stouffer","sumz")) {
    np <- length(plist)
    zz <- lapply(plist, qnorm, lower.tail=FALSE) 
    zz <- Reduce('+', zz) / sqrt(np)
    pv <- pnorm(zz, lower.tail=FALSE)
  } else {
    stop("Invalid method: ",method)
  }
  dimnames(pv) <- dimnames(plist[[1]])
  return(pv)
}


##----------------------------------------------------------------
##-------------------- UTILITIES ---------------------------------
##----------------------------------------------------------------

#' Normalize column medians of matrix
#'
#' This function normalizes the column medians of matrix x. It calls
#' optimized functions from the matrixStats package.
#'
#' @param x Input matrix
#' @param ignore.zero Logical indicating whether to ignore zeros to
#'   exclude for median calculation
#'
#' @export
normalize_medians <- function(x, ignore.zero = NULL) {

  if(is.null(ignore.zero))
    ignore.zero <- (min(x,na.rm = TRUE) == 0)

  x <- as.matrix(x)

  if(ignore.zero) {
    zx <- x
    zx[Matrix::which(x==0)] <- NA
    #medx <- Rfast::colMedians(zx, na.rm=TRUE)
    medx <- matrixStats::colMedians(zx, na.rm = TRUE)    
    medx[is.na(medx)] <- 0
  } else {
    #medx <- Rfast::colMedians(x, na.rm=TRUE)
    medx <- matrixStats::colMedians(x, na.rm = TRUE)    
  }

  nx <- sweep(x, 2, medx, '-') + mean(medx, na.rm = TRUE)
  return(nx)
  
}

#' Compute columnwise ranks of matrix
#'
#' Computes columnwise rank of matrix. Can be sparse. Tries to call
#' optimized functions from Rfast or matrixStats.
#'
#' @param X Input matrix
#' @param sparse Logical indicating to use sparse methods
#' @param signed Logical indicating using signed ranks
#' @param keep.zero Logical indicating whether to keep zero as ranked zero
#' @param ties.method Character Choice of ties.method
#' 
#' @export
colranks <- function(X,
                     sparse = NULL,
                     signed = FALSE,
                     keep.zero = FALSE,
                     ties.method = "average") {

  if(is.null(sparse))
    sparse <- inherits(X,"CsparseMatrix")

  if(sparse) {
    X <- methods::as(X, "CsparseMatrix")
    if(keep.zero) {
      rX <- sparse_colranks(X, signed = signed, ties.method = ties.method)
    } else {
      if(signed) {
        sign.X <- sign(X)
        abs.rX <- Matrix::t(sparseMatrixStats::colRanks(abs(X), ties.method = ties.method))
        rX <- abs.rX * sign.X
      } else {
        rX <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method = ties.method))
      }
    }
  } else {
    if(signed) {
      sign.X <- sign(X)
      abs.rX <- Matrix::t(matrixStats::colRanks(as.matrix(abs(X)), ties.method = ties.method))
      rX <- sign.X * abs.rX
    } else {
      rX <- Matrix::t(matrixStats::colRanks(as.matrix(X), ties.method = ties.method))
    }
  }

  return(rX)

}

#' Compute columm ranks for sparse matrix. Internally used by colranks()
#'
#' @param X Input matrix
#' @param signed Logical: use or not signed ranks
#' @param ties.method Character Choice of ties.method
#' 
sparse_colranks <- function(X, signed = FALSE, ties.method = "average") {
  ## https://stackoverflow.com/questions/41772943
  X <- methods::as(X, "CsparseMatrix")
  n <- diff(X@p)  ## number of non-zeros per column
  lst <- split(X@x, rep.int(1:ncol(X), n))  ## columns to list
  ## column-wise ranking and result collapsing
  if(signed) {
    lst.sign <- lapply(lst, sign)
    lst.rnk  <- lapply(lst, function(x) rank(abs(x),ties.method = ties.method))
    rnk <- unlist(mapply('*', lst.sign, lst.rnk, SIMPLIFY = FALSE))  
  } else {
    rnk <- unlist(lapply(lst, rank, ties.method = ties.method))  
  }

  rX <- X  ## copy sparse matrix
  rX@x <- rnk  ## replace non-zero elements with rank

  return(rX)

}

##-------------------------------------------------------------
##------------------ end of file ------------------------------
##-------------------------------------------------------------
