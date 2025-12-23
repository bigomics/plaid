##----------------------------------------------------------------
##----------------- STATISTICAL TESTS ----------------------------
##----------------------------------------------------------------

#' @importFrom stats p.adjust pt pchisq qnorm pnorm
#' @importFrom Matrix rowMeans colSums crossprod t colScale
#' @importFrom MatrixGenerics rowSds colVars colRanks
#' @importFrom Rfast ttests
#' @importFrom qlcMatrix corSparse
#' @importFrom methods is
NULL

#' Reimplementation of dualGSEA (Bull et al., 2024) but defaults with
#' replaid backend. For the preranked test we still use fgsea. Should
#' be much faster than original using fgsea + GSVA::ssGSEA.
#'
#' @param X Expression matrix with genes on rows and samples ont columns
#' @param y Binary vector (0/1) indicating group membership
#' @param gmt List of gene sets in GMT format
#' @param G Sparse matrix of gene sets. Non-zero entry indicates
#'   gene/feature is part of gene sets. Features on rows, gene sets on
#'   columns.
#' @param gsetX Optional pre-computed matrix of gene set enrichment scores with 
#'   gene sets on rows and samples on columns. If NULL (default), scores will be 
#'   computed using the method specified by `ss.method`. Providing pre-computed 
#'   scores improves efficiency when running multiple analyses.
#' @param fc.method Method for fold change testing ("fgsea", "ztest", "ttest", "rankcor", "cor")
#' @param ss.method Method for single-sample enrichment ("plaid", "replaid.ssgsea", "replaid.gsva", "ssgsea", "gsva")
#' @param pv1 Pre-computed p-values from fold change test. If NULL, will be computed based on fc.test.
#' @param pv2 Pre-computed p-values from single sample test. If NULL, will be computed using gset_ttest.
#' @param metap.method Method for combining p-values ("stouffer", "fisher" or "maxp"). Default "stouffer".
#' @param sort.by Column name to sort results by ("p.dual", "gsetFC", "p.fc", "p.ss"). Default "p.dual".
#' 
#' @return Data frame with columns: gsetFC (gene set fold change), size (gene set size),
#'   p.fc (p-value from fold change test), p.ss (p-value from single sample test),
#'   p.dual (combined p-value), and q.dual (FDR-adjusted combined p-value).
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(1000), nrow = 100, ncol = 20)
#' rownames(X) <- paste0("GENE", 1:100)
#' colnames(X) <- paste0("Sample", 1:20)
#' 
#' # Create binary group vector
#' y <- rep(c(0, 1), each = 10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:20),
#'   "Pathway2" = paste0("GENE", 15:35),
#'   "Pathway3" = paste0("GENE", 30:50)
#' )
#' 
#' # Perform dualGSEA with correlation test (fast method)
#' results_cor <- dualGSEA(X, y, G = NULL, gmt = gmt, fc.method = "cor", ss.method = "replaid.gsva")
#' print(head(results_cor))
#' 
#' \donttest{
#' # Perform dualGSEA with fgsea (requires fgsea package)
#' if (requireNamespace("fgsea", quietly = TRUE)) {
#'   results <- dualGSEA(X, y, G = NULL, gmt = gmt, fc.method = "fgsea", ss.method = "replaid.ssgsea")
#'   print(head(results))
#' }
#' }
#'
#' @export
dualGSEA <- function(X, y, G, gmt=NULL, gsetX=NULL,
                     fc.method = c("fgsea","rankcor","ztest","ttest","cor")[2],
                     ss.method = c('plaid', 'replaid.ssgsea','replaid.gsva', 
                       'ssgsea','gsva')[1],
                     metap.method = c("stouffer","fisher","maxp")[1],
                     pv1 = NULL, pv2 = NULL, sort.by='p.dual') {
  #require(fgsea)
  if (fc.method == "fgsea" && !requireNamespace("fgsea", quietly=TRUE)) {
    stop("The fgsea package must be installed to use this functionality")
  }
  if (ss.method %in% c("ssgsea","gsva") && !requireNamespace("GSVA", quietly=TRUE)) {
    stop("The GSVA package must be installed to use this functionality")
  }
  if(is.null(gmt) && is.null(G)) {
    stop("at least gmt or matrix G must be given")
  }
  if(!is.null(gmt) && !is(gmt, "list")) {
    stop("gmt must be a list")
  }
  
  if(!all(unique(y) %in% c(0,1,NA))) stop("elements of y must be 0 or 1")
  sel <- which(!is.na(y))
  y <- y[sel]
  X <- X[,sel,drop=FALSE] 
  
  if(is.null(G) && !is.null(gmt))  G <- gmt2mat(gmt)
  if(is.null(gmt) && !is.null(G))  gmt <- mat2gmt(G)  
  
  ## pairwise test on logFC
  if(is.null(pv1)) {
    message("FC testing using ", fc.method)
    m1 <- Matrix::rowMeans(X[,which(y==1),drop=FALSE])
    m0 <- Matrix::rowMeans(X[,which(y==0),drop=FALSE])
    fc <- m1 - m0
    if(fc.method == "fgsea") {
      res1 <- fgsea::fgsea(gmt, fc)
      res1 <- data.frame(res1, row.names=res1$pathway)
      pv1 <- res1[,"pval"]
      names(pv1) <- rownames(res1)
    } else     if(fc.method %in% c('ttest','ztest')) {
      ## MatrixGenerics::rowSds dispatches to the right method based on class
      sdx <- MatrixGenerics::rowSds(X, na.rm=TRUE)
      sdx0 <- mean(sdx, na.rm=TRUE)
      zc <- fc / (0.1*sdx0 + sdx)
      if(fc.method == "ttest") {
        res1 <- fc_ttest(zc, G, sort.by="none")
        pv1 <- res1[,'pvalue']
      }
      if(fc.method == "ztest") {
        res1 <- fc_ztest(zc, G, zmat=FALSE)
        pv1 <- res1$p_value
      }
    } else if(fc.method == 'rankcor') {
      res1 <- gset.rankcor(fc, G, compute.p=TRUE, use.rank=TRUE)
      pv1 <- res1$p.value[,1]
    } else if(fc.method == 'cor') {
      res1 <- gset.rankcor(fc, G, compute.p=TRUE, use.rank=FALSE)
      pv1 <- res1$p.value[,1]
    } else {
      stop("invalid fc.method method")
    }
  }
  
  ## single-sample test
  message("single-sample testing using ", ss.method)

  if(is.null(gsetX)) {
    if(ss.method == "plaid") {
      gsetX <- plaid::plaid(X, G)
    }else if(ss.method == "gsva") {
      gsvapar <- GSVA::gsvaParam(X, gmt)
      gsetX <- GSVA::gsva(gsvapar, verbose = FALSE)
    } else if(ss.method == "ssgsea") {
      gsvapar <- GSVA::ssgseaParam(X, gmt)
      gsetX <- GSVA::gsva(gsvapar, verbose = FALSE)
    } else if(ss.method == "replaid.ssgsea") {
      gsetX <- replaid.ssgsea(X, G)
    } else if(ss.method == "replaid.gsva") {
      gsetX <- replaid.gsva(X, G)
    } else {
      stop("invalid ss.method: ",ss.method)
    }
  } else {
    gsetX <- gsetX[colnames(G),]
  }

  gg <- intersect(rownames(G),rownames(X))
  sel <- which(!is.na(y))
  X <- X[gg,sel,drop=FALSE]
  G <- G[gg,]
  y <- y[sel]
  
  if(is.null(gsetX)) {
    message("computing gsetX using plaid. please precompute for efficiency.")
    gsetX <- plaid(X, G)
  }
  pp <- intersect(rownames(gsetX),colnames(G))
  gsetX <- gsetX[pp,colnames(X),drop=FALSE]
  G <- G[,pp,drop=FALSE]  
  
  e1 <- Matrix::rowMeans(gsetX[,y==1,drop=FALSE])
  e0 <- Matrix::rowMeans(gsetX[,y==0,drop=FALSE])
  gsetFC <- e1 - e0    
  gs.size <- Matrix::colSums(G!=0)[rownames(gsetX)]

  if(is.null(pv2)) {
    res2 <- gset_ttest(gsetX, y) 
    pv2 <- res2[,'pvalue']
  }
  
  gs <- rownames(gsetX)
  P <- cbind(pv1[gs], pv2[gs])  
  P[is.na(P)] <- 1
  P <- pmin(pmax(P, 1e-99), 1-1e-99)
  colnames(P) <- c("p.fc","p.ss")

  p.dual <- matrix_metap(P, method=metap.method)   
  q.dual <- stats::p.adjust(p.dual, method="fdr")
  
  res <- cbind(
    gsetFC = gsetFC,
    size = gs.size,
    pvalues = P,
    p.dual = p.dual,
    q.dual = q.dual    
  )
  if(sort.by %in% colnames(res)) {
    osign <- ifelse(sort.by=="gsetFC",-1,1)
    res <- res[order(osign*res[,sort.by]),]
  }
  res
}


#' T-test statistical testing of differentially enrichment
#'
#' This function performs statistical testing for differential
#' enrichment using plaid
#'
#' @param fc Vector of logFC values
#' @param G Sparse matrix of gene sets. Non-zero entry indicates
#'   gene/feature is part of gene sets. Features on rows, gene sets on
#'   columns.
#' @param sort.by Column name to sort results by ("pvalue", "gsetFC", or "none")
#' 
#' @return Data frame with columns: gsetFC (gene set fold change),
#'   pvalue (p-value from one-sample t-test), and qvalue (FDR-adjusted p-value).
#'
fc_ttest <- function(fc, G, sort.by="pvalue") {
  if(is.null(names(fc))) stop("fc must have names")  
  gg <- intersect(rownames(G),names(fc))
  fc <- fc[gg]
  G <- G[gg,]
  
  message("[fc_ttest] computing one-sample t-tests on logFC")    
  mt <- matrix_onesample_ttest(fc, G)
  pv <- mt$p[,1]
  df <- mt$mean[,1]
  qv <- p.adjust(pv, method="fdr")
  gsetFC <- gset_averageCLR(fc, G, center = FALSE)[,1]
  
  res <- cbind(
    gsetFC = gsetFC,
    pvalue = pv,
    qvalue = qv    
  )
  if(!is.null(sort.by) && sort.by %in% colnames(res)) {
    sort.sign <- ifelse(sort.by=="gsetFC",-1,+1)
    res <- res[order(sort.sign*res[,sort.by]),]
  }
  res
}

#' Z-test statistical testing of differentially enrichment
#'
#' This function performs statistical testing for differential
#' enrichment using plaid
#' 
#' @importFrom stats var
#'
#' @param fc Vector of logFC values
#' @param G Sparse matrix of gene sets. Non-zero entry indicates
#'   gene/feature is part of gene sets. Features on rows, gene sets on
#'   columns.
#' @param zmat Logical indicating to return z-matrix
#' @param alpha Scalar weight for SD estimation. Default 0.5.
#' 
#' @return List with element: z_statistic (z-statistic from one-sample z-test),
#'   p_value (p-value from one-sample z-test), and zmat (z-matrix).
#'
fc_ztest <- function(fc, G, zmat=FALSE, alpha=0.5) {
  if(is.null(names(fc))) stop("fc must have names")
  gg <- intersect(rownames(G),names(fc))
  sample_size <- Matrix::colSums(G[gg,]!=0)
  sample_size <- pmax(sample_size, 1) ## avoid div-by-zero
  sample_mean <- (Matrix::t(G[gg,]!=0) %*% fc[gg]) / sample_size
  population_mean <- mean(fc, na.rm=TRUE)
  population_var <- var(fc, na.rm=TRUE)
  gfc <- (G[gg,]!=0) * fc[gg]
  ## MatrixGenerics::colVars dispatches to the right method based on class
  sample_var <- MatrixGenerics::colVars(gfc) * nrow(G) / sample_size
  alpha <- pmin(pmax(alpha,0), 0.999) ## limit
  estim_sd <- sqrt( alpha*sample_var + (1-alpha)*population_var )
  z_statistic <- (sample_mean - population_mean) / (estim_sd / sqrt(sample_size))
  p_value <- 2 * pnorm(abs(z_statistic[,1]), lower.tail = FALSE)  
  if(zmat) {
    zmat <- (Matrix::t(gfc) / estim_sd)
  } else {
    zmat <- NULL
  }
  list(
    z_statistic = z_statistic[,1],
    p_value = p_value,
    zmat = zmat
  )
}

#' Compute geneset expression as the average log-ration of genes in
#' the geneset. Requires log-expression matrix X and (sparse) geneset
#' matrix matG.
#'
#' @param X Log-expression matrix with genes on rows and samples on columns
#' @param matG Sparse gene set matrix with genes on rows and gene sets on columns
#' @param center Logical indicating whether to center the results
#'
#' @return Matrix of gene set expression scores with gene sets on rows and samples on columns.
#'
gset_averageCLR <- function(X, matG, center = TRUE) {
  if (NCOL(X) == 1) X <- cbind(X)
  gg <- intersect(rownames(X), rownames(matG))
  if (length(gg) == 0) {
    message("[gset.averageCLR] no overlapping features")
    return(NULL)
  }
  X <- X[gg, , drop = FALSE]
  matG <- matG[gg, , drop = FALSE]
  sumG <- 1e-8 + Matrix::colSums(matG != 0, na.rm = TRUE)
  nG <- Matrix::colScale(1 * (matG != 0), 1 / sumG)
  gsetX <- Matrix::t(nG) %*% X 
  if(center) gsetX <- gsetX - Matrix::rowMeans(gsetX, na.rm = TRUE)
  as.matrix(gsetX)
}

#' Perform t-test on gene set scores
#'
#' @param gsetX Matrix of gene set scores with gene sets on rows and samples on columns
#' @param y Binary vector (0/1) indicating group membership
#'
#' @return Data frame with columns: diff (difference in means), statistic (t-statistic),
#'   pvalue (p-value), and other t-test results.
#'
gset_ttest <- function(gsetX, y) {
  ii <- which(!is.na(y))
  gsetX <- gsetX[,ii]
  y <- y[ii]
  if(!all(unique(y) %in% c(0,1))) stop("[gset_ttest] elements of y must be 0 or 1")
  res  <- Rfast::ttests(Matrix::t(gsetX), ina=y+1)
  rownames(res) <- rownames(gsetX)
  diff <- rowMeans(gsetX[,y==1]) - rowMeans(gsetX[,y==0])
  res <- cbind( diff=diff, res)
  return(res)
}

##----------------------------------------------------------------
##----------------- FUNCTIONS ------------------------------------
##----------------------------------------------------------------

#' Perform one-sample t-test on matrix with gene sets
#'
#' @param Fm Vector of feature values (e.g., fold changes)
#' @param G Sparse matrix of gene sets with genes on rows and gene sets on columns
#'
#' @return List containing mean, t-statistic, and p-value matrices.
#'
matrix_onesample_ttest <- function(Fm, G) {
  sumG <- Matrix::colSums(G!=0)
  sum_sq  <- Matrix::crossprod(G!=0, Fm^2) 
  meanx <- Matrix::crossprod(G!=0, Fm) / (1e-8 + sumG)
  sdx   <-  sqrt( (sum_sq - meanx^2 * sumG) / (sumG - 1))
  f_stats <- meanx
  t_stats <- meanx / (1e-8 + sdx) * sqrt(sumG)
  p_stats <- apply( abs(t_stats), 2, function(tv)
    2*pt(tv,df=pmax(sumG-1,1),lower.tail=FALSE))
  list(mean = as.matrix(f_stats), t = as.matrix(t_stats), p = p_stats)  
}

#' Matrix version for combining p-values using fisher or stouffer
#' method. Much faster than doing metap::sumlog() and metap::sumz()
#'
#' @param plist List of p-value vectors or matrix of p-values
#' @param method Method for combining p-values ("fisher"/"sumlog" or "stouffer"/"sumz")
#'
#' @return Vector of combined p-values.
#'
matrix_metap <- function(plist, method='stouffer') {
  if(inherits(plist,"matrix")) {
    plist <- as.list(data.frame(plist))
  }
  if(method %in% c("fisher","sumlog")) {
    chisq <- (-2) * Reduce('+', lapply(plist,log))
    df <- 2 * length(plist)
    pv <- pchisq(chisq, df, lower.tail=FALSE)
  } else if(method %in% c("stouffer","sumz")) {
    np <- length(plist)
    zz <- lapply(plist, qnorm, lower.tail=FALSE) 
    zz <- Reduce('+', zz) / sqrt(np)
    pv <- pnorm(zz, lower.tail=FALSE)
  } else if(method %in% c("maxp","pmax","maximump")) {
    pv <- Reduce(pmax, plist)    
  } else {
    stop("Invalid method: ",method)
  }
  dimnames(pv) <- dimnames(plist[[1]])
  return(pv)
}


#' Calculate gene set rank correlation
#'
#' Compute rank correlation between a gene rank vector/matrix and gene sets
#'
#' @param rnk Numeric vector or matrix of gene ranks, with genes as row names
#' @param gset Numeric matrix of gene sets, with genes as row/column names
#' @param compute.p Logical indicating whether to compute p-values
#' @param use.rank Logical indicating whether to rank transform rnk before correlation
#'
#' @return Named list with components:
#' \itemize{
#'  \item rho - Matrix of correlation coefficients between rnk and gset
#'  \item p.value - Matrix of p-values for correlation (if compute.p = TRUE)
#'  \item q.value - Matrix of FDR adjusted p-values (if compute.p = TRUE)
#' }
#'
#' @details This function calculates sparse rank correlation between rnk and each
#' column of gset using \code{qlcMatrix::corSparse()}. It handles missing values in
#' rnk by computing column-wise correlations.
#'
#' P-values are computed from statistical distribution
#'
#' @examples
#' # Create example rank vector
#' set.seed(123)
#' ranks <- rnorm(100)
#' names(ranks) <- paste0("GENE", 1:100)
#' 
#' # Create example gene sets as sparse matrix
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:20),
#'   "Pathway2" = paste0("GENE", 15:35),
#'   "Pathway3" = paste0("GENE", 30:50)
#' )
#' genesets <- gmt2mat(gmt)
#'
#' # Calculate rank correlation
#' result <- gset.rankcor(ranks, genesets, compute.p = TRUE)
#' print(result$rho)
#' print(result$p.value)
#' 
#' @export
gset.rankcor <- function(rnk, gset, compute.p = FALSE, use.rank = TRUE) {
  if (ncol(gset) == 0 || NCOL(rnk) == 0) {
    if (ncol(gset) == 0) message("gset has zero columns")
    if (NCOL(rnk) == 0) message("rnk has zero columns")
    return(NULL)
  }

  #  if (!any(class(gset) %in% c("Matrix", "dgCMatrix", "lgCMatrix", "matrix", "array"))) {
  #    stop("gset must be a matrix")
  #  }
  if (!inherits(gset, "Matrix")) stop("gset must be a matrix")

  is.vec <- (NCOL(rnk) == 1 && !any(class(rnk) %in% c("matrix", "Matrix")))
  if (is.vec && is.null(names(rnk))) stop("rank vector must be named")
  if (!is.vec && is.null(rownames(rnk))) stop("rank matrix must have rownames")
  if (is.vec) rnk <- matrix(rnk, ncol = 1, dimnames = list(names(rnk), "rnk"))
  n1 <- sum(rownames(rnk) %in% colnames(gset), na.rm = TRUE)
  n2 <- sum(rownames(rnk) %in% rownames(gset), na.rm = TRUE)
  if (n1 > n2) gset <- Matrix::t(gset)

  gg <- intersect(rownames(gset), rownames(rnk))
  rnk1 <- rnk[gg, , drop = FALSE]
  gset <- gset[gg, , drop = FALSE]

  if (use.rank) {
    ## MatrixGenerics::colRanks dispatches to the right method based on class
    rnk1 <- MatrixGenerics::colRanks(rnk1, na.last = "keep", ties.method = "random", preserveShape = TRUE)
  }

  ## two cases: (1) in case no missing values, just use corSparse on
  ## whole matrix. (2) in case the rnk matrix has missing values, we
  ## must proceed 1-column at time and do reduced corSparse on
  ## intersection of genes.
  rho1 <- cor_sparse_matrix(gset, rnk1)

  rownames(rho1) <- colnames(gset)
  colnames(rho1) <- colnames(rnk1)
  rho1[is.nan(rho1)] <- NA ## ??

  ## compute p-value
  .cor.pvalue <- function(x, n) 2 * stats::pnorm(-abs(x / ((1 - x**2) / (n - 2))**0.5))
  if (compute.p) {
    pv <- apply(rho1, 2, function(x) .cor.pvalue(x, n = nrow(rnk1)))
    pv[is.nan(pv)] <- NA ## ??
    qv <- apply(pv, 2, stats::p.adjust, method = "fdr")
    df <- list(rho = rho1, p.value = pv, q.value = qv)
  } else {
    df <- list(rho = rho1, p.value = NA, q.value = NA)
  }
  df
}

#' Calculate sparse correlation matrix handling missing values
#'
#' @param G Sparse matrix containing gene sets
#' @param mat Matrix of values
#' @return Correlation matrix between G and mat
#' @details If mat has no missing values, calculates correlation directly using corSparse.
#' Otherwise computes column-wise correlations only using non-missing values.
cor_sparse_matrix <- function(G, mat) {
  if (sum(is.na(mat)) == 0) {
    cor_matrix <- qlcMatrix::corSparse(G, mat)
  } else {
    message("matrix has missing values: computing column-wise reduced cor")
    corSparse.vec <- function(X, y) {
      jj <- which(!is.na(y))
      qlcMatrix::corSparse(X[jj, , drop = FALSE], cbind(y[jj]))
    }
    cor_matrix <- lapply(seq_len(ncol(mat)), function(i) corSparse.vec(G, mat[, i]))
    cor_matrix <- do.call(cbind, cor_matrix)
  }
  return(cor_matrix)
}
