##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @importFrom methods as
#' @importFrom stats ecdf
#' @importFrom Matrix colSums colScale crossprod Diagonal rowMeans t which
#' @importFrom matrixStats colMedians
#' @importFrom MatrixGenerics colRanks rowSds
NULL

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
#' @details PLAID now automatically detects and handles Bioconductor objects.
#' If X is a SummarizedExperiment or SingleCellExperiment, it will extract
#' the appropriate assay. If matG is a BiocSet object or GMT list, it will
#' be converted to sparse matrix format automatically.
#' 
#' @param X Log-transformed expr. matrix. Genes on rows, samples on columns.
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#'   Also accepts BiocSet objects or GMT lists (named list of gene vectors).
#' @param stats Score computation stats: mean or sum of intensity. Default 'mean'.
#' @param chunk Logical: use chunks for large matrices. Default 'NULL' for autodetect.
#' @param normalize Logical: median normalize results or not. Default 'TRUE'.
#' @param nsmooth Smoothing parameter for more stable average when stats="mean". Default 3.
#' @param assay Character: assay name to extract from SummarizedExperiment/SingleCellExperiment. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set (for BiocSet/GMT input). Default 5.
#' @param max.genes Integer: maximum genes per gene set (for BiocSet/GMT input). Default 500.
#'
#' @return Matrix of single-sample enrichment scores.
#' Gene sets on rows, samples on columns.
#' 
#' @examples
#' library(plaid)
#' 
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:100)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:20),
#'   "Pathway2" = paste0("GENE", 15:35),
#'   "Pathway3" = paste0("GENE", 30:50)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute PLAID scores
#' gsetX <- plaid(X, matG)
#' print(dim(gsetX))
#' print(gsetX[1:3, 1:5])
#' 
#' # Use sum statistics instead of mean
#' gsetX_sum <- plaid(X, matG, stats = "sum")
#' 
#' \donttest{
#' # Using real data (if available in package)
#' extdata_path <- system.file("extdata", "pbmc3k-50cells.rda", package = "plaid")
#' if (file.exists(extdata_path)) {
#'   load(extdata_path)
#'   hallmarks <- system.file("extdata", "hallmarks.gmt", package = "plaid")
#'   gmt <- read.gmt(hallmarks)
#'   matG <- gmt2mat(gmt)
#'   gsetX <- plaid(X, matG)
#' }
#' }
#'
#' @export
plaid <- function(X, matG, stats=c("mean","sum"), chunk=NULL, normalize=TRUE,
                  nsmooth=3, assay="logcounts", min.genes=5, max.genes=500) {

  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)

  stats <- stats[1]
  if (NCOL(X) == 1) X <- cbind(X)

  ## make sure X is matrix (not dataframe) and convert to sparse if needed
  if(!inherits(X, "Matrix")) {
    X <- as.matrix(X)
    if(mean(X==0,na.rm=TRUE)>0.5) X <- Matrix::Matrix(X, sparse=TRUE)
  }
  
  gg <- intersect(rownames(X), rownames(matG))
  if (length(gg) == 0) {
    message("[plaid] No overlapping features.")
    return(NULL)
  }

  X <- X[gg, , drop = FALSE]
  matG <- matG[gg, , drop = FALSE]
  G <- 1 * (matG != 0)
  if(stats == "mean") {
    # nsmooth is 'smoothing' parameter for more stable average
    sumG <- 1e-8 + nsmooth + Matrix::colSums(G, na.rm = TRUE)
    G <- Matrix::colScale(G, 1 / sumG)
  }

  ## Calculates PLAID score
  gsetX <- chunked_crossprod(G, X, chunk=NULL)
  gsetX <- as.matrix(gsetX)
  
  if(normalize) {
    if(nrow(gsetX) < 20) {
      ## for few genesets, median-norm is not good
      normfactor <- Matrix::colMeans(gsetX,na.rm=TRUE)
      normfactor <- normfactor - mean(normfactor)
      gsetX <- sweep(gsetX, 2, normfactor, '-') 
    } else {
      gsetX <- normalize_medians(gsetX)
    }
  }
 
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

  for(i in seq_len(k)) {
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
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns. Also accepts BiocSet objects or GMT lists.
#' @param removeLog2 Logical for whether to remove the Log2, i.e. will
#'   apply power transform (base2) on input (default TRUE).
#' @param scoreMean Logical for whether computing sum or mean as score
#'   (default FALSE).
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#'
#' @return Matrix of single-sample scSE enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix (log-transformed)
#' set.seed(123)
#' X <- log2(matrix(rpois(500, lambda = 10) + 1, nrow = 50, ncol = 10))
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute scSE scores (original method)
#' scores <- replaid.scse(X, matG, removeLog2 = TRUE, scoreMean = FALSE)
#' print(scores[1:2, 1:5])
#' 
#' # Compute scSE scores (mean method)
#' scores_mean <- replaid.scse(X, matG, removeLog2 = TRUE, scoreMean = TRUE)
#' print(scores_mean[1:2, 1:5])
#'
#' @export
replaid.scse <- function(X,
                         matG,
                         removeLog2 = NULL,
                         scoreMean = FALSE,
                         assay="logcounts",
                         min.genes=5,
                         max.genes=500) {

  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)

  if(is.null(removeLog2))
    removeLog2 <- min(X, na.rm = TRUE)==0 && max(X, na.rm=TRUE) < 20
  
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
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns. Also accepts BiocSet objects or GMT lists.
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#' 
#' @return Matrix of single-sample singscore enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute singscore
#' scores <- replaid.sing(X, matG)
#' print(scores[1:2, 1:5])
#'
#' @export
replaid.sing <- function(X, matG, assay="logcounts", min.genes=5, max.genes=500) {
  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)
  
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
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns. Also accepts BiocSet objects or GMT lists.
#' @param alpha Weighting factor for exponential weighting of ranks
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#' 
#' @return Matrix of single-sample ssGSEA enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute ssGSEA scores (alpha = 0)
#' scores <- replaid.ssgsea(X, matG, alpha = 0)
#' print(scores[1:2, 1:5])
#' 
#' # Compute ssGSEA scores with weighting (alpha = 0.25)
#' scores_weighted <- replaid.ssgsea(X, matG, alpha = 0.25)
#' print(scores_weighted[1:2, 1:5])
#'
#' @export
replaid.ssgsea <- function(X, matG, alpha = 0, assay="logcounts", min.genes=5, max.genes=500) {
  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)
  
  rX <- colranks(X, keep.zero = FALSE, ties.method = "average")
  if(alpha != 0) {
    ## This is not exactly like original formula. Not sure how to
    ## efficiently implement original rank weighting
    rX <- rX^(1 + alpha)
  }
  ## center and normalize
  rX <- rX - mean(rX,na.rm=TRUE)
  rX <- (rX / max(abs(rX),na.rm=TRUE))
  ## add scaling (to match NES scale)
  qx <- quantile(abs(rX),probs=0.99)
  rX <- 0.9 * rX / qx  
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
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#'   Also accepts BiocSet objects or GMT lists.
#' @param rmax Rank threshold (see Ucell paper). Default rmax = 1500.
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#' 
#' @return Matrix of single-sample UCell enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute UCell scores (default rmax = 1500)
#' scores <- replaid.ucell(X, matG)
#' print(scores[1:2, 1:5])
#' 
#' # Compute UCell scores with custom rmax
#' scores_custom <- replaid.ucell(X, matG, rmax = 1000)
#' print(scores_custom[1:2, 1:5])
#'
#' @export
replaid.ucell <- function(X, matG, rmax = 1500, assay="logcounts", min.genes=5, max.genes=500) {
  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)
  
  rX <- colranks(X, ties.method = "average")
  rX <- pmin( max(rX, na.rm=TRUE) - rX, rmax+1 )
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
#' from the original `AUCell` R package and we showed good concordance
#' of results in the score, logFC and p-values.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on columns.
#'   Also accepts BiocSet objects or GMT lists.
#' @param aucMaxRank Rank threshold (see AUCell paper). Default aucMaxRank = 0.05*nrow(X).
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#' 
#' @return Matrix of single-sample AUCell enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute AUCell scores
#' scores <- replaid.aucell(X, matG)
#' print(scores[1:2, 1:5])
#'
#' @export
replaid.aucell <- function(X, matG, aucMaxRank = NULL, assay="logcounts", min.genes=5, max.genes=500) {
  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)
  
  if (is.null(aucMaxRank)) {
    aucMaxRank <- ceiling(0.05*nrow(X))
  }
  
  rX <- colranks(X, ties.method = "average")
  ww <- 1.08*pmax((rX - (max(rX, na.rm=TRUE) - aucMaxRank)) / aucMaxRank, 0)
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
#' transform expression of each feature to a (0;1) relative expression
#' value. For efficiency reasons, this is here approximated by a
#' z-transform (center+scale) of each row.
#' 
#' @param X Gene or protein expression matrix. Generally log
#'   transformed. See details. Genes on rows, samples on columns.
#'   Also accepts SummarizedExperiment or SingleCellExperiment objects.
#' @param matG Gene sets sparse matrix. Genes on rows, gene sets on
#'   columns. Also accepts BiocSet objects or GMT lists.
#' @param tau Rank weight parameter (see GSVA publication). Default
#'   tau=0.
#' @param rowtf Row transformation method ("z" or "ecdf"). Default "z".
#' @param assay Character: assay name for Bioconductor objects. Default "logcounts".
#' @param min.genes Integer: minimum genes per gene set. Default 5.
#' @param max.genes Integer: maximum genes per gene set. Default 500.
#' 
#' @return Matrix of single-sample GSVA enrichment scores.
#'   Gene sets on rows, samples on columns.
#'
#' @examples
#' # Create example expression matrix
#' set.seed(123)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' rownames(X) <- paste0("GENE", 1:50)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Create example gene sets
#' gmt <- list(
#'   "Pathway1" = paste0("GENE", 1:15),
#'   "Pathway2" = paste0("GENE", 10:25)
#' )
#' matG <- gmt2mat(gmt)
#' 
#' # Compute GSVA scores
#' scores <- replaid.gsva(X, matG)
#' print(scores[1:2, 1:5])
#'
#' @export
replaid.gsva <- function(X, matG, tau = 0, rowtf = c("z", "ecdf")[1], assay="logcounts", min.genes=5, max.genes=500) {
  ## Auto-detect and convert Bioconductor objects
  if (inherits(X, "SummarizedExperiment") || inherits(X, "SingleCellExperiment")) {
    X <- .extract_expression_matrix(X, assay=assay, log.transform=FALSE)
  }
  
  matG <- .convert_geneset_to_matrix(matG, background=rownames(X),
                                     min.genes=min.genes, max.genes=max.genes)
  
  rowtf <- rowtf[1]

  if(rowtf == "z") {
    ## Faster approximation of relative activation
    zX <- (X - Matrix::rowMeans(X)) / (1e-8 + mat.rowsds(X))
  } else if(rowtf=='ecdf') {
    ## this implements original ECDF idea
    zX <- t(apply(X,1,function(x) ecdf(x)(x))) 
  } else {
    stop("unknown row transform",rowtf)
  }

  rX <- colranks(zX, signed = FALSE, ties.method = "average")
  rX <- rX - mean(rX,na.rm=TRUE)
  rX <- rX / max(abs(rX), na.rm=TRUE)
  if(tau > 0) {
    ## Note: This is not exactly like original formula. Not sure how
    ## to efficiently implement original rank weighting
    rX <- sign(rX) * abs(rX)^(1 + tau)
  }
  dimnames(rX) <- dimnames(X)
  gsetX <- plaid(rX, matG)

  return(gsetX)

}

#' Calculate row standard deviations for matrix
#'
#' @param X Input matrix (can be sparse or dense)
#'
#' @return Vector of row standard deviations.
#'
mat.rowsds <- function(X) {
  ## MatrixGenerics::rowSds dispatches to the right method based on class
  MatrixGenerics::rowSds(X)
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
#' @return Matrix with normalized column medians.
#'
#' @examples
#' # Create example matrix
#' set.seed(123)
#' x <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' x[1:3, 1:3] <- 0  # Add some zeros
#' 
#' # Normalize medians
#' x_norm <- normalize_medians(x)
#' head(x_norm)
#'
#' @export
normalize_medians <- function(x, ignore.zero = NULL) {

  if(is.null(ignore.zero))
    ignore.zero <- (min(x,na.rm = TRUE) == 0)

  x <- as.matrix(x)

  if(ignore.zero) {
    zx <- x
    zx[Matrix::which(x==0)] <- NA
    medx <- matrixStats::colMedians(zx, na.rm = TRUE)    
    medx[is.na(medx)] <- 0
  } else {
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
#' @return Matrix of columnwise ranks with same dimensions as input.
#'
#' @examples
#' # Create example matrix
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' rownames(X) <- paste0("Gene", 1:10)
#' colnames(X) <- paste0("Sample", 1:10)
#' 
#' # Compute column ranks
#' ranks <- colranks(X)
#' print(ranks[1:5, 1:5])
#' 
#' # Compute signed ranks
#' signed_ranks <- colranks(X, signed = TRUE)
#' print(signed_ranks[1:5, 1:5])
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
        abs.rX <- Matrix::t(MatrixGenerics::colRanks(abs(X), ties.method = ties.method))
        rX <- abs.rX * sign.X
      } else {
        rX <- Matrix::t(MatrixGenerics::colRanks(X, ties.method = ties.method))
      }
    }
  } else {
    if(signed) {
      sign.X <- sign(X)
      abs.rX <- Matrix::t(MatrixGenerics::colRanks(as.matrix(abs(X)), ties.method = ties.method))
      rX <- sign.X * abs.rX
    } else {
      rX <- Matrix::t(MatrixGenerics::colRanks(as.matrix(X), ties.method = ties.method))
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
#' @return Sparse matrix of columnwise ranks with same dimensions as input.
#'
sparse_colranks <- function(X, signed = FALSE, ties.method = "average") {
  ## https://stackoverflow.com/questions/41772943
  X <- methods::as(X, "CsparseMatrix")
  n <- diff(X@p)  ## number of non-zeros per column
  lst <- split(X@x, rep.int(seq_len(ncol(X)), n))  ## columns to list
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
