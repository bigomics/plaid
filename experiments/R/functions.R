##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##
##

## This code is not striclty needed to run plaid but only for running
## other methods for comparison. No need to declare dependency on
## these packages.
##
##

if(!require(AUCell)) BiocManager::install("AUCell")
if(!require(UCell)) BiocManager::install("UCell")
if(!require(singscore)) BiocManager::install("singscore")
if(!require(GSVA)) BiocManager::install("GSVA")


##----------------------------------------------------------------
##-------------------- RUN OTHER METHODS -------------------------
##----------------------------------------------------------------

run.methods <- function(X, gmt) {
  ## prepare
  message("preparing sparse matrix matG...")
  matG <- gmt2mat(gmt)
  matG <- 1*(matG!=0)
  
  rX  <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="average")) / nrow(X)
  cX <- as.matrix(X - rowMeans(X, na.rm=TRUE))
  rcX <- Matrix::t(matrixStats::colRanks(cX, ties.method="min")) / nrow(X)

  message("preparing scSE folders...")
  prepare.SCSE(X, gmt, path="../scse")
  
  ## run
  message("running methods...")  
  timings <- peakRAM(
    rho    <- gset.rankcor(X, matG, use.rank=FALSE)$rho,
    rankcor <- gset.rankcor(X, matG, use.rank=TRUE)$rho,
    sing   <- gset.singscore(X, gmt, center=FALSE),
    gsva   <- gset.gsva(X, gmt, method="gsva"),
    ssgsea <- run.ssgsea(X, gmt, alpha=0.25),
    ucell  <- t(UCell::ScoreSignatures_UCell(X, gmt)),  ## needs logx
    aucell <- AUCell::getAUC(AUCell::AUCell_run(X, gmt)), ## uses rank
    scse <- run.SCSE(X, gmt, removeLog2=TRUE, scoreMean=FALSE, path="../scse"),
    scse.mean <- run.SCSE(X, gmt, removeLog2=FALSE, scoreMean=TRUE, path="../scse"),
    replaid.scse <- plaid::replaid.scse(X, matG, removeLog2=FALSE, scoreMean=TRUE),
    replaid.sing <- plaid::replaid.sing(X, matG),
    plaid  <- plaid::plaid(X, matG),
    plaid.r  <- plaid::plaid(rX, matG),
    plaid.c  <- plaid::plaid(cX, matG),
    plaid.rc  <- plaid::plaid(rcX, matG)    
  )
  rownames(ucell) <- sub("_UCell","",rownames(ucell))
  res <- list(
    cor = rho[,],
    rankcor = rankcor[,],
    sing = sing[,],
    gsva = gsva[,],
    ssgsea = ssgsea[,],    
    ucell = ucell,
    aucell = aucell,
    scse = scse,
    scse.mean = scse.mean,        
    replaid.scse = replaid.scse,
    replaid.sing = replaid.sing,        
    plaid = plaid,
    plaid.r = plaid.r,
    plaid.c = plaid.c,
    plaid.rc = plaid.rc        
  )
  res <- lapply(res, as.matrix)
  ss <- Reduce(intersect, lapply(res,rownames))
  res <- lapply(res, function(x) x[ss,])
  timings$Function_Call <- names(res)
  list(results = res, timings = timings)
}

## like run.methods but only for timings
run.timings <- function(X, gmt, timeout=3L,
                        methods=c("cor","rankcor","sing","gsva","ssgsea","ucell","aucell",
                                  "scse","scse.mean","replaid.scse","replaid.sing",
                                  "plaid","plaid.r","plaid.c","plaid.rc")) {
  ## prepare
  message("preparing sparse matrix matG...")
  matG <- gmt2mat(gmt)
  matG <- 1*(matG!=0)
  
  rX  <- Matrix::t(sparseMatrixStats::colRanks(X, ties.method="average")) / nrow(X)
  cX <- as.matrix(X - rowMeans(X, na.rm=TRUE))
  rcX <- Matrix::t(matrixStats::colRanks(cX, ties.method="min")) / nrow(X)

  if(any(c("scse","scse.mean") %in% methods)) {
    message("preparing scSE folders...")
    prepare.SCSE(X, gmt, path="../scse")
  }
  
  ## run with timeout
  timeout <- as.integer(timeout)
  wt <- function(m, ...) {
    if(!m %in% methods) return(NA)
    tryCatch({
      res <- R.utils::withTimeout( ... , timeout=timeout)
      return(FALSE)
    }, TimeoutException = function(ex) {
      message("Timeout")
      return(TRUE)
    })
  }
  
  message("running methods...")  
  res <- list()
  timings <- peakRAM(
    res[['cor']] <- wt("cor",gset.rankcor(X, matG, use.rank=FALSE)$rho),
    res[['rankcor']] <- wt("rankcor",rankcor <- gset.rankcor(X, matG, use.rank=TRUE)$rho),
    res[['sing']]   <- wt("sing",gset.singscore(X, gmt, center=FALSE)),
    res[['gsva']]   <- wt("gsva",gset.gsva(X, gmt, method="gsva")),
    res[['ssgsea']] <- wt("ssgsea",run.ssgsea(X, gmt, alpha=0.25)),
    res[['ucell']]  <- wt("ucell",t(UCell::ScoreSignatures_UCell(X, gmt))),  ## needs logx
    res[['aucell']] <- wt("aucell",AUCell::getAUC(AUCell::AUCell_run(X, gmt))), ## uses rank
    res[['scse']] <- wt("scse",run.SCSE(X, gmt, removeLog2=TRUE, scoreMean=FALSE, path="../scse")),
    res[['scse.mean']] <- wt("scse.mean",run.SCSE(X, gmt, removeLog2=FALSE, scoreMean=TRUE, path="../scse")),
    res[['replaid.scse']] <- wt("replaid.scse",plaid::replaid.scse(X, matG, removeLog2=FALSE, scoreMean=TRUE)),
    res[['replaid.sing']] <- wt("replaid.sing",plaid::replaid.sing(X, matG)),
    res[['plaid']]  <- wt("plaid",plaid::plaid(X, matG)),
    res[['plaid.r']]  <- wt("plaid.r",plaid::plaid(rX, matG)),
    res[['plaid.c']]  <- wt("plaid.c",plaid::plaid(cX, matG)),
    res[['plaid.rc']]  <- wt("plaid.rc",plaid::plaid(rcX, matG))    
  )

  res <- unlist(res)
  timings$Function_Call <- names(res)
  timings$Timeout <- res
  timings$Elapsed_Time_sec[is.na(res)] <- NA
  timings$Peak_RAM_Used_MiB[is.na(res)] <- NA  
  timings$Total_RAM_Used_MiB <- NULL

  return(timings)
}



#' Single sample genesets expression scores. See Foroutan 2018. This
#' seems to be a nice and efficient algorithm.
#'
#' @export
gset.gsva <- function(X, geneSets, method = "gsva") {
  if (!method %in% c("gsva", "ssgsea", "plage", "zscore")) {
    stop("ERROR. unknown method", method)
  }
  gsvapar <- switch(method,
    gsva = GSVA::gsvaParam(X, geneSets),
    plage = GSVA::plageParam(X, geneSets),
    zscore = GSVA::zscoreParam(X, geneSets),
    ssgsea = GSVA::ssgseaParam(X, geneSets)
  )
  es <- GSVA::gsva(gsvapar, verbose = FALSE)
  return(es)
}

#' Single sample genesets expression scores. See Foroutan 2018. This
#' seems to be a nice and efficient algorithm.
#'
#' @export
gset.singscore <- function(X, geneSets, center = FALSE, return.score = TRUE) {
  gslist <- lapply(geneSets, GSEABase::GeneSet)
  for (i in 1:length(gslist)) GSEABase::setName(gslist[[i]]) <- names(geneSets)[i]
  gscolc <- GSEABase::GeneSetCollection(gslist)
  if (center) {
    X <- X - rowMeans(X, na.rm = TRUE)
  }
  ranked <- singscore::rankGenes(as.matrix(X)) ## cannot take sparse
  sing <- singscore::multiScore(ranked, upSetColc = gscolc)
  ## sing has 'Scores' and 'Dispersions'
  if (return.score) {
    return(sing$Scores)
  }
  return(sing)
}

#' @export
gset.rankcor <- function(rnk, gset, compute.p = FALSE, use.rank = TRUE) {
  if (ncol(gset) == 0 || NCOL(rnk) == 0) {
    if (ncol(gset) == 0) message("ERROR. gset has zero columns")
    if (NCOL(rnk) == 0) message("ERROR. rnk has zero columns")
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
    if(!inherits(rnk1,"dgCMatrix")) {
      ## this doesnt work for sparse dgCMatrix
      rnk1 <- Matrix::t(matrixStats::colRanks(rnk1, na.last = "keep", ties.method = "random"))
    } else {
      rnk1 <- apply(rnk1, 2, base::rank, na.last = "keep", ties.method="random")
    }
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

run.gsva <- function(X, gmt, tau=1) {
  gsvapar <- GSVA::gsvaParam(X, gmt, tau=tau, maxDiff=TRUE)
  GSVA::gsva(gsvapar, verbose = FALSE)
}

run.ssgsea <- function(X, gmt, alpha=0.25) {
  ssgsea.par = GSVA::ssgseaParam(X, gmt, alpha=alpha, normalize=TRUE)
  GSVA::gsva(ssgsea.par, verbose = FALSE)
}

run.blitzgsea <- function(fc, gmt) {
  reticulate::py_require("blitzgsea")
  blitz <- try(reticulate::import("blitzgsea"))
  if(inherits(blitz,"try-error")) {
    message("Error: blitzgsea not installed")
    return(NULL)
  }
  df <- data.frame(names(fc),fc)
  res <- blitz$gsea(df, gmt)
  return(res[,1:6])
}

##----------------------------------------------------------------
##-------------------- LIMMA WITH METHODS ------------------------
##----------------------------------------------------------------

#' Differential expression analysis using limma
#'
#' @title Differential expression analysis using limma
#'
#' @description Performs differential expression analysis on a gene expression matrix using limma.
#' Handles single sample case by duplicating sample. Auto-detects reference level.
#'
#' @param X Gene expression matrix with genes in rows and samples in columns.
#' @param pheno Phenotype vector or factor for samples.
#' @param B Optional batch covariate matrix.
#' @param remove.na Logical for removing samples with missing phenotype.
#' @param fdr FDR threshold for identifying differentially expressed genes.
#' @param compute.means Logical for computing group means.
#' @param lfc Log fold change cutoff.
#' @param max.na Max proportion of missing values allowed for a gene.
#' @param ref Character vector of possible reference levels.
#' @param trend Logical for fitting a trend model.
#' @param verbose Verbosity level.
#'
#' @details This function performs differential expression analysis on the gene expression matrix \code{X} using limma.
#' It handles filtering, model design matrices, and output formatting.
#' The phenotype vector \code{pheno} is used to define the linear model.
#' Batch effects can be accounted for by providing a \code{B} matrix.
#' The function auto-detects the reference level or the user can provide possible values in \code{ref}.
#' It returns a data frame containing the limma analysis results.
#'
#' @return Data frame with limma results.
#'
#' @export
gx.limma <- function(X, pheno, B = NULL, remove.na = TRUE,
                     fdr = 0.05, compute.means = TRUE, lfc = 0.20,
                     max.na = 0.20, sort.by = "FC",
                     ref = c(
                       "ctrl", "ctr", "control", "ct", "dmso", "nt", "0", "0h", "0hr",
                       "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                       "wt", "wildtype", "untreated", "normal", "false", "healthy"
                     ),
                     trend = FALSE, robust = FALSE, method = 1, verbose = 1) {
  if (sum(duplicated(rownames(X))) > 0) {
    cat("WARNING:: matrix has duplicated rownames\n")
  }

  if (!is.null(B) && NCOL(B) == 1) {
    B <- matrix(B, ncol = 1)
    rownames(B) <- rownames(pheno)
    colnames(B) <- "batch"
  }

  ## detect single sample case
  is.single <- (max(table(pheno), na.rm = TRUE) == 1)
  if (is.single) {
    cat("WARNING:: no replicates, duplicating samples...\n")
    X <- cbind(X, X)
    pheno <- c(pheno, pheno)
    if (!is.null(B)) B <- rbind(B, B)
  }

  ## filter probes and samples??
  ii <- which(rowMeans(is.na(X)) <= max.na)
  jj <- 1:ncol(X)
  if (remove.na && any(is.na(pheno))) {
    jj <- which(!is.na(pheno))
    if (verbose > 0) message(sum(is.na(pheno) > 0), "with missing phenotype\n")
  }
  X0 <- X[ii, jj, drop = FALSE]
  pheno0 <- as.character(pheno[jj])
  X0 <- X0[!(rownames(X0) %in% c(NA, "", "NA")), ]
  B0 <- NULL
  if (!is.null(B)) B0 <- B[jj, , drop = FALSE]

  if (verbose > 0) {
    cat("analyzing", ncol(X0), "samples\n")
    cat("table.pheno: ", table(pheno), "samples\n")
    cat("testing", nrow(X0), "features\n")
    cat("lfc = ", lfc, "\n")
    cat("fdr = ", fdr, "\n")
    cat("max.na = ", max.na, "\n")
    if (!is.null(B0)) cat("including", ncol(B0), "batch covariates\n")
  }

  ## auto-detect reference
  pheno.ref <- c()
  ref.detected <- FALSE
  ref <- toupper(ref)

  is.ref <- (toupper(pheno0) %in% toupper(ref))
  is.ref2 <- grepl(paste(paste0("^", ref), collapse = "|"), pheno0, ignore.case = TRUE)
  if (!any(is.ref) && !all(is.ref2)) {
    is.ref <- is.ref2
  }
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)
  ref.detected

  if (ref.detected) {
    pheno.ref <- unique(pheno0[is.ref])
    if (verbose > 0) cat("setting reference to y=", pheno.ref, "\n")
    levels <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
  } else {
    if (verbose > 0) cat("WARNING: could not auto-detect reference\n")
    levels <- as.character(sort(unique(pheno0)))
    if (verbose > 0) cat("setting reference to first class:", levels[1], "\n")
  }
  if (length(levels) != 2) {
    stop("gx.limma::fatal error:only two class comparisons. Please use gx.limmaF().")
    return(NULL)
  }

  ## setup model and perform LIMMA. See LIMMA userguide p41 ("Two groups").
  ## https://bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
  if (method == 1) {
    ## first method without  contrast matrix
    design <- cbind(1, pheno0 == levels[2])
    colnames(design) <- c("intercept", "main_vs_ref")
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    fit <- limma::eBayes(fit, trend = trend, robust = robust)
    top <- limma::topTable(fit, coef = "main_vs_ref", number = nrow(X0), sort.by = "none")
  } else {
    ## second possible method with explicit contrast matrix
    design <- cbind(1 * (pheno0 == levels[1]), 1 * (pheno0 == levels[2]))
    design <- as.matrix(design)
    colnames(design) <- c("ref", "main")
    if (!is.null(B0)) {
      if (verbose > 0) cat("augmenting design matrix with:", paste(colnames(B0)), "\n")
      sel <- which(colMeans(B0 == 1) < 1) ## take out any constant term
      design <- cbind(design, B0[, sel, drop = FALSE])
    }
    fit <- limma::lmFit(X0, design)
    contr.matrix <- limma::makeContrasts(main_vs_ref = main - ref, levels = design)
    fit2 <- limma::contrasts.fit(fit, contr.matrix)
    fit2 <- limma::eBayes(fit2, trend = trend, robust = robust)
    top <- limma::topTable(fit2, coef = "main_vs_ref", number = Inf, sort.by = "none")
  }

  ## give rownames
  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }
  top <- top[rownames(X0), ]

  ## only significant
  if (!is.null(fdr) && !is.null(lfc)) {
    top <- top[which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc), ]
    if (verbose > 0) cat("found", nrow(top), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")
  }

  if (compute.means && nrow(top) > 0) {
    avg <- t(apply(
      X0[rownames(top), ], 1,
      function(x) tapply(x, pheno0, mean, na.rm = TRUE)
    ))
    avg <- avg[, as.character(levels), drop = FALSE]
    colnames(avg) <- paste0("AveExpr.", colnames(avg))
    top <- cbind(top, avg)
  }
  top$B <- NULL

  if (is.single) {
    top$P.Value <- NA
    top$adj.P.Val <- NA
    top$t <- NA
  }

  if (sort.by == "FC") {
    ## reorder on fold change
    top <- top[order(abs(top$logFC), decreasing = TRUE), ]
  }
  if (sort.by == "p") {
    ## reorder on fold change
    top <- top[order(top$P.Value), ]
  }

  ## unlist???
  return(top)
}

gsva.limma <- function(X, y, gmt, method="gsva") {
  gsetX <- gset.gsva(X, gmt, method=method)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)     
  gx.limma(gsetX, y, fdr=1, lfc=0, sort.by='none')
}

sing.limma <- function(X, y, gmt) {
  gsetX <- gset.singscore(X, gmt, return.score = TRUE)
  gsetX <- normalize_medians(gsetX, ignore.zero=NULL)     
  gx.limma(gsetX, y, fdr=1, lfc=0, sort.by='none')
}

## run.avgFC <- function(fc, matG) {
##   mt <- matrix_ttest(fc, matG)
##   res <- cbind(fc = mt$diff, pv = mt$pvalue)
##   colnames(res) <- c("logFC","p.value")
##   res
## }

prepare.SCSE <- function(X, gmt, path='scse') {
  bin.exe <- "scorer_linux_x64.bin"
  if(!file.exists(file.path(path,bin.exe))) {
    stop("[prepare.SCSE] Error: could not find binary in path")
  }
  message("[prepare.SCSE] preparing scSE files at path = ",path)
  data.dir <- file.path(path,"data")
  sig.dir <- file.path(path,"databases/000")
  res.dir <- file.path(path,"results")    
  if(dir.exists(data.dir)) unlink(paste0(data.dir,"/*.tsv"),force=TRUE)
  if(dir.exists(sig.dir)) unlink(paste0(sig.dir,"/*"),force=TRUE)
  dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
  dir.create(sig.dir, showWarnings=FALSE, recursive=TRUE)
  dfX <- data.frame(colnames(X),as.matrix(t(X)))
  data.table::fwrite(dfX, file=file.path(data.dir,"matrix.csv"), sep='\t', quote=FALSE, col.names=TRUE)  
  i=1
  for(i in 1:length(gmt)) {
    write( gmt[[i]], file=file.path(sig.dir,names(gmt)[i]))
  }
  dir.create(res.dir, showWarnings=FALSE, recursive=TRUE)  
}

run.SCSE <- function(X, gmt, removeLog2=NULL, scoreMean=TRUE, path='scse') {
  #library(jsonlite)
  if(!file.exists(file.path(path,"scorer_linux_x64.bin"))) path <- NULL
  if(is.null(path)) stop("Error: missing binary")
  bin.exe <- "scorer_linux_x64.bin"
  message("[run.SCSE] path = ",path)
  message("[run.SCSE] bin.exe = ",bin.exe)  
  if(is.null(removeLog2)) {
    removeLog2 <- min(X,na.rm=TRUE)==0 && max(X,na.rm=TRUE) < 20
  }
  scoreMean <- ifelse(scoreMean,TRUE,FALSE)
  removeLog2 <- as.integer(removeLog2)
  message("[run.SCSE] removeLog2 = ",removeLog2)
  message("[run.SCSE] scoreMean = ",scoreMean)
  results.dir <- file.path(path,"results")
  if(dir.exists(results.dir)) unlink(paste0(results.dir,"/*tsv"),force=TRUE)  
  conf <- data.frame(removeLog2=1*removeLog2,server=1,DBserver="0",
                     scoreMean=scoreMean)
  json <- jsonlite::toJSON(conf, pretty=TRUE)
  json <- gsub("^\\[\n|\n\\]$","",json)
  json <- gsub("  \\{","{",json)
  json <- gsub("  \\}","}",json)
  json <- gsub("  \\\"","\\\"",json)      
  write(json, file=file.path(path,"conf.json"))
  cmd <- paste0("(cd ",path," && ./",bin.exe,")")
  message("[run.SCSE] cmd = ",cmd)
  system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE )
  fn <- file.path(results.dir,"000_matrix.csv")
  if(!file.exists(fn)) {
    message("Error: failed computation")    
    return(NULL)
  }
  gsetX <- data.table::fread(fn, sep='\t', header=TRUE, check.names=FALSE)
  gsetX <- as.matrix(gsetX, rownames=1)  
  gsetX <- t(gsetX)
  # note: some genesets can be missing if expression is totally zero
  ii <- match(names(gmt), rownames(gsetX))
  jj <- match(colnames(X), colnames(gsetX))  
  if(any(is.na(ii))) message("Warning: some missing genesets in matrix")
  if(any(is.na(jj))) message("Warning: some missing samples in matrix")  
  gsetX <- gsetX[ii,jj]
  gsetX[is.na(gsetX)] <- 0
  rownames(gsetX) <- names(gmt)
  colnames(gsetX) <- colnames(X)
  gsetX
}

read.scse_db <- function(db.dir) {
  ff <- dir(db.dir, recursive=TRUE, full.names=TRUE)
  gmt <- lapply(ff, function(f) scan(f, what='character',quiet=TRUE))
  names(gmt) <- gsub("[.]txt","",basename(ff))
  return(gmt)
}


#' https://gist.github.com/gaoce/39e0907146c752c127728ad74e123b33
gao.ssgsea <- function(X, gene_sets, alpha = 0.25, scale = TRUE, norm = FALSE, single = TRUE) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)

    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos

      ## IK: I think this is wrong for alpha=0
      # rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      rank_alpha  = R_col[gene_ranks]^alpha * indicator_pos

      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

      step_cdf_diff = step_cdf_pos - step_cdf_neg

      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes

      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })

  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))

  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}
