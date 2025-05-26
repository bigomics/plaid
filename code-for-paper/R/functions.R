##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##
##

## These are not striclty needed to run plaid but only for running
## other methods for comparison. No need to declare dependency on
## these packages.
##
##


##----------------------------------------------------------------
##-------------------- RUN OTHER METHODS -------------------------
##----------------------------------------------------------------


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

run.avgFC <- function(fc, matG) {
  mt <- matrix_ttest(fc, matG)
  res <- cbind(fc = mt$diff, pv = mt$pvalue)
  colnames(res) <- c("logFC","p.value")
  res
}

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
    ssgsea <- gset.gsva(X, gmt, method="ssgsea"),
    ucell  <- t(UCell::ScoreSignatures_UCell(X, gmt)),  ## needs logx
    aucell <- AUCell::getAUC(AUCell::AUCell_run(X, gmt)), ## uses rank
    ## scse.sum  <- run.SCSE(X, gmt, removeLog2=TRUE, scoreMean=FALSE),
    scse <- run.SCSE(X, gmt, removeLog2=FALSE, scoreMean=TRUE, path="../scse"),
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
  system(paste("(cd",path,"&&",bin.exe,")"), ignore.stdout=TRUE, ignore.stderr=TRUE )
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
  if(any(is.na(ii))) message("Warning: missing genesets in matrix")
  if(any(is.na(jj))) message("Warning: missing samples in matrix")  
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
