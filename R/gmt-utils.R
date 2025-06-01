## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.

#' @description Convert a GMT file (Gene Matrix Transposed) to a binary matrix,
#' where rows represent genes and columns represent gene sets.
#' The binary matrix indicates presence or absence of genes in a gene set.
#'
#' @param gmt List representing the GMT file: each element is a character vector representing a gene set.
#' @param max.genes Max number of genes to include in the binary matrix. Default = -1 to include all genes.
#' @param ntop Number of top genes to consider for each gene set. Default = -1 to include all genes.
#' @param sparse Logical: create a sparse matrix. Default `TRUE`. If `FALSE` creates a dense matrix.
#' @param bg Character vector of background gene set. Default `NULL` to consider all unique genes.
#' @param use.multicore Logical: use parallel processing ('parallel' R package). Default `TRUE`.
#'
#' @export
#'
#' @return A binary matrix representing the presence or absence of genes in each gene set.
#'         Rows correspond to genes, and columns correspond to gene sets.
gmt2mat <- function(gmt,
                    max.genes = -1,
                    ntop = -1, sparse = TRUE,
                    bg = NULL,
                    use.multicore = TRUE) {

  gmt <- gmt[order(-sapply(gmt, length))]
  gmt <- gmt[!duplicated(names(gmt))]
  if (ntop > 0) gmt <- lapply(gmt, utils::head, n = ntop)
  
  if (is.null(names(gmt))) names(gmt) <- paste("gmt.", 1:length(gmt), sep = "")
  if (is.null(bg))
    bg <- names(sort(table(unlist(gmt)), decreasing = TRUE))
  
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

  return(D)

}

#' @description Convert binary matrix to a GMT (Gene Matrix Transposed) list.
#' The binary matrix indicates presence or absence of genes in each gene set.
#' Rows represent genes and columns represent gene sets.
#'
#' @param mat Matrix with non-zero entries representing genes in each gene set.
#' Rows represent genes and columns represent gene sets.
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
  return(gmt)
}


#' @description Read data from a GMT file (Gene Matrix Transposed).
#' The GMT format is commonly used to store gene sets or gene annotations.
#' @param gmt.file Path to GMT file.
#' @param dir (Optional) The directory where the GMT file is located.
#' @param add.source (optional) Include the source information in the gene sets' names.
#' @param nrows (optional) Number of rows to read from the GMT file.
#'
#' @export
#'
#' @return A list of gene sets: each gene set is represented as a character vector of gene names.
#' 
read.gmt <- function(gmt.file,
                     dir = NULL,
                     add.source = FALSE,
                     nrows = -1) {
  f0 <- gmt.file
  if (strtrim(gmt.file, 1) == "/") dir <- NULL
  if (!is.null(dir)) f0 <- paste(sub("/$", "", dir), "/", gmt.file, sep = "")
  gmt <- utils::read.csv(f0, sep = "!", header = FALSE, comment.char = "#", nrows = nrows)[, 1]
  gmt <- as.character(gmt)
  gmt <- sapply(gmt, strsplit, split = "\t")
  names(gmt) <- NULL
  gmt.name <- sapply(gmt, "[", 1)
  gmt.source <- sapply(gmt, "[", 2)
  gmt.genes <- sapply(gmt, function(x) {
    if (length(x) < 3) return("");
    paste(x[3:length(x)], collapse = " ")
  })
  gset <- strsplit(gmt.genes, split = "[ \t]")
  gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
  names(gset) <- gmt.name

  if (add.source)
    names(gset) <- paste0(names(gset), " (", gmt.source, ")")
  
  return(gset)

}


#' @description Write gene sets to GMT file (Gene Matrix Transposed).
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt A list of gene sets in GMT format: each gene set is represented as a vector of gene names.
#' @param file The file path to write the GMT file.
#' @param source A character vector specifying the source of each gene set.
#' If not provided, the names of the gene sets are used as the source.
#' 
#' @export
#' @return NULL
#' 
write.gmt <- function(gmt, file, source = NA) {
  gg <- lapply(gmt, paste, collapse = "\t")
  if (is.na(source)) source <- names(gmt)
  ee <- paste(names(gmt), "\t", source, "\t", gg, sep = "")
  write(ee, file = file)
}
