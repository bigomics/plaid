##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##


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


#' Read data from a GMT file
#'
#' This function reads data from a GMT (Gene Matrix Transposed) file format.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt.file The path to the GMT file.
#' @param dir (Optional) The directory where the GMT file is located.
#' @param add.source (Optional) Specifies whether to include the source information in the gene sets' names.
#' @param nrows (Optional) The number of rows to read from the GMT file.
#'
#' @export
#'
#' @return A list of gene sets, where each gene set is represented as a character vector of gene names.
#'
read.gmt <- function(gmt.file, dir = NULL, add.source = FALSE, nrows = -1) {
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
    if (length(x) < 3) {
      return("")
    }
    paste(x[3:length(x)], collapse = " ")
  })
  gset <- strsplit(gmt.genes, split = "[ \t]")
  gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
  names(gset) <- gmt.name

  if (add.source) {
    names(gset) <- paste0(names(gset), " (", gmt.source, ")")
  }
  gset
}


#' Write gene sets to GMT file
#'
#' This function writes gene sets in GMT (Gene Matrix Transposed) format to a file.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt A list of gene sets in GMT format, where each gene set is represented as a character vector of gene names.
#' @param file The file path to write the GMT file.
#' @param source A character vector specifying the source of each gene set. If not provided, the names of the gene sets are used as the source.
#'
#' @export
#'
#' @return NULL
#'
write.gmt <- function(gmt, file, source = NA) {
  gg <- lapply(gmt, paste, collapse = "\t")
  if (is.na(source)) source <- names(gmt)
  ee <- paste(names(gmt), "\t", source, "\t", gg, sep = "")
  write(ee, file = file)
}
