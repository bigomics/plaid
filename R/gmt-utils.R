## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.

#' @importFrom utils head read.csv
#' @importFrom Matrix Matrix rowSums which
#' @importFrom collapse fduplicated fmatch fsubset funique qtable vec vlengths whichNA whichv
NULL

#' Convert GMT to Binary Matrix
#'
#' @description Convert a GMT file (Gene Matrix Transposed) to a binary matrix,
#' where rows represent genes and columns represent gene sets.
#' The binary matrix indicates presence or absence of genes in a gene set.
#'
#' @param gmt List representing the GMT file: each element is a character vector representing a gene set.
#' @param max.genes Max number of genes to include in the binary matrix. Default = -1 to include all genes.
#' @param ntop Number of top genes to consider for each gene set. Default = -1 to include all genes.
#' @param sparse Logical: create a sparse matrix. Default `TRUE`. If `FALSE` creates a dense matrix.
#' @param bg Character vector of background gene set. Default `NULL` to consider all unique genes.
#' @param use.multicore Logical: use parallel processing ('parallel' R package). Default `FALSE`. Deprecated.
#'
#' @export
#'
#' @return A binary matrix representing the presence or absence of genes in each gene set.
#'         Rows correspond to genes, and columns correspond to gene sets.
#'
#' @examples
#' # Create example GMT data
#' gmt <- list(
#'   "Pathway1" = c("GENE1", "GENE2", "GENE3"),
#'   "Pathway2" = c("GENE2", "GENE4", "GENE5"),
#'   "Pathway3" = c("GENE1", "GENE5", "GENE6")
#' )
#' 
#' # Convert to binary matrix
#' mat <- gmt2mat(gmt)
#' print(mat)
#' 
#' # Create dense matrix instead of sparse
#' mat_dense <- gmt2mat(gmt, sparse = FALSE)
#' print(mat_dense)
gmt2mat <- function(gmt,
                    max.genes = -1,
                    ntop = -1,
                    sparse = TRUE,
                    bg = NULL,
                    use.multicore = FALSE) {

  if (use.multicore) {
    warning(
      "`use.multicore` is deprecated and will be removed in a future version."
    )
  }

  gmt <- gmt[order(vlengths(gmt), decreasing = TRUE)]
  gmt <- gmt[!duplicated(names(gmt))]
  if (ntop > 0) gmt <- lapply(gmt, utils::head, n = ntop)
  if (is.null(names(gmt))) names(gmt) <- paste0("gmt.", seq_along(gmt))

  genes <- vec(gmt)
  if (is.null(bg)) {
    bg <- names(sort(qtable(genes), decreasing = TRUE))
  } else {
    bg <- funique(bg)
  }
  
  if (max.genes < 0) max.genes <- length(bg)
  gg <- utils::head(bg, n = max.genes)
  kk <- names(gmt)

  NR <- length(gg)
  NC <- length(kk)
  idx_row <- fmatch(genes, gg)
  idx_col <- rep.int(seq_len(NC), vlengths(gmt))

  idx_not_NA <- whichNA(idx_row, invert = TRUE)
  if (length(idx_not_NA) != length(idx_row)) {
    idx_row <- fsubset(idx_row, idx_not_NA)
    idx_col <- fsubset(idx_col, idx_not_NA)
  }

  # Matrices are in column-major order
  idx <- NR * (idx_col - 1L) + idx_row

  # Identify duplicates
  idx_not_dup <- whichv(fduplicated(idx, all = FALSE), FALSE)
  if (length(idx_not_dup) != length(idx)) {
    idx_row <- fsubset(idx_row, idx_not_dup)
    idx_col <- fsubset(idx_col, idx_not_dup)
    idx <- fsubset(idx, idx_not_dup)
  }

  if (sparse) {
    D <- Matrix::sparseMatrix(
      i = idx_row,
      j = idx_col,
      x = 1,
      dims = c(NR, NC),
      check = FALSE,
      use.last.ij = FALSE
    )
  } else {
    D <- matrix(0, nrow = NR, ncol = NC)
    D[idx] <- 1
  }
  rownames(D) <- gg
  colnames(D) <- kk
  D <- D[order(-Matrix::rowSums(D != 0, na.rm = TRUE)), ,drop=FALSE]

  return(D)

}

#' Convert Binary Matrix to GMT
#'
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
#' @examples
#' # Create example binary matrix
#' mat <- matrix(0, nrow = 6, ncol = 3)
#' rownames(mat) <- paste0("GENE", 1:6)
#' colnames(mat) <- paste0("Pathway", 1:3)
#' mat[1:3, 1] <- 1  # Pathway1: GENE1, GENE2, GENE3
#' mat[c(2,4,5), 2] <- 1  # Pathway2: GENE2, GENE4, GENE5
#' mat[c(1,5,6), 3] <- 1  # Pathway3: GENE1, GENE5, GENE6
#' 
#' # Convert to GMT list
#' gmt <- mat2gmt(mat)
#' print(gmt)
#'
mat2gmt <- function(mat) {
  idx <- Matrix::which(mat != 0, arr.ind = TRUE)
  gmt <- tapply(rownames(idx), idx[, 2], list)
  names(gmt) <- colnames(mat)[as.integer(names(gmt))]
  return(gmt)
}


#' Read GMT File
#'
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
#' @examples
#' \donttest{
#' # Read GMT file (requires file to exist)
#' gmt_file <- system.file("extdata", "hallmarks.gmt", package = "plaid")
#' if (file.exists(gmt_file)) {
#'   gmt <- read.gmt(gmt_file)
#'   print(names(gmt))
#'   print(head(gmt[[1]]))
#'   
#'   # Read with source information
#'   gmt_with_source <- read.gmt(gmt_file, add.source = TRUE)
#'   print(head(names(gmt_with_source)))
#' }
#' }
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
  gmt <- lapply(gmt, function(s) strsplit(s, split = "\t")[[1]])
  names(gmt) <- NULL
  gmt.name <- vapply(gmt, "[", character(1), 1)
  gmt.source <- vapply(gmt, "[", character(1), 2)
  gmt.genes <- vapply(gmt, function(x) {
    if (length(x) < 3) return("");
    paste(x[3:length(x)], collapse = " ")
  }, character(1))
  gset <- strsplit(gmt.genes, split = "[ \t]")
  gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
  names(gset) <- gmt.name

  if (add.source)
    names(gset) <- paste0(names(gset), " (", gmt.source, ")")
  
  return(gset)

}


#' Write GMT File
#'
#' @description Write gene sets to GMT file (Gene Matrix Transposed).
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt A list of gene sets in GMT format: each gene set is represented as a vector of gene names.
#' @param file The file path to write the GMT file.
#' @param source A character vector specifying the source of each gene set.
#' If not provided, the names of the gene sets are used as the source.
#' 
#' @export
#' @return Does not return anything.
#' 
#' @examples
#' # Create example GMT data
#' gmt <- list(
#'   "Pathway1" = c("GENE1", "GENE2", "GENE3"),
#'   "Pathway2" = c("GENE2", "GENE4", "GENE5"),
#'   "Pathway3" = c("GENE1", "GENE5", "GENE6")
#' )
#' 
#' \donttest{
#' # Write to GMT file (creates file in temp directory)
#' temp_file <- tempfile(fileext = ".gmt")
#' write.gmt(gmt, temp_file)
#' 
#' # Write with custom source information
#' temp_file2 <- tempfile(fileext = ".gmt")
#' write.gmt(gmt, temp_file2, source = c("DB1", "DB2", "DB3"))
#' 
#' # Clean up
#' unlink(c(temp_file, temp_file2))
#' }
#'
write.gmt <- function(gmt, file, source = NA) {
  gg <- lapply(gmt, paste, collapse = "\t")
  if (length(source) == 1 && is.na(source[1])) source <- names(gmt)
  ee <- paste(names(gmt), "\t", source, "\t", gg, sep = "")
  write(ee, file = file)
}
