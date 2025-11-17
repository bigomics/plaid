##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom BiocSet es_elementset
NULL

#' Extract expression matrix from Bioconductor objects
#'
#' @description Internal function to extract expression matrix from various
#' Bioconductor objects including SummarizedExperiment and SingleCellExperiment.
#' 
#' @param object A Bioconductor object (SummarizedExperiment, SingleCellExperiment)
#'   or a regular matrix.
#' @param assay Character or integer specifying which assay to extract. Default "logcounts".
#'   For SummarizedExperiment objects, common options are "counts", "logcounts", 
#'   "normcounts". The function will try several common names if the specified
#'   assay is not found.
#' @param log.transform Logical. If TRUE and the data appears to be counts (not log-transformed),
#'   will apply log2(x+1) transformation. Default FALSE.
#'
#' @return A matrix with genes/features on rows and samples on columns.
#'
#' @details This function is designed to handle various Bioconductor object types:
#' \itemize{
#'   \item SummarizedExperiment: extracts from assays
#'   \item SingleCellExperiment: extracts from assays (prefers logcounts)
#'   \item Matrix/matrix: returns as-is
#' }
#'
#' The function will attempt to find an appropriate assay by trying common names
#' in the following order: specified assay, "logcounts", "normcounts", "counts",
#' or the first available assay.
#'
#' @noRd
.extract_expression_matrix <- function(object, assay = "logcounts", log.transform = FALSE) {
  
  # If already a matrix, return as-is
  if (is.matrix(object) || inherits(object, "Matrix")) {
    return(object)
  }
  
  # Check for SummarizedExperiment or SingleCellExperiment
  if (inherits(object, "SummarizedExperiment") || 
      inherits(object, "SingleCellExperiment")) {
    
    # Try to extract the requested assay
    avail_assays <- SummarizedExperiment::assayNames(object)
    
    if (length(avail_assays) == 0) {
      stop("[.extract_expression_matrix] No assays found in object.")
    }
    
    # Try to find the best assay
    if (assay %in% avail_assays) {
      X <- SummarizedExperiment::assay(object, assay)
    } else {
      # Try common alternatives in order of preference
      preferred_order <- c("logcounts", "normcounts", "counts")
      found_assay <- NULL
      
      for (pref in preferred_order) {
        if (pref %in% avail_assays) {
          found_assay <- pref
          break
        }
      }
      
      if (is.null(found_assay)) {
        # Just use first available assay
        found_assay <- avail_assays[1]
        message("[.extract_expression_matrix] Using first available assay: ", found_assay)
      } else {
        message("[.extract_expression_matrix] Requested assay '", assay, 
                "' not found. Using: ", found_assay)
      }
      
      X <- SummarizedExperiment::assay(object, found_assay)
    }
    
    # Apply log transformation if requested
    if (log.transform) {
      # Check if data looks like counts (non-negative integers or large values)
      if (min(X, na.rm = TRUE) >= 0 && max(X, na.rm = TRUE) > 100) {
        message("[.extract_expression_matrix] Applying log2(x+1) transformation.")
        if (inherits(X, "dgCMatrix")) {
          X@x <- log2(X@x + 1)
        } else {
          X <- log2(X + 1)
        }
      }
    }
    
    return(X)
  }
  
  # If we get here, unsupported object type
  stop("[.extract_expression_matrix] Unsupported object type: ", class(object)[1],
       "\nSupported types: matrix, Matrix, SummarizedExperiment, SingleCellExperiment")
}


#' Convert BiocSet to sparse matrix format for plaid
#'
#' @description Internal function to convert BiocSet objects to the sparse
#' matrix format required by plaid (genes on rows, gene sets on columns).
#'
#' @param geneset A BiocSet object, list (GMT format), or already a matrix.
#' @param background Character vector of background genes to include in the matrix.
#'   Default NULL uses all genes from the gene sets.
#' @param min.genes Minimum number of genes required for a gene set to be included.
#'   Default 5.
#' @param max.genes Maximum number of genes allowed for a gene set to be included.
#'   Default 500.
#'
#' @return A sparse matrix (dgCMatrix) with genes on rows and gene sets on columns.
#'
#' @details This function handles conversion from BiocSet objects to the
#' sparse matrix format used by plaid. It can also handle GMT lists or matrices
#' directly. Gene sets are filtered by size (min/max genes).
#'
#' For BiocSet objects, the function extracts element-set mappings and converts
#' them to a binary sparse matrix indicating gene membership.
#'
#' @noRd
.convert_geneset_to_matrix <- function(geneset, background = NULL, 
                                       min.genes = 5, max.genes = 500) {
  
  # If already a matrix, return as-is
  if (is.matrix(geneset) || inherits(geneset, "Matrix")) {
    return(geneset)
  }
  
  # Handle BiocSet objects
  if (inherits(geneset, "BiocSet")) {
    message("[.convert_geneset_to_matrix] Converting BiocSet to matrix format...")
    
    # Extract the element-set mapping from BiocSet
    # BiocSet uses a tibble structure with 'element' and 'set' columns
    es_tbl <- BiocSet::es_elementset(geneset)
    
    if (nrow(es_tbl) == 0) {
      stop("[.convert_geneset_to_matrix] BiocSet object is empty.")
    }
    
    # Convert to GMT list format first
    gmt <- split(es_tbl$element, es_tbl$set)
    gmt <- lapply(gmt, as.character)
    
  } else if (is.list(geneset)) {
    # Already in GMT format
    gmt <- geneset
  } else {
    stop("[.convert_geneset_to_matrix] Unsupported geneset type: ", class(geneset)[1],
         "\nSupported types: BiocSet, list (GMT), matrix, Matrix")
  }
  
  # Filter gene sets by size
  gset_sizes <- vapply(gmt, length, integer(1))
  valid_sets <- gset_sizes >= min.genes & gset_sizes <= max.genes
  
  if (sum(valid_sets) == 0) {
    stop("[.convert_geneset_to_matrix] No gene sets passed size filters (min=", 
         min.genes, ", max=", max.genes, ")")
  }
  
  if (sum(!valid_sets) > 0) {
    message("[.convert_geneset_to_matrix] Filtered out ", sum(!valid_sets), 
            " gene sets (size filters: ", min.genes, "-", max.genes, " genes)")
  }
  
  gmt <- gmt[valid_sets]
  
  # Convert to sparse matrix using existing gmt2mat function
  matG <- gmt2mat(gmt, bg = background, sparse = TRUE)
  
  return(matG)
}
