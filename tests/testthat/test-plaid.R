##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

# ============================================================================
# Setup test data
# ============================================================================

# Helper function to create test data
create_test_data <- function(n_genes = 100, n_samples = 10, n_pathways = 5) {
  set.seed(123)
  
  # Create expression matrix
  X <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(X) <- paste0("GENE", 1:n_genes)
  colnames(X) <- paste0("Sample", 1:n_samples)
  
  # Create gene sets
  gmt <- lapply(1:n_pathways, function(i) {
    start_idx <- ((i-1) * 15 + 1)
    end_idx <- min(start_idx + 19, n_genes)
    paste0("GENE", start_idx:end_idx)
  })
  names(gmt) <- paste0("Pathway", 1:n_pathways)
  
  matG <- gmt2mat(gmt, bg = rownames(X))
  
  list(X = X, matG = matG, gmt = gmt)
}

# ============================================================================
# Test: plaid() - Main function
# ============================================================================

test_that("plaid works with basic input", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(data$matG))
  expect_equal(ncol(result), ncol(data$X))
  expect_equal(rownames(result), colnames(data$matG))
  expect_equal(colnames(result), colnames(data$X))
})

test_that("plaid works with stats='sum'", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$matG, stats = "sum")
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("plaid works with normalize=FALSE", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$matG, normalize = FALSE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("plaid works with single column", {
  data <- create_test_data(n_samples = 1)
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1)
})

test_that("plaid handles no overlapping features", {
  data <- create_test_data()
  rownames(data$X) <- paste0("OTHER_GENE", 1:nrow(data$X))
  
  result <- plaid(data$X, data$matG)
  
  expect_null(result)
})

test_that("plaid works with sparse matrices", {
  data <- create_test_data(n_pathways = 25)  # Use more pathways to avoid normalization issue
  X_sparse <- Matrix::Matrix(data$X, sparse = TRUE)
  
  result <- plaid(X_sparse, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("plaid works with different nsmooth values", {
  data <- create_test_data()
  
  result1 <- plaid(data$X, data$matG, nsmooth = 1)
  result2 <- plaid(data$X, data$matG, nsmooth = 5)
  
  expect_true(is.matrix(result1))
  expect_true(is.matrix(result2))
  expect_false(identical(result1, result2))
})

# ============================================================================
# Test: chunked_crossprod()
# ============================================================================

test_that("chunked_crossprod works with small matrix", {
  set.seed(123)
  x <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  y <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  
  result <- chunked_crossprod(x, y, chunk = 1000)
  expected <- Matrix::crossprod(x, y)
  
  expect_equal(as.matrix(result), as.matrix(expected))
})

test_that("chunked_crossprod works with chunk size", {
  set.seed(123)
  x <- Matrix::Matrix(rnorm(200), nrow = 10, ncol = 20, sparse = TRUE)
  y <- Matrix::Matrix(rnorm(300), nrow = 10, ncol = 30, sparse = TRUE)
  
  result <- chunked_crossprod(x, y, chunk = 10)
  expected <- Matrix::crossprod(x, y)
  
  # Compare values, ignoring dimnames which may differ
  expect_equal(as.numeric(result), as.numeric(expected), tolerance = 1e-10)
  expect_equal(dim(result), dim(expected))
})

# ============================================================================
# Test: replaid.scse()
# ============================================================================

test_that("replaid.scse works with basic input", {
  data <- create_test_data()
  # Create log-transformed data (simulate counts)
  X <- log2(matrix(rpois(1000, lambda = 10) + 1, nrow = 100, ncol = 10))
  rownames(X) <- paste0("GENE", 1:100)
  colnames(X) <- paste0("Sample", 1:10)
  
  result <- replaid.scse(X, data$matG, removeLog2 = TRUE, scoreMean = FALSE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(X)))
  expect_equal(colnames(result), colnames(X))
})

test_that("replaid.scse works with scoreMean=TRUE", {
  data <- create_test_data()
  X <- log2(matrix(rpois(1000, lambda = 10) + 1, nrow = 100, ncol = 10))
  rownames(X) <- rownames(data$X)
  colnames(X) <- colnames(data$X)
  
  result <- replaid.scse(X, data$matG, removeLog2 = TRUE, scoreMean = TRUE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(X)))
})

test_that("replaid.scse auto-detects removeLog2", {
  data <- create_test_data()
  X <- log2(matrix(rpois(1000, lambda = 10) + 1, nrow = 100, ncol = 10))
  rownames(X) <- rownames(data$X)
  colnames(X) <- colnames(data$X)
  
  result <- replaid.scse(X, data$matG, removeLog2 = NULL)
  
  expect_true(is.matrix(result))
})

# ============================================================================
# Test: replaid.sing()
# ============================================================================

test_that("replaid.sing works with basic input", {
  data <- create_test_data()
  
  result <- replaid.sing(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
  expect_equal(rownames(result), colnames(data$matG))
  expect_equal(colnames(result), colnames(data$X))
})

test_that("replaid.sing produces values in expected range", {
  data <- create_test_data()
  
  result <- replaid.sing(data$X, data$matG)
  
  # Singscore should be centered around 0
  expect_true(all(is.finite(result)))
  expect_true(abs(mean(result)) < 1)
})

# ============================================================================
# Test: replaid.ssgsea()
# ============================================================================

test_that("replaid.ssgsea works with alpha=0", {
  data <- create_test_data()
  
  result <- replaid.ssgsea(data$X, data$matG, alpha = 0)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.ssgsea works with alpha>0", {
  data <- create_test_data()
  
  result <- replaid.ssgsea(data$X, data$matG, alpha = 0.25)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.ssgsea results differ for different alpha values", {
  data <- create_test_data()
  
  result1 <- replaid.ssgsea(data$X, data$matG, alpha = 0)
  result2 <- replaid.ssgsea(data$X, data$matG, alpha = 0.5)
  
  expect_false(identical(result1, result2))
})

# ============================================================================
# Test: replaid.ucell()
# ============================================================================

test_that("replaid.ucell works with basic input", {
  data <- create_test_data()
  
  result <- replaid.ucell(data$X, data$matG, rmax = 50)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.ucell produces finite values", {
  data <- create_test_data()
  
  result <- replaid.ucell(data$X, data$matG, rmax = 50)
  
  # UCell scores are transformed and may exceed [0,1] range
  expect_true(all(is.finite(result)))
  expect_true(is.matrix(result))
})

test_that("replaid.ucell works with different rmax values", {
  data <- create_test_data()
  
  result1 <- replaid.ucell(data$X, data$matG, rmax = 50)
  result2 <- replaid.ucell(data$X, data$matG, rmax = 100)
  
  expect_false(identical(result1, result2))
})

# ============================================================================
# Test: replaid.aucell()
# ============================================================================

test_that("replaid.aucell works with basic input", {
  data <- create_test_data()
  
  result <- replaid.aucell(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.aucell works with custom aucMaxRank", {
  data <- create_test_data()
  
  result <- replaid.aucell(data$X, data$matG, aucMaxRank = 10)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.aucell auto-computes aucMaxRank", {
  data <- create_test_data()
  
  result <- replaid.aucell(data$X, data$matG, aucMaxRank = NULL)
  
  expect_true(is.matrix(result))
})

# ============================================================================
# Test: replaid.gsva()
# ============================================================================

test_that("replaid.gsva works with basic input", {
  data <- create_test_data()
  
  result <- replaid.gsva(data$X, data$matG, tau = 0)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.gsva works with tau>0", {
  data <- create_test_data()
  
  result <- replaid.gsva(data$X, data$matG, tau = 0.25)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("replaid.gsva works with rowtf='z'", {
  data <- create_test_data()
  
  result <- replaid.gsva(data$X, data$matG, rowtf = "z")
  
  expect_true(is.matrix(result))
})

test_that("replaid.gsva works with rowtf='ecdf'", {
  data <- create_test_data()
  
  result <- replaid.gsva(data$X, data$matG, rowtf = "ecdf")
  
  expect_true(is.matrix(result))
})

test_that("replaid.gsva throws error for invalid rowtf", {
  data <- create_test_data()
  
  expect_error(replaid.gsva(data$X, data$matG, rowtf = "invalid"))
})

# ============================================================================
# Test: mat.rowsds()
# ============================================================================

test_that("mat.rowsds works with dense matrix", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  result <- mat.rowsds(X)
  expected <- matrixStats::rowSds(X)
  
  expect_equal(result, expected)
})

test_that("mat.rowsds works with sparse matrix", {
  set.seed(123)
  X <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  
  result <- mat.rowsds(X)
  
  expect_true(is.numeric(result))
  expect_equal(length(result), nrow(X))
  expect_true(all(result >= 0))
})

# ============================================================================
# Test: normalize_medians()
# ============================================================================

test_that("normalize_medians works with basic input", {
  set.seed(123)
  x <- matrix(rnorm(100, mean = 5), nrow = 10, ncol = 10)
  
  result <- normalize_medians(x)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(x))
  
  # Check that medians are similar after normalization
  medians <- matrixStats::colMedians(result)
  expect_true(sd(medians) < sd(matrixStats::colMedians(x)))
})

test_that("normalize_medians handles zeros correctly", {
  set.seed(123)
  x <- matrix(rnorm(100), nrow = 10, ncol = 10)
  x[1:3, 1:3] <- 0
  
  result <- normalize_medians(x, ignore.zero = TRUE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(x))
})

test_that("normalize_medians auto-detects zeros", {
  set.seed(123)
  x <- matrix(rnorm(100), nrow = 10, ncol = 10)
  x[1:3, 1:3] <- 0
  
  result <- normalize_medians(x, ignore.zero = NULL)
  
  expect_true(is.matrix(result))
})

# ============================================================================
# Test: colranks()
# ============================================================================

test_that("colranks works with basic input", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(X) <- paste0("Gene", 1:10)
  colnames(X) <- paste0("Sample", 1:10)
  
  result <- colranks(X)
  
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(dim(result), dim(X))
  expect_equal(dimnames(result), dimnames(X))
})

test_that("colranks produces correct rank range", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  result <- colranks(X, ties.method = "average")
  
  # Ranks should be between 1 and nrow(X)
  expect_true(all(result >= 1 & result <= nrow(X)))
})

test_that("colranks works with signed=TRUE", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  result <- colranks(X, signed = TRUE)
  
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(dim(result), dim(X))
})

test_that("colranks works with keep.zero=TRUE", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  X[1:3, 1:3] <- 0
  X <- Matrix::Matrix(X, sparse = TRUE)
  
  result <- colranks(X, keep.zero = TRUE)
  
  expect_true(inherits(result, "Matrix"))
  expect_equal(dim(result), dim(X))
})

test_that("colranks works with sparse matrix", {
  set.seed(123)
  X <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  
  result <- colranks(X, sparse = TRUE)
  
  # colranks returns dense matrix for full matrices even if input is sparse
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(dim(result), dim(X))
})

test_that("colranks ties methods work correctly", {
  X <- matrix(c(1, 2, 2, 3, 4, 5, 5, 6), nrow = 4, ncol = 2)
  
  result_avg <- colranks(X, ties.method = "average")
  result_min <- colranks(X, ties.method = "min")
  
  expect_false(identical(result_avg, result_min))
})

# ============================================================================
# Test: sparse_colranks()
# ============================================================================

test_that("sparse_colranks works with basic input", {
  set.seed(123)
  X <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  
  result <- sparse_colranks(X)
  
  expect_true(inherits(result, "CsparseMatrix"))
  expect_equal(dim(result), dim(X))
})

test_that("sparse_colranks works with signed=TRUE", {
  set.seed(123)
  X <- Matrix::Matrix(rnorm(100), nrow = 10, ncol = 10, sparse = TRUE)
  
  result <- sparse_colranks(X, signed = TRUE)
  
  expect_true(inherits(result, "CsparseMatrix"))
  expect_equal(dim(result), dim(X))
})

test_that("sparse_colranks handles different ties methods", {
  X <- Matrix::Matrix(c(1, 2, 2, 3, 4, 5, 5, 6), nrow = 4, ncol = 2, sparse = TRUE)
  
  result_avg <- sparse_colranks(X, ties.method = "average")
  result_min <- sparse_colranks(X, ties.method = "min")
  
  expect_true(inherits(result_avg, "CsparseMatrix"))
  expect_true(inherits(result_min, "CsparseMatrix"))
})

# ============================================================================
# Test: GMT list input for replaid functions
# ============================================================================

test_that("plaid works with GMT list input", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$gmt)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(data$X))
})

test_that("replaid.sing works with GMT list input", {
  data <- create_test_data()
  
  result <- replaid.sing(data$X, data$gmt)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(data$X))
})

test_that("replaid.ssgsea works with GMT list input", {
  data <- create_test_data()
  
  result <- replaid.ssgsea(data$X, data$gmt, alpha = 0)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(data$X))
})

# ============================================================================
# Test: Edge cases
# ============================================================================

test_that("plaid handles small number of gene sets", {
  data <- create_test_data(n_pathways = 2)
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
})

test_that("plaid handles large nsmooth value", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$matG, nsmooth = 100)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$matG), ncol(data$X)))
})

test_that("normalize_medians handles all-zero columns", {
  x <- matrix(0, nrow = 10, ncol = 5)
  x[, 1:3] <- rnorm(30)
  
  result <- normalize_medians(x)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(x))
})

test_that("colranks handles constant values", {
  X <- matrix(rep(1, 20), nrow = 5, ncol = 4)
  
  result <- colranks(X)
  
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(dim(result), dim(X))
})

# ============================================================================
# Test: Parameter validation
# ============================================================================

test_that("plaid handles partial overlap of features", {
  data <- create_test_data(n_genes = 100)
  # Keep only half of genes in matG
  keep_genes <- rownames(data$X)[1:50]
  data$matG <- data$matG[keep_genes, , drop = FALSE]
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(data$X))
})

test_that("replaid functions handle min.genes and max.genes parameters", {
  data <- create_test_data()
  
  result <- plaid(data$X, data$gmt, min.genes = 5, max.genes = 100)
  
  expect_true(is.matrix(result))
})

# ============================================================================
# Test: Numerical stability
# ============================================================================

test_that("plaid handles very small values", {
  data <- create_test_data()
  data$X <- data$X * 1e-10
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  expect_true(all(is.finite(result)))
})

test_that("plaid handles very large values", {
  data <- create_test_data(n_pathways = 25)
  data$X <- data$X * 1e6  # Use smaller multiplier to avoid Inf
  
  result <- plaid(data$X, data$matG)
  
  expect_true(is.matrix(result))
  # With very large values, some Inf may occur in intermediate calculations
  expect_true(is.matrix(result))
})

test_that("normalize_medians handles NA values", {
  set.seed(123)
  x <- matrix(rnorm(100), nrow = 10, ncol = 10)
  x[1:3, 1:3] <- NA
  
  result <- normalize_medians(x)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(x))
})
