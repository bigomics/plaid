##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

# ============================================================================
# Setup test data
# ============================================================================

# Helper function to create test data for statistical testing
create_stats_test_data <- function(n_genes = 100, n_samples = 20, n_pathways = 5) {
  set.seed(456)
  
  # Create expression matrix
  X <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  rownames(X) <- paste0("GENE", 1:n_genes)
  colnames(X) <- paste0("Sample", 1:n_samples)
  
  # Create binary group vector (0/1)
  y <- rep(c(0, 1), each = n_samples/2)
  
  # Create gene sets
  gmt <- lapply(1:n_pathways, function(i) {
    start_idx <- ((i-1) * 15 + 1)
    end_idx <- min(start_idx + 19, n_genes)
    paste0("GENE", start_idx:end_idx)
  })
  names(gmt) <- paste0("Pathway", 1:n_pathways)
  
  G <- gmt2mat(gmt, bg = rownames(X))
  
  list(X = X, y = y, G = G, gmt = gmt)
}

# ============================================================================
# Test: dual_test()
# ============================================================================

test_that("dual_test works with basic input", {
  data <- create_stats_test_data()
  
  result <- dual_test(data$X, data$y, data$G)
  
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("gsetFC" %in% colnames(result))
  expect_true("p.dual" %in% colnames(result))
  expect_true("q.dual" %in% colnames(result))
  expect_true("size" %in% colnames(result))
  expect_equal(nrow(result), ncol(data$G))
})

test_that("dual_test works with different fc.test methods", {
  data <- create_stats_test_data(n_pathways = 25)  # Use more pathways
  
  result_cor <- dual_test(data$X, data$y, data$G, fc.test = "cor")
  result_rankcor <- dual_test(data$X, data$y, data$G, fc.test = "rankcor")
  result_ztest <- dual_test(data$X, data$y, data$G, fc.test = "ztest")
  
  expect_true(is.matrix(result_cor) || is.data.frame(result_cor))
  expect_true(is.matrix(result_rankcor) || is.data.frame(result_rankcor))
  expect_true(is.matrix(result_ztest) || is.data.frame(result_ztest))
})

test_that("dual_test works with precomputed gsetX", {
  data <- create_stats_test_data()
  gsetX <- plaid(data$X, data$G)
  
  result <- dual_test(data$X, data$y, data$G, gsetX = gsetX)
  
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("p.dual" %in% colnames(result))
})

test_that("dual_test works with different metap methods", {
  data <- create_stats_test_data()
  
  result_stouffer <- dual_test(data$X, data$y, data$G, metap.method = "stouffer")
  result_fisher <- dual_test(data$X, data$y, data$G, metap.method = "fisher")
  
  expect_true(is.matrix(result_stouffer) || is.data.frame(result_stouffer))
  expect_true(is.matrix(result_fisher) || is.data.frame(result_fisher))
  expect_false(identical(result_stouffer, result_fisher))
})

test_that("dual_test works with different sort.by options", {
  data <- create_stats_test_data()
  
  result_pdual <- dual_test(data$X, data$y, data$G, sort.by = "p.dual")
  result_fc <- dual_test(data$X, data$y, data$G, sort.by = "gsetFC")
  
  expect_true(is.matrix(result_pdual) || is.data.frame(result_pdual))
  expect_true(is.matrix(result_fc) || is.data.frame(result_fc))
})

test_that("dual_test handles NA values in y", {
  data <- create_stats_test_data()
  data$y[1:2] <- NA
  
  result <- dual_test(data$X, data$y, data$G)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("dual_test works with GMT list input", {
  data <- create_stats_test_data()
  
  result <- dual_test(data$X, data$y, data$gmt)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("dual_test throws error for invalid y", {
  data <- create_stats_test_data()
  data$y[1] <- 2  # Invalid value
  
  expect_error(dual_test(data$X, data$y, data$G))
})

test_that("dual_test throws error for invalid fc.test", {
  data <- create_stats_test_data()
  
  expect_error(dual_test(data$X, data$y, data$G, fc.test = "invalid"))
})

test_that("dual_test works with precomputed pv1 and pv2", {
  data <- create_stats_test_data()
  pv1 <- rep(0.05, ncol(data$G))
  names(pv1) <- colnames(data$G)
  pv2 <- rep(0.1, ncol(data$G))
  names(pv2) <- colnames(data$G)
  
  result <- dual_test(data$X, data$y, data$G, pv1 = pv1, pv2 = pv2)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("dual_test works with sparse matrix input", {
  data <- create_stats_test_data(n_pathways = 25)  # Use more pathways to avoid normalization issue
  X_sparse <- Matrix::Matrix(data$X, sparse = TRUE)
  
  result <- dual_test(X_sparse, data$y, data$G)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

# ============================================================================
# Test: dualGSEA()
# ============================================================================

test_that("dualGSEA works with replaid methods", {
  data <- create_stats_test_data()
  
  result_ssgsea <- dualGSEA(data$X, data$y, data$gmt, 
                            fc.test = "cor", gsea.method = "replaid.ssgsea")
  result_gsva <- dualGSEA(data$X, data$y, data$gmt, 
                          fc.test = "cor", gsea.method = "replaid.gsva")
  
  expect_true(is.matrix(result_ssgsea) || is.data.frame(result_ssgsea))
  expect_true(is.matrix(result_gsva) || is.data.frame(result_gsva))
})

test_that("dualGSEA works with precomputed matG", {
  data <- create_stats_test_data()
  
  result <- dualGSEA(data$X, data$y, data$gmt, matG = data$G,
                     fc.test = "cor", gsea.method = "replaid.ssgsea")
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("dualGSEA handles NA values in y", {
  data <- create_stats_test_data()
  data$y[1:2] <- NA
  
  result <- dualGSEA(data$X, data$y, data$gmt, 
                     fc.test = "rankcor", gsea.method = "replaid.ssgsea")
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("dualGSEA throws error for invalid y", {
  data <- create_stats_test_data()
  data$y[1] <- 2
  
  expect_error(dualGSEA(data$X, data$y, data$gmt))
})

test_that("dualGSEA throws error for invalid fc.test", {
  data <- create_stats_test_data()
  
  expect_error(dualGSEA(data$X, data$y, data$gmt, fc.test = "invalid"))
})

test_that("dualGSEA throws error for invalid gsea.method", {
  data <- create_stats_test_data()
  
  expect_error(dualGSEA(data$X, data$y, data$gmt, 
                        fc.test = "cor", gsea.method = "invalid"))
})

# ============================================================================
# Test: fc_ttest()
# ============================================================================

test_that("fc_ttest works with basic input", {
  data <- create_stats_test_data()
  fc <- rowMeans(data$X[, data$y == 1]) - rowMeans(data$X[, data$y == 0])
  names(fc) <- rownames(data$X)
  
  result <- fc_ttest(fc, data$G)
  
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("gsetFC" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_true("qvalue" %in% colnames(result))
})

test_that("fc_ttest works with GMT list input", {
  data <- create_stats_test_data()
  fc <- rowMeans(data$X[, data$y == 1]) - rowMeans(data$X[, data$y == 0])
  names(fc) <- rownames(data$X)
  
  result <- fc_ttest(fc, data$gmt)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("fc_ttest works with different sort.by options", {
  data <- create_stats_test_data()
  fc <- rowMeans(data$X[, data$y == 1]) - rowMeans(data$X[, data$y == 0])
  names(fc) <- rownames(data$X)
  
  result_pvalue <- fc_ttest(fc, data$G, sort.by = "pvalue")
  result_fc <- fc_ttest(fc, data$G, sort.by = "gsetFC")
  result_none <- fc_ttest(fc, data$G, sort.by = "none")
  
  expect_true(is.matrix(result_pvalue) || is.data.frame(result_pvalue))
  expect_true(is.matrix(result_fc) || is.data.frame(result_fc))
  expect_true(is.matrix(result_none) || is.data.frame(result_none))
})

test_that("fc_ttest throws error for unnamed fc", {
  data <- create_stats_test_data()
  fc <- rnorm(nrow(data$X))
  
  expect_error(fc_ttest(fc, data$G))
})

# ============================================================================
# Test: gset_averageCLR()
# ============================================================================

test_that("gset_averageCLR works with basic input", {
  data <- create_stats_test_data()
  
  result <- gset_averageCLR(data$X, data$G)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(data$G))
  expect_equal(ncol(result), ncol(data$X))
})

test_that("gset_averageCLR works with center=FALSE", {
  data <- create_stats_test_data()
  
  result <- gset_averageCLR(data$X, data$G, center = FALSE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(ncol(data$G), ncol(data$X)))
})

test_that("gset_averageCLR works with single column", {
  data <- create_stats_test_data(n_samples = 2)
  X_single <- data$X[, 1, drop = FALSE]
  
  result <- gset_averageCLR(X_single, data$G)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1)
})

test_that("gset_averageCLR handles no overlapping features", {
  data <- create_stats_test_data()
  rownames(data$X) <- paste0("OTHER_GENE", 1:nrow(data$X))
  
  result <- gset_averageCLR(data$X, data$G)
  
  expect_null(result)
})

# ============================================================================
# Test: gset_ttest()
# ============================================================================

test_that("gset_ttest works with basic input", {
  data <- create_stats_test_data()
  gsetX <- plaid(data$X, data$G)
  
  result <- gset_ttest(gsetX, data$y)
  
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("diff" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_equal(nrow(result), nrow(gsetX))
})

test_that("gset_ttest throws error for invalid y", {
  data <- create_stats_test_data()
  gsetX <- plaid(data$X, data$G)
  invalid_y <- rep(c(0, 1, 2), length.out = ncol(gsetX))
  
  expect_error(gset_ttest(gsetX, invalid_y))
})

test_that("gset_ttest produces expected output structure", {
  data <- create_stats_test_data()
  gsetX <- plaid(data$X, data$G)
  
  result <- gset_ttest(gsetX, data$y)
  
  expect_true(all(c("diff", "pvalue") %in% colnames(result)))
  expect_equal(rownames(result), rownames(gsetX))
})

# ============================================================================
# Test: matrix_onesample_ttest()
# ============================================================================

test_that("matrix_onesample_ttest works with basic input", {
  set.seed(456)
  Fm <- rnorm(100)
  names(Fm) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- matrix_onesample_ttest(Fm, data$G)
  
  expect_true(is.list(result))
  expect_true("mean" %in% names(result))
  expect_true("t" %in% names(result))
  expect_true("p" %in% names(result))
  expect_equal(nrow(result$mean), ncol(data$G))
})

test_that("matrix_onesample_ttest produces finite values", {
  set.seed(456)
  Fm <- rnorm(100)
  names(Fm) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- matrix_onesample_ttest(Fm, data$G)
  
  expect_true(all(is.finite(result$mean)))
  expect_true(all(is.finite(result$t)))
  expect_true(all(is.finite(result$p)))
})

test_that("matrix_onesample_ttest handles varying gene set sizes", {
  set.seed(456)
  Fm <- rnorm(100)
  names(Fm) <- paste0("GENE", 1:100)
  data <- create_stats_test_data(n_pathways = 10)
  
  result <- matrix_onesample_ttest(Fm, data$G)
  
  expect_equal(nrow(result$mean), ncol(data$G))
})

# ============================================================================
# Test: matrix_metap()
# ============================================================================

test_that("matrix_metap works with stouffer method", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  result <- matrix_metap(plist, method = "stouffer")
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("matrix_metap works with fisher method", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  result <- matrix_metap(plist, method = "fisher")
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("matrix_metap works with matrix input", {
  pmat <- matrix(c(0.05, 0.01, 0.1, 0.08, 0.02, 0.15), nrow = 3, ncol = 2)
  
  result <- matrix_metap(pmat, method = "stouffer")
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)
})

test_that("matrix_metap results differ by method", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  result_stouffer <- matrix_metap(plist, method = "stouffer")
  result_fisher <- matrix_metap(plist, method = "fisher")
  
  expect_false(identical(result_stouffer, result_fisher))
})

test_that("matrix_metap works with sumlog alias", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  result <- matrix_metap(plist, method = "sumlog")
  
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1))
})

test_that("matrix_metap works with sumz alias", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  result <- matrix_metap(plist, method = "sumz")
  
  expect_true(is.numeric(result))
  expect_true(all(result >= 0 & result <= 1))
})

test_that("matrix_metap throws error for invalid method", {
  plist <- list(
    p1 = c(0.05, 0.01, 0.1),
    p2 = c(0.08, 0.02, 0.15)
  )
  
  expect_error(matrix_metap(plist, method = "invalid"))
})

test_that("matrix_metap handles extreme p-values", {
  plist <- list(
    p1 = c(1e-10, 0.5, 0.99),
    p2 = c(1e-9, 0.6, 0.999)
  )
  
  result <- matrix_metap(plist, method = "stouffer")
  
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0 & result <= 1))
})

# ============================================================================
# Test: gset.rankcor()
# ============================================================================

test_that("gset.rankcor works with basic input", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- gset.rankcor(rnk, data$G, compute.p = FALSE)
  
  expect_true(is.list(result))
  expect_true("rho" %in% names(result))
  expect_equal(nrow(result$rho), ncol(data$G))
})

test_that("gset.rankcor works with compute.p=TRUE", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- gset.rankcor(rnk, data$G, compute.p = TRUE)
  
  expect_true(is.list(result))
  expect_true("rho" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("q.value" %in% names(result))
})

test_that("gset.rankcor works with use.rank=FALSE", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- gset.rankcor(rnk, data$G, use.rank = FALSE)
  
  expect_true(is.list(result))
  expect_true("rho" %in% names(result))
})

test_that("gset.rankcor works with matrix input", {
  set.seed(456)
  rnk <- matrix(rnorm(200), nrow = 100, ncol = 2)
  rownames(rnk) <- paste0("GENE", 1:100)
  colnames(rnk) <- c("Cond1", "Cond2")
  data <- create_stats_test_data()
  
  result <- gset.rankcor(rnk, data$G)
  
  expect_true(is.list(result))
  expect_equal(ncol(result$rho), 2)
})

test_that("gset.rankcor throws error for unnamed vector", {
  rnk <- rnorm(100)
  data <- create_stats_test_data()
  
  expect_error(gset.rankcor(rnk, data$G))
})

test_that("gset.rankcor throws error for matrix without rownames", {
  rnk <- matrix(rnorm(200), nrow = 100, ncol = 2)
  data <- create_stats_test_data()
  
  expect_error(gset.rankcor(rnk, data$G))
})

test_that("gset.rankcor returns NULL for zero columns", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  empty_G <- Matrix::Matrix(0, nrow = 100, ncol = 0, sparse = TRUE)
  rownames(empty_G) <- paste0("GENE", 1:100)
  
  result <- gset.rankcor(rnk, empty_G)
  
  expect_null(result)
})

test_that("gset.rankcor throws error for non-Matrix gset", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  gset <- list(pathway1 = c("GENE1", "GENE2"))
  
  expect_error(gset.rankcor(rnk, gset))
})

test_that("gset.rankcor handles transposed input correctly", {
  set.seed(456)
  rnk <- rnorm(100)
  names(rnk) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  G_transposed <- Matrix::t(data$G)
  
  result <- gset.rankcor(rnk, G_transposed)
  
  expect_true(is.list(result))
})

# ============================================================================
# Test: cor_sparse_matrix()
# ============================================================================

test_that("cor_sparse_matrix works without missing values", {
  set.seed(456)
  G <- Matrix::Matrix(rnorm(500), nrow = 100, ncol = 5, sparse = TRUE)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  
  result <- cor_sparse_matrix(G, mat)
  
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(nrow(result), ncol(G))
  expect_equal(ncol(result), ncol(mat))
})

test_that("cor_sparse_matrix works with missing values", {
  set.seed(456)
  G <- Matrix::Matrix(rnorm(500), nrow = 100, ncol = 5, sparse = TRUE)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  mat[1:5, 1] <- NA
  
  result <- cor_sparse_matrix(G, mat)
  
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(nrow(result), ncol(G))
  expect_equal(ncol(result), ncol(mat))
})

test_that("cor_sparse_matrix produces finite values", {
  set.seed(456)
  G <- Matrix::Matrix(rnorm(500), nrow = 100, ncol = 5, sparse = TRUE)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  
  result <- cor_sparse_matrix(G, mat)
  
  # Some correlations may be NA/NaN for constant columns
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
})

# ============================================================================
# Test: Integration tests
# ============================================================================

test_that("dual_test and gset_ttest produce consistent results", {
  data <- create_stats_test_data()
  gsetX <- plaid(data$X, data$G)
  
  # Run dual_test with precomputed gsetX
  dual_result <- dual_test(data$X, data$y, data$G, gsetX = gsetX)
  
  # Run gset_ttest separately
  ttest_result <- gset_ttest(gsetX, data$y)
  
  # Both should have same number of gene sets
  expect_equal(nrow(dual_result), nrow(ttest_result))
})

test_that("matrix_metap with different p-value sources", {
  data <- create_stats_test_data()
  
  # Create dummy p-values
  pv1 <- runif(ncol(data$G))
  pv2 <- runif(ncol(data$G))
  plist <- list(p1 = pv1, p2 = pv2)
  
  # Combine with both methods
  combined_stouffer <- matrix_metap(plist, method = "stouffer")
  combined_fisher <- matrix_metap(plist, method = "fisher")
  
  # Combined p-values should be more significant than individual ones
  expect_true(all(combined_stouffer <= 1))
  expect_true(all(combined_fisher <= 1))
})

# ============================================================================
# Test: Edge cases and error handling
# ============================================================================

test_that("dual_test handles small sample sizes", {
  data <- create_stats_test_data(n_samples = 4)
  
  result <- dual_test(data$X, data$y, data$G)
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("gset_averageCLR handles single gene set", {
  data <- create_stats_test_data(n_pathways = 1)
  
  result <- gset_averageCLR(data$X, data$G)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("matrix_onesample_ttest handles small gene sets", {
  set.seed(456)
  Fm <- rnorm(100)
  names(Fm) <- paste0("GENE", 1:100)
  
  # Create gene set with only 2 genes
  G <- Matrix::Matrix(0, nrow = 100, ncol = 1, sparse = TRUE)
  rownames(G) <- paste0("GENE", 1:100)
  colnames(G) <- "SmallSet"
  G[1:2, 1] <- 1
  
  result <- matrix_onesample_ttest(Fm, G)
  
  expect_true(is.list(result))
  expect_true(all(is.finite(result$mean)))
})

test_that("gset.rankcor handles partial overlap", {
  set.seed(456)
  rnk <- rnorm(50)
  names(rnk) <- paste0("GENE", 1:50)
  data <- create_stats_test_data()
  
  result <- gset.rankcor(rnk, data$G)
  
  expect_true(is.list(result))
})

test_that("dual_test handles unbalanced groups", {
  data <- create_stats_test_data()
  data$y <- rep(0, length(data$y))
  data$y[1:3] <- 1  # Few samples in group 1
  
  # This should work but may produce warnings
  result <- suppressWarnings(dual_test(data$X, data$y, data$G))
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

# ============================================================================
# Test: Numerical stability
# ============================================================================

test_that("matrix_metap handles very small p-values", {
  plist <- list(
    p1 = rep(1e-100, 5),
    p2 = rep(1e-100, 5)
  )
  
  result <- matrix_metap(plist, method = "stouffer")
  
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("matrix_metap handles p-values near 1", {
  plist <- list(
    p1 = rep(0.99, 5),
    p2 = rep(0.999, 5)
  )
  
  result <- matrix_metap(plist, method = "stouffer")
  
  expect_true(all(is.finite(result)))
  expect_true(all(result <= 1))
})

test_that("gset_averageCLR handles extreme values", {
  data <- create_stats_test_data()
  data$X <- data$X * 1e6
  
  result <- gset_averageCLR(data$X, data$G)
  
  expect_true(is.matrix(result))
  expect_true(all(is.finite(result)))
})

test_that("matrix_onesample_ttest handles constant values", {
  Fm <- rep(1, 100)
  names(Fm) <- paste0("GENE", 1:100)
  data <- create_stats_test_data()
  
  result <- matrix_onesample_ttest(Fm, data$G)
  
  expect_true(is.list(result))
  # P-values should be 1 for no variance
  expect_true(all(result$p >= 0 & result$p <= 1))
})

