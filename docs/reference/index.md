# Package index

## All functions

- [`chunked_crossprod()`](https://bigomics.github.io/plaid/reference/chunked_crossprod.md)
  : Chunked computation of cross product
- [`colranks()`](https://bigomics.github.io/plaid/reference/colranks.md)
  : Compute columnwise ranks of matrix
- [`cor_sparse_matrix()`](https://bigomics.github.io/plaid/reference/cor_sparse_matrix.md)
  : Calculate sparse correlation matrix handling missing values
- [`dualGSEA()`](https://bigomics.github.io/plaid/reference/dualGSEA.md)
  : Reimplementation of dualGSEA (Bull et al., 2024) but defaults with
  replaid backend. For the preranked test we still use fgsea. Should be
  much faster than original using fgsea + GSVA::ssGSEA.
- [`fc_ttest()`](https://bigomics.github.io/plaid/reference/fc_ttest.md)
  : T-test statistical testing of differentially enrichment
- [`fc_ztest()`](https://bigomics.github.io/plaid/reference/fc_ztest.md)
  : Z-test statistical testing of differentially enrichment
- [`gmt2mat()`](https://bigomics.github.io/plaid/reference/gmt2mat.md) :
  Convert GMT to Binary Matrix
- [`gset.rankcor()`](https://bigomics.github.io/plaid/reference/gset.rankcor.md)
  : Calculate gene set rank correlation
- [`gset_averageCLR()`](https://bigomics.github.io/plaid/reference/gset_averageCLR.md)
  : Compute geneset expression as the average log-ration of genes in the
  geneset. Requires log-expression matrix X and (sparse) geneset matrix
  matG.
- [`gset_ttest()`](https://bigomics.github.io/plaid/reference/gset_ttest.md)
  : Perform t-test on gene set scores
- [`mat.rowsds()`](https://bigomics.github.io/plaid/reference/mat.rowsds.md)
  : Calculate row standard deviations for matrix
- [`mat2gmt()`](https://bigomics.github.io/plaid/reference/mat2gmt.md) :
  Convert Binary Matrix to GMT
- [`matrix_metap()`](https://bigomics.github.io/plaid/reference/matrix_metap.md)
  : Matrix version for combining p-values using fisher or stouffer
  method. Much faster than doing metap::sumlog() and metap::sumz()
- [`matrix_onesample_ttest()`](https://bigomics.github.io/plaid/reference/matrix_onesample_ttest.md)
  : Perform one-sample t-test on matrix with gene sets
- [`normalize_medians()`](https://bigomics.github.io/plaid/reference/normalize_medians.md)
  : Normalize column medians of matrix
- [`plaid()`](https://bigomics.github.io/plaid/reference/plaid.md) :
  Compute PLAID single-sample enrichment score
- [`read.gmt()`](https://bigomics.github.io/plaid/reference/read.gmt.md)
  : Read GMT File
- [`replaid.aucell()`](https://bigomics.github.io/plaid/reference/replaid.aucell.md)
  : Fast calculation of AUCell
- [`replaid.gsva()`](https://bigomics.github.io/plaid/reference/replaid.gsva.md)
  : Fast approximation of GSVA
- [`replaid.scse()`](https://bigomics.github.io/plaid/reference/replaid.scse.md)
  : Fast calculation of scSE score
- [`replaid.sing()`](https://bigomics.github.io/plaid/reference/replaid.sing.md)
  : Fast calculation of singscore
- [`replaid.ssgsea()`](https://bigomics.github.io/plaid/reference/replaid.ssgsea.md)
  : Fast calculation of ssGSEA
- [`replaid.ucell()`](https://bigomics.github.io/plaid/reference/replaid.ucell.md)
  : Fast calculation of UCell
- [`sparse_colranks()`](https://bigomics.github.io/plaid/reference/sparse_colranks.md)
  : Compute columm ranks for sparse matrix. Internally used by
  colranks()
- [`write.gmt()`](https://bigomics.github.io/plaid/reference/write.gmt.md)
  : Write GMT File
