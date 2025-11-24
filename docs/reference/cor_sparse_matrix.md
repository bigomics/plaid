# Calculate sparse correlation matrix handling missing values

Calculate sparse correlation matrix handling missing values

## Usage

``` r
cor_sparse_matrix(G, mat)
```

## Arguments

- G:

  Sparse matrix containing gene sets

- mat:

  Matrix of values

## Value

Correlation matrix between G and mat

## Details

If mat has no missing values, calculates correlation directly using
corSparse. Otherwise computes column-wise correlations only using
non-missing values.
