# Chunked computation of cross product

Compute crossprod (t(x) %\*% y) for very large y by computing in chunks.

## Usage

``` r
chunked_crossprod(x, y, chunk = NULL)
```

## Arguments

- x:

  Matrix First matrix for multiplication. Can be sparse.

- y:

  Matrix Second matrix for multiplication. Can be sparse.

- chunk:

  Integer Chunk size (max number of columns) for computation.

## Value

Matrix. Result of matrix cross product.
