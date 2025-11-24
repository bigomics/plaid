# Matrix version for combining p-values using fisher or stouffer method. Much faster than doing metap::sumlog() and metap::sumz()

Matrix version for combining p-values using fisher or stouffer method.
Much faster than doing metap::sumlog() and metap::sumz()

## Usage

``` r
matrix_metap(plist, method = "stouffer")
```

## Arguments

- plist:

  List of p-value vectors or matrix of p-values

- method:

  Method for combining p-values ("fisher"/"sumlog" or "stouffer"/"sumz")

## Value

Vector of combined p-values.
