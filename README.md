[![codecov](https://codecov.io/github/bigomics/plaid/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/plaid)

# PLAID: ultrafast single-sample gene set enrichment scoring

Plaid (Pathway Level Average Intensity Detection) is an ultra-fast
method to compute single-sample enrichment scores for gene expression
or proteomics data. For each sample, plaid computes the gene set score
as the average intensity of the genes/proteins in the gene set. The
output is a gene set score matrix suitable for further analyses.

A distinctive feature of PLAID is that it can simulate few of the most
widely used gene set scoring algorithms (GSVA, ssGSEA, scSE, ucell, sing),
enabling researchers to replace those functions and gaining much improved
runtime efficiency and memory requirement. Typically, PLAID can be more than 
100 times faster than the original algorithm.

Plaid is freely available on GitHub. It's a main gene sets scoring
algorithm in OmicsPlayground, our Bioinformatics platform at BigOmics
Analytics. In OmicsPlayground, you can perform Plaid without coding
needs.

## Installation

You can install plaid from Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("plaid")
```

You can also install the development version from GitHub:

```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("bigomics/plaid")
```

## Usage

For detailed usage examples and tutorials, please see our vignettes:

- [Getting Started with PLAID](vignettes/plaid-vignette.Rmd)
- [Comparing Methods](vignettes/compare-vignette.Rmd)

**Key features:**
- Ultra-fast single-sample gene set enrichment scoring
- Automatically detects and handles Bioconductor objects (`SummarizedExperiment`, `SingleCellExperiment`, `BiocSet`)
- Works with regular matrices, sparse matrices, and Bioconductor data structures
- Includes multiple scoring methods (plaid, sing, ssgsea, scSE, ucell, gsva)
- Built-in differential enrichment testing

## References

- Zito A., et al. PLAID: ultrafast single-sample gene set enrichment scoring. [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.06.14.659661v1)

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch
