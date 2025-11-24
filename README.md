# PLAID: ultrafast single-sample gene set enrichment  <img src='man/figures/logo.png' align="right" width="150"/>

[![codecov](https://codecov.io/github/bigomics/plaid/graph/badge.svg?token=66J6W41C0G)](https://codecov.io/github/bigomics/plaid)

[PLAID](https://bigomics.github.io/plaid) (Pathway Level Average Intensity Detection) is an ultrafast
method to compute single-sample enrichment scores for gene expression
or proteomics data. For each sample, PLAID computes the gene set score
as the average intensity of the genes/proteins in the gene set. The
output is a gene set score matrix suitable for further analyses.

A distinctive feature of PLAID is that it can simulate few of the most
widely used single-sample gene set scoring algorithms 
([ssGSEA](https://github.com/rcastelo/GSVA),
[GSVA](https://github.com/rcastelo/GSVA),
[AUCell](https://github.com/aertslab/AUCell),
[singscore](https://github.com/DavisLaboratory/singscore),
[scSE](https://doi.org/10.1093/nar/gkz601),
[UCell](https://github.com/carmonalab/UCell)), enabling researchers
to replace those functions and gain much improved runtime efficiency
and memory requirement. Typically, PLAID can be more than 100 times
faster and requiring 10 times less memory than the original algorithm.

#### Key features

- Ultra-fast single-sample gene set enrichment scoring
- Includes multiple scoring methods (plaid, singscore, ssGSEA, GSVA, scSE, UCell, AUCell)
- Works with regular matrices, sparse matrices, and Bioconductor data structures
- Automatically detects and handles Bioconductor objects 
([`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html), 
[`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html), [`BiocSet`](https://bioconductor.org/packages/release/bioc/html/BiocSet.html))
- Built-in differential enrichment testing

#### Warning

PLAID is fast. Ludicrously fast. Please fasten your seatbelts before usage.

## Installation

You can install PLAID from Bioconductor:

```r
BiocManager::install("plaid")
```

You can also install the development version from GitHub:

```r
remotes::install_github("bigomics/plaid")
```

## Usage

For detailed usage examples and tutorials, please see our vignettes:

- [Getting Started with PLAID](https://bigomics.github.io/plaid/articles/01_plaid-vignette.html)
- [Comparing PLAID with other methods](https://bigomics.github.io/plaid/articles/02_compare-vignette.html)

PLAID is the main single-sample gene set scoring algorithm in OmicsPlayground, our 
Bioinformatics platform at [BigOmics](https://bigomics.ch). In OmicsPlayground, you 
can perform PLAID without coding needs.

## References

For more technical details please refer to our papers. Please cite us when you use
PLAID as part of your research. 

- Zito A., et al. PLAID: ultrafast single-sample gene set enrichment scoring. Bioinformatics, 2025, [btaf621](https://doi.org/10.1093/bioinformatics/btaf621).
- Akhmedov M., et al., Omics Playground: a comprehensive self-service platform for visualization, analytics and exploration of Big Omics Data, NAR Genomics and Bioinformatics, 2020, [lqz019](https://doi.org/10.1093/nargab/lqz019).

## Support

For support feel free to reach our Bioinformatics Data Science Team at
BigOmics Analytics: help@bigomics.ch

If you like PLAID, please recommend us to your friends, buy us [coffee](https://buymeacoffee.com/bigomics)
and brag about PLAID on your social media. 
