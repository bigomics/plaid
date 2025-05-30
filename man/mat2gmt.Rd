% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmt-utils.R
\name{mat2gmt}
\alias{mat2gmt}
\title{Convert binary matrix back to GMT list}
\usage{
mat2gmt(mat)
}
\arguments{
\item{mat}{A matrix with non-zero entries representing genes in
each gene set where rows represent genes and columns represent
gene sets.}
}
\value{
A list of vector representing each gene set. Each list
element correspond to a gene set and is a vector of genes
}
\description{
This function converts binary matrix to a GMT (Gene Matrix
Transposed) list, The binary matrix indicates the presence or
absence of genes in each gene set, where rows represent genes and
columns represent gene sets.
}
