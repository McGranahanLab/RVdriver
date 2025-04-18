% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdriver_functions.R
\name{get_rna_max_weights}
\alias{get_rna_max_weights}
\title{Compute RNA Weights based on the max coverage across all variants}
\usage{
get_rna_max_weights(coverage)
}
\arguments{
\item{coverage}{A numeric vector of RNA depths across variants}
}
\value{
A numeric vector of weights with a maximum value of 1 (for the variant
that had the highest coverage).
}
\description{
\code{get_rna_max_weights} computes a per-variant weight by dividing coverage
by the maximum coverage found among variants.
}
