% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdriver_functions.R
\name{get_rna_vaf_res_lm}
\alias{get_rna_vaf_res_lm}
\title{Fit a Linear Model for RNA VAF}
\usage{
get_rna_vaf_res_lm(df, g, weighted)
}
\arguments{
\item{df}{A data frame with columns \code{RNA_VAF, func_new, rna_vaf_max_weight}.}

\item{g}{A character string, the gene name.}

\item{weighted}{Logical; if \code{TRUE}, use weights = \code{rna_vaf_max_weight}.}
}
\value{
A data frame from \code{get_lm_summary}.
}
\description{
\code{get_rna_vaf_res_lm} fits a linear model (optionally weighted) of \code{RNA_VAF}
~ \code{func_new}, and extracts summary information via \code{get_lm_summary}.
}
