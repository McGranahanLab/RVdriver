% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdriver_functions.R
\name{get_lm_summary}
\alias{get_lm_summary}
\title{Summarising functon for an lm}
\usage{
get_lm_summary(lm, g)
}
\arguments{
\item{lm}{An object of class \code{\link{lm}}.}

\item{g}{A character string, the gene name.}
}
\value{
A data frame with columns \code{method, gene, df, effect_size, pval}.
\code{pval} is halved if \code{tval > 0} or set to \code{1 - pval/2} otherwise.
}
\description{
\code{get_lm_summary} extracts the t-value and p-value for the mutation_func
from an \code{\link{lm}} model fit, adjusts the p-value
for a one-sided test, and returns a data frame.
}
