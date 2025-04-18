% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdriver_functions.R
\name{parallel_function}
\alias{parallel_function}
\title{Parallel Execution Function for RVdriver}
\usage{
parallel_function(
  iter,
  gene_df,
  syn_df,
  synon_threshold,
  g,
  seed_list,
  weighted = FALSE,
  non_synon_muts_list,
  sampling_method = sampling_method,
  include_gene_synons,
  rounding_method
)
}
\arguments{
\item{iter}{Integer index of iteration.}

\item{gene_df}{Data-frame with the gene's non-synonymous mutations.}

\item{syn_df}{Data-frame with the mutated samples' synonymous mutations (sampled).}

\item{synon_threshold}{Number of synonymous mutations to sample (dependent on cohort size)}

\item{g}{Gene name.}

\item{seed_list}{Vector of integer seeds for random sampling.}

\item{weighted}{Logical; passed to \code{get_results_df_lm}.}

\item{non_synon_muts_list}{Named list of data frames of non-synonymous mutations per patient.}

\item{sampling_method}{A character specifying the sampling approach (\code{"strict"} or \code{"relaxed"}).}

\item{include_gene_synons}{Logical; if \code{TRUE}, try sampling synonymous from the same gene.}

\item{rounding_method}{Method for determining the number of synonymous mutations to be samples per patient (default: \code{"expo"}).}
}
\value{
A list with elements \code{results} (data frame) and
\code{proportion_sampled} (data frame).
}
\description{
Internal function used by \code{rvdriver} to handle a single iteration of parallel
sampling and linear modelling.
}
