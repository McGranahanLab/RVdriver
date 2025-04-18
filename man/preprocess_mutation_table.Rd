% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvdriver_functions.R
\name{preprocess_mutation_table}
\alias{preprocess_mutation_table}
\title{Preprocess Mutation Table}
\usage{
preprocess_mutation_table(
  mutation_table,
  non_synonymous_key,
  depth_filt,
  num_muts_filt,
  max_coding_muts_per_sample,
  max_muts_per_gene_per_sample
)
}
\arguments{
\item{mutation_table}{A data frame or tibble containing the mutation data.
Expected to have at least the columns \code{patient_id}, \code{RNA_alt_count}, \code{func},
\code{gene}, \code{RNA_depth}, \code{pos}, and \code{chrom}.}

\item{non_synonymous_key}{A character vector of mutation types considered as non-synonymous.}

\item{depth_filt}{Numeric. The minimum RNA depth at a mutated position for the mutation to be considered}

\item{num_muts_filt}{Numeric. Minimum number of mutations (per gene) for the gene to be considered by RVdriver}

\item{max_coding_muts_per_sample}{Numeric. Maximum number of mutations allowed per sample.}

\item{max_muts_per_gene_per_sample}{Numeric. Maximum number of mutations allowed per gene for each sample.}
}
\value{
Processed mutation table:
\describe{
\item{\code{mutation_table}}{The filtered mutation table.}
}
}
\description{
Performs preprocessing on a mutation table. This includes:
\itemize{
\item Removing samples with no expressed mutations.
\item Excluding samples with more than a specified number of mutations per sample.
\item Removing dinucleotide/trinucleotide (adjacent) mutations.
\item Excluding mutations in a gene when a sample has more than a maximum allowed number.
}

For each step the function records how many mutations and/or samples were removed.
}
