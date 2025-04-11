#' @title Compute RNA Weights
#' @description
#' \code{get_rna_weights} computes a per-variant weight based on the fraction of
#' reads for the variant, against all other variants
#'
#' @param coverage A numeric vector of RNA depth of coverage
#'
#' @return A numeric vector of weights
#' @export
get_rna_weights <- function(coverage){
  # coverage: RNA depth at mutated positions
  weighted_expr <- coverage / sum(coverage)
  return(weighted_expr)
}


#' @title Compute RNA Weights based on the max coverage across all variants
#' @description
#' \code{get_rna_max_weights} computes a per-variant weight by dividing coverage
#' by the maximum coverage found among variants.
#'
#' @param coverage A numeric vector of RNA depths across variants
#'
#' @return A numeric vector of weights with a maximum value of 1 (for the variant
#'   that had the highest coverage).
#' @export
get_rna_max_weights <- function(coverage){
  weighted_expr <- coverage / max(coverage)
  return(weighted_expr)
}


#' @title Summarising functon for an lm
#' @description
#' \code{get_lm_summary} extracts the t-value and p-value for the mutation_func
#' from an \code{\link{lm}} model fit, adjusts the p-value
#' for a one-sided test, and returns a data frame.
#'
#' @param lm An object of class \code{\link{lm}}.
#' @param g A character string, the gene name.
#'
#' @return A data frame with columns \code{method, gene, df, effect_size, pval}.
#'   \code{pval} is halved if \code{tval > 0} or set to \code{1 - pval/2} otherwise.
#' @export
get_lm_summary <- function(lm, g) {
  coef_table <- summary(lm)$coefficients

  # Extract the t-value and p-value from row "func_newnon_synonymous"
  tval <- coef_table["func_newnon_synonymous", "t value"]
  pval <- coef_table["func_newnon_synonymous", "Pr(>|t|)"]

  # Adjust p-value for a one-sided test
  if(tval > 0) {
    pval <- pval / 2
  } else {
    pval <- 1 - (pval / 2)
  }

  return(data.frame(
    method = "RVdriver_lm",
    gene = g,
    df = NA,
    effect_size = tval,
    pval = pval
  ))
}


#' @title Run get_rna_vaf_res_lm using a TryCatch function
#' @description
#' \code{get_results_df_lm} runs \code{get_rna_vaf_res_lm} in a try-catch block
#' (via \code{myTryCatch}), checks for warnings, and returns a result data frame.
#'
#' @param comp_df A data frame with columns \code{RNA_VAF, func_new, rna_vaf_max_weight}.
#' @param g A character string, the gene name.
#' @param weighted Logical; if \code{TRUE}, use weighted LM (weighted by mutation depth).
#'
#' @return A data frame with columns \code{method, gene, df, effect_size, pval}.
#' @export
get_results_df_lm <- function(comp_df, g, weighted = FALSE){
  rna_vaf_res <- myTryCatch(get_rna_vaf_res_lm(comp_df, g, weighted))
  if (is.null(rna_vaf_res$value)) {
    rna_vaf_res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA)
  } else if (!is.null(rna_vaf_res$warning)) {
    if(rna_vaf_res$warning$message == "essentially perfect fit: summary may be unreliable") {
      rna_vaf_res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA)
    } else {
      rna_vaf_res <- rna_vaf_res$value
    }
  } else {
    rna_vaf_res <- rna_vaf_res$value
  }
  res_df <- do.call(rbind, list(rna_vaf_res))
  return(res_df)
}


#' @title Fit a Linear Model for RNA VAF
#' @description
#' \code{get_rna_vaf_res_lm} fits a linear model (optionally weighted) of \code{RNA_VAF}
#' ~ \code{func_new}, and extracts summary information via \code{get_lm_summary}.
#'
#' @param df A data frame with columns \code{RNA_VAF, func_new, rna_vaf_max_weight}.
#' @param g A character string, the gene name.
#' @param weighted Logical; if \code{TRUE}, use weights = \code{rna_vaf_max_weight}.
#'
#' @return A data frame from \code{get_lm_summary}.
#' @export
get_rna_vaf_res_lm <- function(df, g, weighted){
  unique_cluster <- (length(unique(df$cluster)) == 1)

  if (weighted) {
    if (unique_cluster) {
      model <- lm(RNA_VAF ~ func_new, weights = rna_vaf_max_weight, data = df)
    } else {
      model <- lm(RNA_VAF ~ func_new, weights = rna_vaf_max_weight, data = df)
    }
  } else {
    if (unique_cluster) {
      model <- lm(RNA_VAF ~ func_new, data = df)
    } else {
      model <- lm(RNA_VAF ~ func_new, data = df)
    }
  }

  return(get_lm_summary(model, g))
}


#' @title Try-Catch Helper
#' @description
#' \code{myTryCatch} is a helper function that captures warnings, errors, and
#' returns them in a list.
#'
#' @param expr Expression being evaluated
#'
#' @return A list with elements \code{value}, \code{warning}, and \code{error}.
#' @export
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }),
    warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warning = warn, error = err)
}


#' @title Parallel Execution Function for RVdriver
#' @description
#' Internal function used by \code{rvdriver} to handle a single iteration of parallel
#' sampling and linear modelling.
#'
#' @param iter Integer index of iteration.
#' @param gene_df Data-frame with the gene's non-synonymous mutations.
#' @param syn_df Data-frame with the mutated samples' synonymous mutations (sampled).
#' @param synon_threshold Number of synonymous mutations to sample (dependent on cohort size)
#' @param g Gene name.
#' @param seed_list Vector of integer seeds for random sampling.
#' @param weighted Logical; passed to \code{get_results_df_lm}.
#' @param non_synon_muts_list Named list of data frames of non-synonymous mutations per patient.
#' @param sampling_method A character specifying the sampling approach (\code{"strict"} or \code{"relaxed"}).
#' @param include_gene_synons Logical; if \code{TRUE}, try sampling synonymous from the same gene.
#' @param rounding_method Method for determining the number of synonymous mutations to be samples per patient (default: \code{"expo"}).
#'
#' @return A list with elements \code{results} (data frame) and
#'   \code{proportion_sampled} (data frame).
#' @export
parallel_function <- function(iter,
                              gene_df,
                              syn_df,
                              synon_threshold,
                              g,
                              seed_list,
                              weighted = FALSE,
                              non_synon_muts_list,
                              sampling_method = sampling_method,
                              include_gene_synons,
                              rounding_method) {

  set.seed(seed_list[iter])
  res_all_cats <- data.frame()
  proportion_sampled_per_patient <- data.frame()

  func_above_1 <- gene_df %>%
    dplyr::group_by(func) %>%
    dplyr::filter(RNA_depth >= depth_filt) %>%
    dplyr::summarise(num = length(unique(patient_id))) %>%
    dplyr::filter(num > num_muts_filt)

  for (mut_func in unique(func_above_1$func)){
    mut_func_df <- gene_df %>%
      dplyr::filter(func == mut_func)

    num_mutations <- nrow(mut_func_df)
    num_samples <- length(unique(mut_func_df$patient_id))

    # dynamically computed threshold
    if(rounding_method == "ceiling") {
      synon_threshold_filt <- sqrt_ceiling(num_samples)
    } else if(rounding_method == "floor") {
      synon_threshold_filt <- sqrt_floor(num_samples)
    } else if (rounding_method == "round") {
      synon_threshold_filt <- sqrt_round(num_samples)
    } else if (rounding_method == "expo") {
      synon_threshold_filt <- exponential_growth_limited_2x(num_samples, rate = rate)
    }

    num_muts_greater_than_7_depth <- nrow(mut_func_df %>% dplyr::filter(RNA_depth > 7))
    proportion_greater_than_7_depth <- num_muts_greater_than_7_depth / num_mutations

    mut_func_df$func_new <- "non_synonymous"

    syn_df_filt <- syn_df %>%
      dplyr::filter(patient_id %in% mut_func_df$patient_id)

    # Sample synonyms
    syn_df_sampled_list <- create_synonymous_background_including_mut_type(
      syn_df_filt,
      mut_func_df,
      synon_threshold_filt,
      synonymous_background,
      seed_list[iter],
      non_synon_muts_list,
      sampling_method,
      include_gene_synons
    )

    proportion_sampled_per_patient_df <- syn_df_sampled_list$proportion_sampled_per_patient
    proportion_sampled_per_patient <- rbind(proportion_sampled_per_patient,
                                            proportion_sampled_per_patient_df)

    sampling_summary_df <- syn_df_sampled_list$sampling_summary_df
    syn_df_sampled <- syn_df_sampled_list$final_synonymous_background

    if (nrow(syn_df_sampled) == 0 ) {
      res <- data.frame(method = "RVdriver_lm", gene = g, df = NA,
                        effect_size = NA, pval = NA)
    } else {
      syn_df_sampled$func <- "synonymous"
      syn_df_sampled$func_new <- syn_df_sampled$func
      syn_df_sampled_filt <- syn_df_sampled

      df_combined <- rbind(mut_func_df, syn_df_sampled_filt)
      df_combined$func_new <- factor(df_combined$func_new,
                                     levels = c("synonymous",
                                                unique(df_combined$func_new[df_combined$func_new != "synonymous"])))

      res_weighted <- get_results_df_lm(df_combined, g, weighted = TRUE) %>%
        dplyr::mutate(method = "weighted_RVdriver_lm")
      res_non_weighted <- get_results_df_lm(df_combined, g, weighted = FALSE)
      res <- rbind(res_weighted, res_non_weighted)

      res <- cbind(res, sampling_summary_df)
      res$func <- mut_func
      res$num_mutations <- num_mutations
      res$num_mutated_samples <- num_samples
      res$num_muts_greater_than_7_depth <- num_muts_greater_than_7_depth
      res$proportion_greater_than_7_depth <- proportion_greater_than_7_depth
      res$num_sampled <- synon_threshold_filt

      res_all_cats <- rbind(res_all_cats, res)
    }
  }

  # if more than one mutation category in gene_df, do a combined analysis
  if (length(unique(gene_df$func)) > 1) {
    mut_func_df <- gene_df %>%
      dplyr::mutate(func_new = "non_synonymous")

    syn_df_filt <- syn_df %>%
      dplyr::filter(patient_id %in% mut_func_df$patient_id)

    # Sample synonyms (full threshold)
    syn_df_sampled_list <- create_synonymous_background_including_mut_type(
      syn_df_filt,
      mut_func_df,
      synon_threshold,
      synonymous_background,
      seed_list[iter],
      non_synon_muts_list,
      sampling_method,
      include_gene_synons
    )

    proportion_sampled_per_patient_df <- syn_df_sampled_list$proportion_sampled_per_patient
    proportion_sampled_per_patient <- rbind(proportion_sampled_per_patient, proportion_sampled_per_patient_df)

    sampling_summary_df <- syn_df_sampled_list$sampling_summary_df
    syn_df_sampled <- syn_df_sampled_list$final_synonymous_background

    if (nrow(syn_df_sampled) == 0 ) {
      res <- data.frame(method = "RVdriver_lm", gene = g, df = NA,
                        effect_size = NA, pval = NA)
    } else {
      syn_df_sampled$func <- "synonymous"
      syn_df_sampled$func_new <- syn_df_sampled$func

      num_mutations <- nrow(mut_func_df)
      num_samples <- length(unique(mut_func_df$patient_id))
      num_muts_greater_than_7_depth <- nrow(mut_func_df %>% dplyr::filter(RNA_depth > 7))
      proportion_greater_than_7_depth <- num_muts_greater_than_7_depth / num_mutations

      df_combined <- rbind(mut_func_df, syn_df_sampled)
      df_combined$func_new <- factor(df_combined$func_new,
                                     levels = c("synonymous",
                                                unique(df_combined$func_new[df_combined$func_new != "synonymous"])))

      res_weighted <- get_results_df_lm(df_combined, g, weighted = TRUE) %>%
        dplyr::mutate(method = "weighted_RVdriver_lm")
      res_non_weighted <- get_results_df_lm(df_combined, g, weighted = FALSE)
      res <- rbind(res_weighted, res_non_weighted)

      res <- cbind(res, sampling_summary_df)
      res$func <- "all_mutation_types"
      res$num_mutations <- num_mutations
      res$num_mutated_samples <- num_samples
      res$num_muts_greater_than_7_depth <- num_muts_greater_than_7_depth
      res$proportion_greater_than_7_depth <- proportion_greater_than_7_depth
      res$num_sampled <- synon_threshold

      res_all_cats <- rbind(res_all_cats, res)
    }
  }

  return(list(results = res_all_cats, proportion_sampled = proportion_sampled_per_patient))
}


#' @title Main RVdriver Function
#' @description
#' \code{rvdriver} is the main function that samples synonymous mutations per sample,
#' computing an lm(RNA_VAF ~ mutation_function, weights = mutation_depth),
#' (via \code{parallel_function}) across multiple iterations, returning summary statistics
#' per gene.
#'
#' @param comp_df Data frame including all mutations to be utilised in each iteration.
#' @param seed_list A vector of seeds for reproducibility across iterations.
#' @param synon_threshold Number of synonymous mutations to be sampled per sample
#' @param num_iter Integer: number of iterations.
#' @param g The gene name.
#' @param ncores Number of cores for mclapply parellelisation
#' @param non_synon_muts_list Named list of data frames of nonsynonymous mutations by patient.
#' @param sampling_method \code{"strict"} or \code{"relaxed"} approach. default: relaxed
#' @param include_gene_synons Logical; if \code{TRUE}, try sampling synons from same gene. default: TRUE
#' @param rounding_method Method for determining number of synonymous mutations to be sampled (e.g. \code{"expo"}).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{res_df}: a data frame summarising the results
#'   \item \code{proportion_sampled_per_patient}: a data frame summarising how many
#'   synonyous mutations got sampled per patient
#' }
#' @export
rvdriver <- function(comp_df,
                     seed_list,
                     synon_threshold,
                     num_iter,
                     g,
                     ncores,
                     non_synon_muts_list,
                     sampling_method,
                     include_gene_synons,
                     rounding_method) {

  syn_df <- comp_df %>% dplyr::filter(func == "synonymous")
  gene_df <- comp_df %>% dplyr::filter(gene == g, func != "synonymous")

  # run parallel_function in parallel
  results_list_all <- parallel::mclapply(
    X = 1:num_iter,
    FUN = parallel_function,
    gene_df = gene_df,
    syn_df = syn_df,
    synon_threshold = synon_threshold,
    g = g,
    weighted = FALSE,
    seed_list = seed_list,
    mc.cores = ncores,
    non_synon_muts_list = non_synon_muts_list,
    sampling_method = sampling_method,
    include_gene_synons = include_gene_synons,
    rounding_method = rounding_method
  )

  # Extract results and proportion sampled
  results_list <- lapply(results_list_all, function(x) x$results)
  proportion_sampled_per_patient <- lapply(results_list_all, function(x) x$proportion_sampled)

  # Combine into data frames
  res_df <- do.call(rbind, results_list)
  proportion_sampled_per_patient <- do.call(rbind, proportion_sampled_per_patient)
  proportion_sampled_per_patient <- unique(proportion_sampled_per_patient)

  res_df_lm <- res_df %>%
    dplyr::filter(method == "RVdriver_lm")
  res_df_weighted_lm <- res_df %>%
    dplyr::filter(method == "weighted_RVdriver_lm")

  # Summaries
  weighted_df <- res_df_lm %>%
    dplyr::group_by(func) %>%
    dplyr::summarise(
      method = "non_weighted_RVdriver",
      gene = g,
      median_tval = median(effect_size, na.rm = TRUE),
      median_df = median(df, na.rm = TRUE),
      pval_sd = sd(pval, na.rm = TRUE),
      num_signif = sum(pval < 0.05, na.rm = TRUE),
      pval = exp(mean(log(pval), na.rm = TRUE)),
      number_synonymous_filtered = mean(num_sampled),
      num_muts_greater_than_7_depth = mean(num_muts_greater_than_7_depth),
      proportion_greater_than_7_depth = mean(proportion_greater_than_7_depth),
      num_mutations = mean(num_mutations),
      num_mutated_samples = mean(num_mutated_samples),
      num_samples_w_more_than_required = mean(num_samples_w_more_than_required),
      num_samples_sampling_from_mutated_samples = mean(num_samples_sampling_from_mutated_samples),
      num_samples_sampling_from_global_synon = mean(num_samples_sampling_from_global_synon),
      num_samples_sampling_from_mut_type_within_patient = mean(num_samples_sampling_from_mut_type_within_patient),
      num_samples_sampling_incomplete = mean(num_samples_sampling_incomplete),
      .groups = "drop"
    ) %>%
    dplyr::select(method, gene, func, median_tval, median_df, pval, pval_sd, num_signif,
                  num_mutations, num_mutated_samples, number_synonymous_filtered,
                  num_muts_greater_than_7_depth, proportion_greater_than_7_depth,
                  num_samples_w_more_than_required, num_samples_sampling_from_mutated_samples,
                  num_samples_sampling_from_global_synon,
                  num_samples_sampling_from_mut_type_within_patient,
                  num_samples_sampling_incomplete)

  weighted_df_lm <- res_df_weighted_lm %>%
    dplyr::group_by(func) %>%
    dplyr::summarise(
      method = "weighted_RVdriver",
      gene = g,
      median_tval = median(effect_size, na.rm = TRUE),
      median_df = median(df, na.rm = TRUE),
      pval_sd = sd(pval, na.rm = TRUE),
      num_signif = sum(pval < 0.05, na.rm = TRUE),
      pval = exp(mean(log(pval), na.rm = TRUE)),
      number_synonymous_filtered = mean(num_sampled),
      num_muts_greater_than_7_depth = mean(num_muts_greater_than_7_depth),
      proportion_greater_than_7_depth = mean(proportion_greater_than_7_depth),
      num_mutations = mean(num_mutations),
      num_mutated_samples = mean(num_mutated_samples),
      num_samples_w_more_than_required = mean(num_samples_w_more_than_required),
      num_samples_sampling_from_mutated_samples = mean(num_samples_sampling_from_mutated_samples),
      num_samples_sampling_from_global_synon = mean(num_samples_sampling_from_global_synon),
      num_samples_sampling_from_mut_type_within_patient = mean(num_samples_sampling_from_mut_type_within_patient),
      num_samples_sampling_incomplete = mean(num_samples_sampling_incomplete),
      .groups = "drop"
    ) %>%
    dplyr::select(method, gene, func, median_tval, median_df, pval, pval_sd, num_signif,
                  num_mutations, num_mutated_samples, number_synonymous_filtered,
                  num_muts_greater_than_7_depth, proportion_greater_than_7_depth,
                  num_samples_w_more_than_required, num_samples_sampling_from_mutated_samples,
                  num_samples_sampling_from_global_synon,
                  num_samples_sampling_from_mut_type_within_patient,
                  num_samples_sampling_incomplete)

  res_df <- rbind(weighted_df, weighted_df_lm)
  return(list(res_df = res_df, proportion_sampled_per_patient = proportion_sampled_per_patient))
}


#' @title Sample and Update Synonymous Mutations
#' @description
#' Internal helper function that tries to sample additional synonymous mutations from
#' different sources (patient, global, etc.) until the threshold is reached.
#'
#' @param patient_ids Vector of patient IDs.
#' @param threshold Numeric threshold of how many mutations to sample.
#' @param source_df_name Character name of the data frame object in the environment.
#' @param source_name Label for tracking in \code{sampling_summary_df}.
#' @param missing Integer number of how many mutations are still missing for the patient.
#'
#' @return Updated \code{missing} 0 how many missing mutations still remain after running the function
#' @export
sample_and_update <- function(patient_ids, threshold, source_df_name, source_name, missing) {
  for (patient_id_to_check in patient_ids) {
    while (missing > 0 && nrow(get(source_df_name)) > 0) {
      sampled_priority <- get(source_df_name) %>%
        dplyr::sample_n(min(missing, nrow(get(source_df_name))), replace = FALSE) %>%
        dplyr::select(colnames(syn_df)) %>%
        dplyr::mutate(pseudo_sample = patient_id_to_check)

      final_synonymous_background <<- rbind(final_synonymous_background, sampled_priority)
      updated_df <- get(source_df_name) %>%
        dplyr::anti_join(sampled_priority %>% dplyr::select(-pseudo_sample),
                         by = colnames(get(source_df_name)))
      assign(source_df_name, updated_df, envir = .GlobalEnv)

      missing <- threshold - nrow(final_synonymous_background %>%
                                    dplyr::filter(pseudo_sample == patient_id_to_check))
      sampling_summary_df[[source_name]] <<- sampling_summary_df[[source_name]] + 1
    }
  }
  return(missing)
}


#' @title Create Synonymous Background (with mutation-type sampling)
#' @description
#' \code{create_synonymous_background_including_mut_type} attempts to sample enough
#' synonymous mutations per tumour - in a heirarchical order:
#' Consider synonymous mutations only within the tumour harbouring a mutation within the gene of interest
#' Consider synonymous mutations within all tumours harbouring a mutation within the gene of interest
#' Consider all synonymous mutations from tumours within the cancer type being considered
#' Consider mutations of the same type from within that sample.
#' If sufficient mutations have not been sampled at this point, proceed with the incompletely sampled background.
#'
#' @param syn_df Data frame of synonymous mutations for relevant patients.
#' @param gene_df Data frame of non-synonymous mutations for the gene of interest
#' @param synon_threshold Integer specifying how many synonymous mutations to try to sample.
#' @param synonymous_background A data frame of all synonymous mutations across the cohort
#' @param seed An integer seed for reproducibility.
#' @param non_synon_muts_list Named list of data frames of non-synonymous mutations by patient.
#' @param sampling_method "strict" or "relaxed" - Default: relaxed
#' @param include_gene_synons Logical; if \code{TRUE}, Synonymous mutations from the gene of interest are also sampled
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{final_synonymous_background}: the sampled synonymous mutations
#'   \item \code{sampling_summary_df}: a data frame with summary of from where synonymous mutations were sampled
#'   \item \code{proportion_sampled_per_patient}: fraction of required synonymous mutations that actually got sampled for each patient
#' }
#' @export
create_synonymous_background_including_mut_type <- function(
    syn_df,
    gene_df,
    synon_threshold,
    synonymous_background_filtered,
    seed,
    non_synon_muts_list,
    sampling_method,
    include_gene_synons = FALSE
) {
  set.seed(seed)

  syn_df_summary <- syn_df %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(synon_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(patient_id = gene_df$patient_id, fill = list(synon_count = 0))

  patients_above <- syn_df_summary %>%
    dplyr::filter(synon_count >= synon_threshold) %>%
    dplyr::pull(patient_id)
  patients_below <- syn_df_summary %>%
    dplyr::filter(synon_count > 0 & synon_count < synon_threshold) %>%
    dplyr::pull(patient_id)
  patients_with_0 <- syn_df_summary %>%
    dplyr::filter(synon_count == 0) %>%
    dplyr::pull(patient_id)

  final_synonymous_background <- data.frame(matrix(ncol = ncol(syn_df) + 1, nrow = 0))
  colnames(final_synonymous_background) <- c(colnames(syn_df), "pseudo_sample")

  # Force matching classes
  for (col in colnames(syn_df)) {
    class(final_synonymous_background[[col]]) <- class(syn_df[[col]])
  }
  final_synonymous_background$pseudo_sample <- character()

  sampling_summary_df <- data.frame(
    num_samples_w_more_than_required = length(patients_above),
    num_samples_sampling_from_mutated_samples = 0,
    num_samples_sampling_from_global_synon = 0,
    num_samples_sampling_from_mut_type_within_patient = 0,
    num_samples_sampling_incomplete = 0
  )

  synonymous_background_filtered <- synonymous_background_filtered %>%
    dplyr::select(colnames(syn_df)) %>%
    dplyr::mutate(func = "synonymous") %>%
    dplyr::filter(!patient_id %in% unique(gene_df$patient_id))

  if(include_gene_synons){
    # first sample any synonymous mutations from our gene of interest
    gene_synons <- syn_df %>%
      dplyr::filter(gene == unique(gene_df$gene)) %>%
      dplyr::mutate(pseudo_sample = patient_id)
    gene_synons_all <- synonymous_background_filtered %>%
      dplyr::filter(gene == unique(gene_df$gene)) %>%
      dplyr::mutate(pseudo_sample = patient_id)
    all_gene_synons <- rbind(gene_synons, gene_synons_all)

    remaining_synonymous <- syn_df %>%
      dplyr::anti_join(all_gene_synons, by = colnames(syn_df))
    synonymous_background_filtered <- synonymous_background_filtered %>%
      dplyr::anti_join(all_gene_synons, by = colnames(syn_df))

    final_synonymous_background <- rbind(final_synonymous_background, all_gene_synons)
  }

  if(!include_gene_synons){
    remaining_synonymous <- syn_df # default if gene synonymous mutations not forced
  }

  if (length(patients_above) > 0) {
    # sample enough synonymous mutations from patients that have more than required
    for (patient_id_to_check in patients_above) {
      patient_syn_df <- syn_df %>%
        dplyr::filter(patient_id == patient_id_to_check)
      sampled_synonymous <- patient_syn_df %>%
        dplyr::sample_n(synon_threshold, replace = FALSE) %>%
        dplyr::mutate(pseudo_sample = patient_id_to_check)
      final_synonymous_background <- rbind(final_synonymous_background, sampled_synonymous)
    }
    # recalc remaining
    remaining_synonymous <- syn_df %>%
      dplyr::anti_join(final_synonymous_background, by = colnames(syn_df))
  }

  # handle patients_below
  for (patient_id_to_check in patients_below) {
    patient_syn_df <- syn_df %>% dplyr::filter(patient_id == patient_id_to_check)
    sampled_synonymous <- patient_syn_df %>% dplyr::mutate(pseudo_sample = patient_id_to_check)
    final_synonymous_background <- rbind(final_synonymous_background, sampled_synonymous)
  }
  remaining_synonymous <- syn_df %>%
    dplyr::anti_join(final_synonymous_background, by = colnames(syn_df))

  # define local sample_and_update
  sample_and_update <- function(patient_ids, threshold, source_df, source_name, missing) {
    for (patient_id_to_check in patient_ids) {
      while (missing > 0 && nrow(source_df) > 0) {
        source_df_before <- nrow(source_df)
        sampled_priority <- source_df %>%
          dplyr::sample_n(min(missing, nrow(source_df)), replace = FALSE) %>%
          dplyr::select(colnames(syn_df)) %>%
          dplyr::mutate(pseudo_sample = patient_id_to_check)

        final_synonymous_background <<- rbind(final_synonymous_background, sampled_priority)
        source_df <- source_df %>%
          dplyr::anti_join(sampled_priority %>% dplyr::select(-pseudo_sample),
                           by = colnames(source_df))

        source_df_after <- nrow(source_df)
        if (source_df_after >= source_df_before) {
          stop(paste("Error: Rows not removed from", source_name, "for patient", patient_id_to_check))
        }

        missing <- threshold - nrow(final_synonymous_background %>%
                                      dplyr::filter(pseudo_sample == patient_id_to_check))
      }
    }
    return(list(missing = missing, source_df = source_df))
  }

  # now fill in from smaller pools
  for (patient_id_to_check in patients_below) {
    num_sampled <- nrow(final_synonymous_background %>%
                          dplyr::filter(pseudo_sample == patient_id_to_check))
    missing <- synon_threshold - num_sampled

    if (sampling_method == "strict") {
      missing <- 0
      next
    }

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mutated_samples <-
        sampling_summary_df$num_samples_sampling_from_mutated_samples + 1
      next
    }

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                remaining_synonymous, "remaining_synonymous", missing)
    missing <- result$missing
    remaining_synonymous <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mutated_samples <-
        sampling_summary_df$num_samples_sampling_from_mutated_samples + 1
      next
    }

    if (sampling_method == "strict") {
      missing <- 0
      next
    }

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                synonymous_background_filtered, "synonymous_background_filtered", missing)
    missing <- result$missing
    synonymous_background_filtered <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_global_synon <-
        sampling_summary_df$num_samples_sampling_from_global_synon + 1
      next
    }

    mut_type_avail <- non_synon_muts_list[[patient_id_to_check]]$patient_df %>%
      dplyr::filter(gene != unique(gene_df$gene),
                    func %in% unique(gene_df$func[gene_df$patient_id == patient_id_to_check])) %>%
      dplyr::filter(RNA_depth >= depth_filt) %>%
      dplyr::select(colnames(syn_df))

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                mut_type_avail, "mut_type_avail", missing)
    missing <- result$missing
    mut_type_avail <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mut_type_within_patient <-
        sampling_summary_df$num_samples_sampling_from_mut_type_within_patient + 1
      next
    }
    if (missing > 0) {
      sampling_summary_df$num_samples_sampling_incomplete <-
        sampling_summary_df$num_samples_sampling_incomplete + 1
    }
  }

  # handle patients_with_0 similarly
  for (patient_id_to_check in patients_with_0) {
    num_sampled <- nrow(final_synonymous_background %>%
                          dplyr::filter(pseudo_sample == patient_id_to_check))
    missing <- synon_threshold - num_sampled

    if (sampling_method == "strict") {
      missing <- 0
      next
    }

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mutated_samples <-
        sampling_summary_df$num_samples_sampling_from_mutated_samples + 1
      next
    }

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                remaining_synonymous, "remaining_synonymous", missing)
    missing <- result$missing
    remaining_synonymous <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mutated_samples <-
        sampling_summary_df$num_samples_sampling_from_mutated_samples + 1
      next
    }

    if (sampling_method == "strict") {
      missing <- 0
      next
    }

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                synonymous_background_filtered, "synonymous_background_filtered", missing)
    missing <- result$missing
    synonymous_background_filtered <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_global_synon <-
        sampling_summary_df$num_samples_sampling_from_global_synon + 1
      next
    }

    mut_type_avail <- non_synon_muts_list[[patient_id_to_check]]$patient_df %>%
      dplyr::filter(gene != unique(gene_df$gene),
                    func %in% unique(gene_df$func[gene_df$patient_id == patient_id_to_check])) %>%
      dplyr::filter(RNA_depth >= depth_filt) %>%
      dplyr::select(colnames(syn_df))

    result <- sample_and_update(c(patient_id_to_check), synon_threshold,
                                mut_type_avail, "mut_type_avail", missing)
    missing <- result$missing
    mut_type_avail <- result$source_df

    if(missing == 0) {
      sampling_summary_df$num_samples_sampling_from_mut_type_within_patient <-
        sampling_summary_df$num_samples_sampling_from_mut_type_within_patient + 1
      next
    }
    if (missing > 0) {
      sampling_summary_df$num_samples_sampling_incomplete <-
        sampling_summary_df$num_samples_sampling_incomplete + 1
    }
  }

  # Proportion of synonymous mutations that got sampled
  proportion_sampled_per_patient <- syn_df %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      total_synon_mutations = dplyr::n(),
      sampled_synon_mutations = sum(
        paste(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func) %in%
          paste(final_synonymous_background$patient_id,
                final_synonymous_background$RNA_VAF,
                final_synonymous_background$gene,
                final_synonymous_background$RNA_depth,
                final_synonymous_background$rna_vaf_max_weight,
                final_synonymous_background$func)
      )
    ) %>%
    dplyr::mutate(
      proportion_sampled = sampled_synon_mutations / total_synon_mutations,
      gene_tested = unique(gene_df$gene),
      synon_threshold = synon_threshold
    )

  return(list(
    final_synonymous_background = final_synonymous_background %>%
      dplyr::select(-pseudo_sample),
    sampling_summary_df = sampling_summary_df,
    proportion_sampled_per_patient = proportion_sampled_per_patient
  ))
}

#' @title Preprocess Mutation Table
#'
#' @description
#' Performs preprocessing on a mutation table. This includes:
#' \itemize{
#'   \item Removing samples with no expressed mutations.
#'   \item Excluding samples with more than a specified number of mutations per sample.
#'   \item Removing dinucleotide/trinucleotide (adjacent) mutations.
#'   \item Excluding mutations in a gene when a sample has more than a maximum allowed number.
#' }
#'
#' For each step the function records how many mutations and/or samples were removed.
#'
#' @param mutation_table A data frame or tibble containing the mutation data.
#'   Expected to have at least the columns \code{patient_id}, \code{RNA_alt_count}, \code{func},
#'   \code{gene}, \code{RNA_depth}, \code{pos}, and \code{chrom}.
#' @param non_synonymous_key A character vector of mutation types considered as non-synonymous.
#' @param depth_filt Numeric. The minimum RNA depth at a mutated position for the mutation to be considered
#' @param num_muts_filt Numeric. Minimum number of mutations (per gene) for the gene to be considered by RVdriver
#' @param max_coding_muts_per_sample Numeric. Maximum number of mutations allowed per sample.
#' @param max_muts_per_gene_per_sample Numeric. Maximum number of mutations allowed per gene for each sample.
#'
#' @return Processed mutation table:
#' \describe{
#'   \item{\code{mutation_table}}{The filtered mutation table.}
#' }
#' @export
preprocess_mutation_table <- function(mutation_table,
                                      non_synonymous_key,
                                      depth_filt,
                                      num_muts_filt,
                                      max_coding_muts_per_sample,
                                      max_muts_per_gene_per_sample) {
  summary_list <- list()

  # Record initial counts
  init_samples <- length(unique(mutation_table$patient_id))
  init_mutations <- nrow(mutation_table)
  summary_list$initial <- list(samples = init_samples, mutations = init_mutations)

  #### Step 1: Remove samples with no expressed mutations - low quality ####
  # Calculate the proportion of expressed mutations (RNA_alt_count > 2) per sample.
  prop_expressed <- mutation_table %>%
    group_by(patient_id) %>%
    summarise(prop_expressed = sum(RNA_alt_count > 2) / n(),
              num_expressed = sum(RNA_alt_count > 2),
              num = n(),
              .groups = "drop") %>%
    filter(num_expressed > 0)

  before_expressed <- nrow(mutation_table)
  mutation_table <- mutation_table %>%
    filter(patient_id %in% prop_expressed$patient_id)
  after_expressed <- nrow(mutation_table)
  post_samples <- length(unique(mutation_table$patient_id))

  summary_list$expressed_filter <- list(
    removed_mutations = before_expressed - after_expressed,
    removed_samples = init_samples - post_samples
  )

  #### Step 2: Define genes of interest ####
  genes <- mutation_table %>%
    filter(func %in% non_synonymous_key) %>%
    select(patient_id, gene, RNA_depth) %>%
    distinct() %>%  # do not count duplicate patient-gene pairs
    group_by(gene) %>%
    filter(RNA_depth >= depth_filt) %>%
    distinct(patient_id, gene) %>%
    summarise(num_muts = n(), .groups = "drop") %>%
    filter(num_muts > num_muts_filt) %>%
    select(gene, num_muts_above_depth_filter = num_muts)

  genes_of_interest <- genes$gene
  summary_list$genes_of_interest_count <- length(genes_of_interest)

  #### Step 3: Exclude samples exceeding the maximum coding mutations per sample ####
  nsampl <- sort(table(mutation_table$patient_id))
  to_exclude <- names(nsampl[nsampl > max_coding_muts_per_sample])
  before_coding <- nrow(mutation_table)
  mutation_table <- mutation_table[!(mutation_table$patient_id %in% to_exclude), ]
  after_coding <- nrow(mutation_table)

  summary_list$max_coding_filter <- list(
    removed_samples = length(to_exclude),
    removed_mutations = before_coding - after_coding
  )

  #### Step 4: Remove di/trinucleotide mutations ####
  # Create a unique key for each mutation.
  mutation_table$key <- paste(mutation_table$patient_id, mutation_table$gene, mutation_table$pos, sep = "_")

  # Identify multi-nucleotide events
  multi_muts <- mutation_table %>%
    group_by(patient_id, gene) %>%
    mutate(num_muts = n()) %>%
    filter(num_muts > 1) %>%
    mutate(pos = as.numeric(pos)) %>%
    arrange(patient_id, gene, pos, chrom) %>%
    mutate(
      diff_pos = c(Inf, diff(pos)),
      group_start = cumsum(diff_pos > 1)
    ) %>%
    filter(diff_pos == 1 | lead(diff_pos) == 1) %>%
    group_by(patient_id, gene, group_start) %>%
    mutate(group_key = paste(patient_id, gene, min(pos), sep = "_"),
           key = paste(patient_id, gene, pos, sep = "_")) %>%
    ungroup()

  before_multi <- nrow(mutation_table)
  mutation_table <- mutation_table %>%
    filter(!(key %in% multi_muts$key))
  after_multi <- nrow(mutation_table)

  summary_list$dinucleotide_filter <- list(
    removed_mutations = before_multi - after_multi
  )

  #### Step 5: Limit mutaations per gene per sample ####
  # Rank mutations within each patient-gene combination.
  mutrank <- ave(mutation_table$pos, paste(mutation_table$patient_id, mutation_table$gene), FUN = function(x) rank(x))
  before_gene_limit <- nrow(mutation_table)
  removed_gene_mut <- 0
  if (any(mutrank > max_muts_per_gene_per_sample)) {
    removed_gene_mut <- sum(mutrank > max_muts_per_gene_per_sample)
    mutation_table <- mutation_table[mutrank <= max_muts_per_gene_per_sample, ]
  }
  after_gene_limit <- nrow(mutation_table)

  summary_list$max_gene_filter <- list(
    removed_mutations = removed_gene_mut,
    total_removed = before_gene_limit - after_gene_limit
  )

  #### Summary Message ####
  message("Preprocessing Summary:")
  message(sprintf("Initial: %d mutations across %d samples", init_mutations, init_samples))
  message(sprintf("After expressed filter: removed %d mutations, %d samples removed",
                  summary_list$expressed_filter$removed_mutations,
                  summary_list$expressed_filter$removed_samples))
  message(sprintf("After max coding filter: removed %d mutations, %d samples removed",
                  summary_list$max_coding_filter$removed_mutations,
                  summary_list$max_coding_filter$removed_samples))
  message(sprintf("After dinucleotide/trinucleotide filter: removed %d mutations",
                  summary_list$dinucleotide_filter$removed_mutations))
  message(sprintf("After gene mutation limit filter: removed %d mutations",
                  summary_list$max_gene_filter$removed_mutations))
  message(sprintf("Number of genes of interest based on filter: %d", summary_list$genes_of_interest_count))

  return(list(mutation_table = mutation_table, genes_of_interest = genes_of_interest))
}


#' @title Plot NMD scaling and return the scaled RNA VAFs for nonsense mutations
#'
#' @description
#' This function performs NMD scaling on nonsense mutations by computing a scaling
#' factor as \code{mean(RNA_VAF[Missense]) / mean(RNA_VAF[Nonsense])} (forced to be at
#' least 1) and applies it to the nonsense mutations (with values capped at 1).
#' It creates a boxplot that compares:
#' \itemize{
#'   \item \strong{Nonsense (Pre-scaling)}: Original RNA VAF for nonsense mutations.
#'   \item \strong{Nonsense (Post-scaling)}: Scaled RNA VAF for nonsense mutations.
#'   \item \strong{Missense}: RNA VAF for missense mutations.
#' }
#'
#' @param mutation_table A data frame containing mutation data. Expected columns include
#'   \code{gene}, \code{func}, \code{RNA_VAF}, \code{RNA_depth}, \code{pos}, and optionally
#'   \code{chrom} and \code{is_cgc}.
#' @param cancer_type Character. Cancer type label used in plotting and saving
#' @param outdir Character. Directory in which to save the output plot.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{plot}: The ggplot object for the NMD scaling comparison.
#'   \item \code{mutation_table}: The updated mutation table with the added
#'         \code{RNA_VAF_scaled} column.
#' }
#' @export
plot_nmd_scaling_comparison <- function(mutation_table, cancer_type, outdir) {

  # Cap RNA_VAF at 1
  mutation_table <- mutation_table %>%
    mutate(RNA_VAF = ifelse(RNA_VAF > 1, 1, RNA_VAF))

  # Restrict analysis to nonsense and missense mutations
  data_subset <- mutation_table %>%
    filter(func %in% c("Nonsense_Mutation", "Missense_Mutation"))

  # Compute mean RNA_VAF for each type and derive scaling factor
  scaling_stats <- data_subset %>%
    group_by(func) %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = func, values_from = mean_RNA_VAF)

  nmd_scaling_factor <- scaling_stats$Missense_Mutation / scaling_stats$Nonsense_Mutation
  cat("NMD Scaling Factor for Nonsense Mutations: ", nmd_scaling_factor, "\n")

  if (is.na(nmd_scaling_factor) || nmd_scaling_factor < 1) {
    nmd_scaling_factor <- 1
  }

  # Update mutation_table: for nonsense, compute scaled value (capped at 1)
  mutation_table <- mutation_table %>%
    mutate(RNA_VAF_scaled = ifelse(func == "Nonsense_Mutation",
                                   pmin(RNA_VAF * nmd_scaling_factor, 1),
                                   RNA_VAF))

  # Prepare data for plotting
  nonsense_data <- data_subset %>%
    filter(func == "Nonsense_Mutation") %>%
    mutate(mutation_id = paste(gene, pos, sep = "_"),
           RNA_VAF_post = pmin(RNA_VAF * nmd_scaling_factor, 1))

  missense_data <- data_subset %>%
    filter(func == "Missense_Mutation") %>%
    mutate(mutation_id = paste(gene, pos, sep = "_"))

  nonsense_long <- nonsense_data %>%
    select(mutation_id, RNA_VAF, RNA_VAF_post) %>%
    pivot_longer(cols = c("RNA_VAF", "RNA_VAF_post"),
                 names_to = "scaling",
                 values_to = "value") %>%
    mutate(group = ifelse(scaling == "RNA_VAF", "Nonsense (Pre-scaling)", "Nonsense (Post-scaling)"))

  missense_long <- missense_data %>%
    mutate(group = "Missense") %>%
    select(mutation_id, value = RNA_VAF, group)

  plot_data <- bind_rows(nonsense_long, missense_long) %>%
    mutate(group = factor(group, levels = c("Nonsense (Pre-scaling)",
                                            "Nonsense (Post-scaling)",
                                            "Missense")))

  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    theme_bw() +
    labs(x = "Mutation Type", y = "RNA VAF",
         title = paste(cancer_type, "NMD Scaling Comparison"),
         subtitle = paste("Scaling Factor =", round(nmd_scaling_factor, 3))) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_jitter(aes(x = group, y = value),
                width = 0.15, size = 1, alpha = 0.7)

  # Save the plot to outdir
  ggsave(filename = file.path(outdir, paste0(cancer_type, "_nmd_scaling_comparison.png")),
         plot = p, width = 8, height = 6, dpi = 300)

  return(list(plot = p, mutation_table = mutation_table))
}

#' @title Compute RNA depth weights for all mutations
#'
#' @description
#' Computes mutation-specific weights based on RNA depth. The function:
#' \enumerate{
#'   \item Caps RNA depth at the 90th percentile of RNA depths across the cohort
#'   \item Adds a new column (\code{RNA_depth_capped}) for the capped depth.
#'   \item Generates and saves a histogram of the RNA depth distribution.
#'   \item Computes \code{log_cov} as \code{log2(RNA_depth_capped + 1)}.
#'   \item Computes weights (\code{rna_vaf_weight} and \code{rna_vaf_max_weight})
#'         using helper functions.
#'   \item Replaces zeros in \code{rna_vaf_max_weight} as if they had a coverage of 1
#' }
#'
#' @param mutation_table Input mutation table
#' @param cancer_type Character. Cancer type label used in plotting and saving
#' @param outdir Character. Directory in which to save the depth plot.
#'
#' @return The updated mutation table with added columns:
#'   \code{RNA_depth_capped}, \code{log_cov}, \code{rna_vaf_weight}, and \code{rna_vaf_max_weight}.
#' @export
compute_rna_depth_weights <- function(mutation_table, cancer_type, outdir) {

  # Cap RNA depth at the 90th percentile
  upper_quartile <- quantile(mutation_table$RNA_depth, 0.9)
  mutation_table <- mutation_table %>%
    mutate(RNA_depth_capped = ifelse(RNA_depth >= upper_quartile, upper_quartile, RNA_depth))

  # Create and save the RNA depth distribution plot
  depth_plot <- mutation_table %>%
    ggplot(aes(x = RNA_depth)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "RNA Depth",
         y = "Number of Mutations",
         title = paste0(cancer_type, " RNA Depth Distribution")) +
    scale_y_log10() +
    scale_x_log10()
  ggsave(filename = file.path(outdir, paste0(cancer_type, "_rna_depth_distribution.png")),
         plot = depth_plot, width = 8, height = 6, dpi = 300)

  # Compute log_cov and weights
  mutation_table <- mutation_table %>%
    ungroup() %>%
    mutate(log_cov = log2(RNA_depth_capped + 1),
           rna_vaf_weight = get_rna_weights(log_cov),
           rna_vaf_max_weight = get_rna_max_weights(log_cov))
  rna_weight_at_1 <- log2(2) / max(mutation_table$log_cov)
  mutation_table$rna_vaf_max_weight[mutation_table$rna_vaf_max_weight == 0] <- rna_weight_at_1

  return(mutation_table)
}

#' @title Perform synonymous mutation scaling and output plots
#'
#' @description
#' Performs scaling of synonymous mutations to match the average RNA VAFs of missense mutations
#'
#' @param synonymous_background A data frame containing the synonymous mutations
#' @param mutation_table A data frame with all mutations.
#' @param non_synonymous_key A character vector of non-synonymous mutation types.
#' @param genes A data frame of genes of interest
#' @param cgc_list A data frame of CGC genes.
#' @param synon_depth_filt Numeric. The depth threshold for synonymous mutations.
#' @param cancer_type Character. Cancer type label
#' @param outdir Character. Directory in which to save the scaling plots.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{plot}: A combined ggplot object comparing results before and after scaling.
#'   \item \code{synonymous_background}: The updated synonymous background table with scaled \code{RNA_VAF}.
#' }
#' @export
perform_synonymous_scaling <- function(synonymous_background,
                                       mutation_table,
                                       non_synonymous_key,
                                       genes,
                                       cgc_list,
                                       synon_depth_filt,
                                       cancer_type,
                                       outdir) {

  ############################
  # Pre-scaling summary
  ############################
  synon_scaling <- synonymous_background %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF, na.rm = TRUE),
              median_RNA_VAF = median(RNA_VAF, na.rm = TRUE),
              num_muts = n(),
              .groups = "drop")
  synon_scaling_all <- synonymous_background %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF, na.rm = TRUE),
              median_RNA_VAF = median(RNA_VAF, na.rm = TRUE),
              num_muts = n())
  gene_non_synon <- mutation_table %>%
    filter(func %in% non_synonymous_key) %>%
    filter(gene %in% genes$gene) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene, patient_id) %>%
    summarise(mean_rna_vaf = mean(RNA_VAF, na.rm = TRUE), .groups = "drop") %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf, na.rm = TRUE),
              median_RNA_VAF = median(mean_rna_vaf, na.rm = TRUE),
              num_muts = n(), .groups = "drop")
  gene_non_synon_all <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene) %>%
    summarise(mean_rna_vaf = mean(RNA_VAF, na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf, na.rm = TRUE),
              median_RNA_VAF = median(mean_rna_vaf, na.rm = TRUE),
              num_muts = n())

  combined_data <- bind_rows(
    synon_scaling %>% mutate(mutation_type = "Synonymous"),
    gene_non_synon %>% mutate(mutation_type = "Non-Synonymous")
  )

  comparison_plot_patient <- combined_data %>%
    group_by(patient_id) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    ggplot(aes(x = mutation_type, y = mean_RNA_VAF)) +
    geom_line(aes(group = patient_id), color = "gray", linetype = "dashed", alpha = 0.6) +
    geom_boxplot(aes(fill = mutation_type), alpha = 0.5) +
    geom_sina(aes(size = num_muts), alpha = 0.7) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "RNA VAF per Patient\nfor Synonymous vs\nNon-Synonymous Mutations\n(Before Scaling)",
         x = "Mutation Type", y = "RNA VAF") +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       label = "p.signif", label.y = max(combined_data$mean_RNA_VAF, na.rm = TRUE) + 0.05)

  gene_non_synon_w_cgc <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    group_by(gene) %>%
    summarise(mean_rna_vaf = mean(RNA_VAF, na.rm = TRUE), .groups = "drop") %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf, na.rm = TRUE),
              median_RNA_VAF = median(mean_rna_vaf, na.rm = TRUE),
              num_muts = n())
  gene_non_synon_all <- rbind(gene_non_synon_all, gene_non_synon_w_cgc)
  gene_non_synon_all$include_cgc <- c(FALSE, TRUE)
  gene_non_synon_all$cancer_type <- cancer_type

  target_mean <- gene_non_synon_all$mean_RNA_VAF[1]
  target_median <- gene_non_synon_all$median_RNA_VAF[1]
  original_mean <- synon_scaling_all$mean_RNA_VAF[1]
  original_median <- synon_scaling_all$median_RNA_VAF[1]

  find_scaling_factor <- function(k) {
    scaled_mean <- original_mean * k
    scaled_median <- original_median * k
    mean_diff <- abs(scaled_mean - target_mean)
    median_diff <- abs(scaled_median - target_median)
    return(max(mean_diff, median_diff))
  }
  result <- optimize(f = find_scaling_factor, interval = c(0, 2), maximum = FALSE)
  optimal_k <- result$minimum

  cat("Scaled mean:", original_mean * optimal_k, "\n")
  cat("Scaled median:", original_median * optimal_k, "\n")
  cat("Original mean:", original_mean, "\n")
  cat("Original median:", original_median, "\n")
  cat("Target mean:", target_mean, "\n")
  cat("Target median:", target_median, "\n")

  # Apply scaling to synonymous_background
  scaling_to_choose <- optimal_k
  message(sprintf("Scaling Factor applied to synonymous mutations is %s", scaling_to_choose))
  if(scaling_to_choose < 1){
    message(sprintf("WARNING: synonymous mutation scaling factor: %s, is < 1, results should be checked carefully, consider capping at 1", scaling_to_choose))
  }
  synonymous_background <- synonymous_background %>%
    group_by(patient_id) %>%
    mutate(RNA_VAF = RNA_VAF * scaling_to_choose) %>%
    ungroup() %>%
    mutate(RNA_VAF = ifelse(RNA_VAF > 1, 1, RNA_VAF))

  synon_scaling_after <- synonymous_background %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF, na.rm = TRUE),
              median_RNA_VAF = median(RNA_VAF, na.rm = TRUE),
              num_muts = n(), .groups = "drop")
  combined_data_alt <- bind_rows(
    synon_scaling_after %>% mutate(mutation_type = "Synonymous"),
    gene_non_synon %>% mutate(mutation_type = "Non-Synonymous")
  )

  comparison_plot_patient_post <- combined_data_alt %>%
    group_by(patient_id) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    ggplot(aes(x = mutation_type, y = mean_RNA_VAF)) +
    geom_line(aes(group = patient_id), color = "gray", linetype = "dashed", alpha = 0.6) +
    geom_boxplot(aes(fill = mutation_type), alpha = 0.5) +
    geom_sina(aes(size = num_muts), alpha = 0.7) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "RNA VAF per Patient\nfor Synonymous vs\nNon-Synonymous Mutations\n(After Scaling)",
         x = "Mutation Type", y = "RNA VAF") +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       label = "p.signif",
                       label.y = max(combined_data_alt$mean_RNA_VAF, na.rm = TRUE) + 0.05)

  comp_plot_all <- comparison_plot_patient + comparison_plot_patient_post

  # plot comparing non-synonymous vs. silent mutations
  silent <- synonymous_background %>%
    filter(RNA_depth >= synon_depth_filt) %>%
    filter(func == "Silent") %>%
    select(gene, func, RNA_VAF)
  non_synon_in_genes <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    mutate(func = "non_synonymous") %>%
    select(gene, func, RNA_VAF) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene, func) %>%
    summarise(RNA_VAF = mean(RNA_VAF, na.rm = TRUE), .groups = "drop")
  non_synon_in_genes <- bind_rows(non_synon_in_genes, silent)

  comparison_plot <- non_synon_in_genes %>%
    mutate(func = ifelse(func == "Silent", "synonymous", "non_synonymous")) %>%
    ggplot(aes(x = func, y = RNA_VAF)) +
    geom_sina(alpha = 0.3) +
    geom_boxplot() +
    theme_bw() +
    stat_compare_means() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., digits = 2)),
                 size = 6, vjust = 1) +
    stat_summary(fun.y = mean, geom = "text", aes(label = round(..y.., digits = 2)),
                 size = 6, vjust = -1, color = "red") +
    stat_summary(fun.y = mean, geom = "point", size = 3,
                 color = "red", shape = 18) +
    labs(title = paste0(cancer_type, " RNA VAF across mutation categories"),
         subtitle = paste0("Median scaling factor = ", round(original_median * optimal_k, digits = 3),
                           "\nMean scaling factor = ", round(original_mean * optimal_k, digits = 3),
                           "\nOptimised scaling factor = ", round(optimal_k, digits = 3)))

  # Combine
  comp_plot_all2 <- (comparison_plot) / (comp_plot_all)

  # Save
  ggsave(filename = file.path(outdir, paste0(cancer_type, "_synonymous_scaling_comparison.png")),
         plot = comp_plot_all2, width = 10, height = 12, dpi = 300)

  return(list(plot = comp_plot_all2,
              synonymous_background = synonymous_background))
}


#' @title Calculate sampling rate based on cohort size
#'
#' @description
#' Computes a rate such that at 5% of the total cohort size, the desired number
#' of synonymous mutations is achieved. The rate is chosen as the minimum of the rate
#' computed from the desired number of synonymous mutations and 0.202 which ensures at least
#' one synonymous mutation is sampled per patient when only 2 nonsynonymous mutations are present
#'
#' @param total_cohort Numeric. Total number of samples in the cohort.
#' @param synonymous_sampled Numeric. Desired number of synonymous mutations sampled
#'   at 5% of the cohort (10).
#' @param max_rate Numeric. Maximum rate (0.202).
#'
#' @return A numeric rate.
#' @export
calculate_rate <- function(total_cohort, synonymous_sampled = 10, max_rate = 0.202) {
  # Calculate 5% of the total cohort size
  five_percent_cohort <- 0.05 * total_cohort

  rate_general <- log(synonymous_sampled) / five_percent_cohort

  # Return the minimum between the computed rate and max_rate.
  rate <- min(rate_general, max_rate)
  return(rate)
}

#' @title Exponential growth sampling capped at 2x sample size
#'
#' @description
#' Computes the number of synonymous mutations to be sampled per test.
#'
#' @param n_samples Numeric. Number of samples.
#' @param base_num Numeric. The base number for the exponential calculation (default is 1).
#' @param rate Numeric. The growth rate.
#'
#' @return An integer sample size.
#' @export
exponential_growth_limited_2x <- function(n_samples, base_num = 1, rate) {
  # Compute the exponential growth result, rounded to nearest integer.
  result <- round(base_num * exp(rate * n_samples))

  # Limit the result to be no more than twice the number of samples.
  result <- pmin(result, 2 * n_samples)

  # Ensure that if the result is less than 1, it is set to 1. - this should never be the case
  result[result < 1] <- 1

  return(result)
}

#' @title Plot QQ Plot for RVdriver Results
#'
#' @description
#' Generates a QQ plot comparing the expected and observed -log10 p-values for RVdriver results
#' The plot is saved to the output directory.
#'
#' @param data A data frame containing the RVdriver results
#' @param method_name Character. The name of the method to plot.
#' @param outdir Character. Directory in which to save the QQ plot.
#'
#' @return QQ plot
#' @export
plot_qq <- function(data, method_name, outdir) {
  library(dplyr)
  library(ggplot2)

  data <- data %>%
    filter(method == method_name) %>%
    filter(!is.na(pval) & !is.na(p.adj))
  n <- nrow(data)
  data <- data %>%
    mutate(pval = ifelse(pval < 1e-6, 1e-6, pval)) %>%
    arrange(pval) %>%
    mutate(expected = -log10(ppoints(n)),
           observed = -log10(pval))
  lambda <- median(data$observed) / -log10(0.5)
  p <- ggplot(data, aes(x = expected, y = observed)) +
    geom_point(aes(color = ifelse(p.adj < 0.25 & is_cgc, "CGC Gene",
                                  ifelse(p.adj < 0.25, "Significant", "Not Significant"))), size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_text(data = subset(data, p.adj < 0.25),
              aes(label = gene),
              vjust = -1, hjust = -0.5, size = 3, check_overlap = TRUE) +
    labs(title = paste("QQ Plot for", method_name),
         x = "-log10(Expected P-value)",
         y = "-log10(Observed P-value)") +
    scale_color_manual(values = c("CGC Gene" = "red", "Significant" = "blue", "Not Significant" = "gray")) +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(filename = file.path(outdir, paste0(method_name, "_qq_plot.png")),
         plot = p, width = 8, height = 6, dpi = 300)

  return(p)
}

