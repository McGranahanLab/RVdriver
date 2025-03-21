############### FUNCTIONS ####################

get_rna_weights <- function(coverage){
	# get weights for each variant based on RNA coverage
	# coverage: RNA coverage
	# return: weights
	weighted_expr <- coverage / sum(coverage)
	return(weighted_expr)
}

get_rna_max_weights <- function(coverage){
	# get weights for each variant based on RNA coverage
	# coverage: RNA coverage
	# return: weights
	weighted_expr <- coverage / max(coverage)
	return(weighted_expr)
}

get_lm_summary <- function(lm, g) {
	# get p-value and coefficient from linear model
	# lm: linear model
	# return: data frame with p value and coefficient
	# Extract coefficient table
	coef_table <- summary(lm)$coefficients

	# If 'func' is a factor, ensure you’re extracting the correct t-value.
	# The second row might not always correspond to the level of interest.
	# Alternatively, specify the level by name, e.g., "funcLevelName".
	tval <- coef_table["func_newnon_synonymous", "t value"]
	pval <- coef_table["func_newnon_synonymous", "Pr(>|t|)"]

	# Check the direction of the effect (t-value) and adjust p-value for a one-sided test
	if(tval > 0) {
	pval <- pval / 2
	} else {
	pval <- 1 - (pval / 2)
	}

	return(data.frame(method = "RVdriver_lm",
	gene = g,
	df = NA,
	effect_size = tval,
	pval = pval))
}


get_results_df_lm <- function(comp_df, g, weighted = FALSE){
  rna_vaf_res <- myTryCatch(get_rna_vaf_res_lm(comp_df, g, weighted))
  if(is.null(rna_vaf_res$value)) {
    rna_vaf_res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA) 
  } 
  else if(!is.null(rna_vaf_res$warning)){
    if(rna_vaf_res$warning$message == "essentially perfect fit: summary may be unreliable") {
      rna_vaf_res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA) 
    }
  } else {
    rna_vaf_res <- rna_vaf_res$value
  }
  res_df <- do.call(rbind, list(rna_vaf_res))
}


get_rna_vaf_res_lm <- function(df, g, weighted){

   if(length(unique(df$cluster)) == 1){
        unique_cluster <- TRUE
    } else {
        unique_cluster <- FALSE
    }


  if(weighted){
        if(unique_cluster){
            model <- lm(RNA_VAF ~ func_new, weights = rna_vaf_max_weight, data = df)
        } else {
            # Note: In a simple linear model (lm), we generally can't include random effects like in lmer
            # If cluster effect is substantial, consider adding it as a fixed effect
            model <- lm(RNA_VAF ~ func_new , weights = rna_vaf_max_weight, data = df)
        }
    } else {
        if(unique_cluster){
            model <- lm(RNA_VAF ~ func_new, data = df)
        } else {
            model <- lm(RNA_VAF ~ func_new,  data = df)
        }
    }

    return(get_lm_summary(model, g))

}

# from https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}



# Define the function for parallel execution
parallel_function <- function(iter, gene_df, syn_df, synon_threshold, g, seed_list, weighted = FALSE, non_synon_muts_list, sampling_method = sampling_method, include_gene_synons, rounding_method) {
  
  # set seed for this iteration #
  set.seed(seed_list[iter])
  
  res_all_cats <- data.frame()
  proportion_sampled_per_patient <- data.frame()
  
  func_above_1 <- gene_df %>%
    group_by(func) %>%
    filter(RNA_depth >= depth_filt) %>%
    summarise(num = length(unique(patient_id))) %>%
    filter(num > num_muts_filt) 
  
  # loop through each of the types of func in gene df
  for(mut_func in unique(func_above_1$func)){
    mut_func_df <- gene_df %>%
      filter(func == mut_func)
    
    num_mutations <- nrow(mut_func_df)
    num_samples <- length(unique(mut_func_df$patient_id))
    synon_threshold_filt <- if(rounding_method == "ceiling") {
      sqrt_ceiling(num_samples)
    } else if(rounding_method == "floor") {
      sqrt_floor(num_samples)
    } else if (rounding_method == "round") {
      sqrt_round(num_samples)
    } else if (rounding_method == "expo") {
      exponential_growth_limited_2x(num_samples, rate = rate) 
    }
    num_muts_greater_than_7_depth <- nrow(mut_func_df %>% filter(RNA_depth > 7))
    proportion_greater_than_7_depth <- num_muts_greater_than_7_depth / num_mutations
    
    mut_func_df$func_new <- "non_synonymous"
    
    syn_df_filt <- syn_df %>%
      filter(patient_id %in% mut_func_df$patient_id)
    
    # Sample a specified number of synonymous mutations for each patient
    syn_df_sampled <- create_synonymous_background_including_mut_type(syn_df_filt, mut_func_df, synon_threshold_filt, synonymous_background, seed = seed_list[iter], non_synon_muts_list, sampling_method = sampling_method, include_gene_synons = include_gene_synons)
    
    proportion_sampled_per_patient_df <- syn_df_sampled$proportion_sampled_per_patient
    
    proportion_sampled_per_patient <- rbind(proportion_sampled_per_patient, proportion_sampled_per_patient_df)
    
    sampling_summary_df <- syn_df_sampled$sampling_summary_df
    
    syn_df_sampled <- syn_df_sampled$final_synonymous_background
      
      if(nrow(syn_df_sampled) == 0 ) {
        res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA)
      }
      
      else {  
      syn_df_sampled$func <- "synonymous"
      
      syn_df_sampled$func_new <- syn_df_sampled$func
      
      syn_df_sampled_filt <- syn_df_sampled 
      
      # Combine the data for the gene and the sampled synonymous mutations
      df_combined <- rbind(mut_func_df, syn_df_sampled_filt)
      
      df_combined$func_new <- factor(df_combined$func_new, levels = c("synonymous", unique(df_combined$func_new[df_combined$func_new != "synonymous"])))
      
      res_weighted <- get_results_df_lm(df_combined, g, weighted = TRUE) %>% 
        mutate(method = "weighted_RVdriver_lm")
      res_non_weighted <- get_results_df_lm(df_combined, g, weighted = FALSE)
      res <- rbind(res_weighted, res_non_weighted)
      
      # Run the analysis on the combined data
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
  
  # if there is more than one category of mutation, then run the analysis on all mutations
  if(length(unique(gene_df$func)) > 1) {
    # Combine the data for the gene and the sampled synonymous mutations
    
    mut_func_df <- gene_df %>%
      mutate(func_new = "non_synonymous")
    
    syn_df_filt <- syn_df %>%
      filter(patient_id %in% mut_func_df$patient_id)
    
    # Sample a specified number of synonymous mutations for each patient
    syn_df_sampled <- create_synonymous_background_including_mut_type(syn_df_filt, mut_func_df, synon_threshold, synonymous_background, seed = seed_list[iter], non_synon_muts_list, sampling_method = sampling_method, include_gene_synons = include_gene_synons)
    
    proportion_sampled_per_patient_df <- syn_df_sampled$proportion_sampled_per_patient
    
    proportion_sampled_per_patient <- rbind(proportion_sampled_per_patient, proportion_sampled_per_patient_df)
    
    sampling_summary_df <- syn_df_sampled$sampling_summary_df
    
    syn_df_sampled <- syn_df_sampled$final_synonymous_background

    if(nrow(syn_df_sampled) == 0 ) {
        res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA)
    }
      
    else {  
    
        syn_df_sampled$func <- "synonymous"
        
        syn_df_sampled$func_new <- syn_df_sampled$func
        
        num_mutations <- nrow(mut_func_df)
        num_samples <- length(unique(mut_func_df$patient_id))
        
        num_muts_greater_than_7_depth <- nrow(mut_func_df %>% filter(RNA_depth > 7))
        proportion_greater_than_7_depth <- num_muts_greater_than_7_depth / num_mutations
        
        df_combined <- rbind(mut_func_df, syn_df_sampled)
        
        df_combined$func_new <- factor(df_combined$func_new, levels = c("synonymous", unique(df_combined$func_new[df_combined$func_new != "synonymous"])))
        
        # if(all(df_combined$RNA_VAF == 0)){
        #   res <- data.frame(method = "RVdriver_lm", gene = g, df = NA, effect_size = NA, pval = NA)
        # }
        

        res_weighted <- get_results_df_lm(df_combined, g, weighted = TRUE) %>% 
            mutate(method = "weighted_RVdriver_lm")
        res_non_weighted <- get_results_df_lm(df_combined, g, weighted = FALSE)
        res <- rbind(res_weighted, res_non_weighted)
        
        # Run the analysis on the combined data
        res <- cbind(res, sampling_summary_df)
        # res_lm <- get_results_df_lm(df_combined, g, weighted)
        
        # res <- rbind(res, res_lm)
        # res <- res_lm
        
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



rvdriver <- function(comp_df, seed_list, synon_threshold, num_iter, g, ncores, non_synon_muts_list, sampling_method, include_gene_synons, rounding_method){
  # comp_df: dataframe of all comparisons
  # seed_list: list of seed genes
  # synon_threshold: threshold for synonymous genes
  # num_iter: number of iterations
  # gene: gene of interest
  # returns: dataframe of results
  
  # Get the synonymous and gene-specific dataframes once
  syn_df <- comp_df %>% filter(func == "synonymous")
  gene_df <- comp_df %>% filter(gene == g, func != "synonymous")
  
  num_samples <- length(unique(gene_df$patient_id))
  
  
  # If possible, parallelize these using mclapply
  results_list_all <- mclapply(1:num_iter, parallel_function,
                               gene_df = gene_df, syn_df = syn_df,
                               synon_threshold = synon_threshold, g = g,
                               weighted = FALSE, seed_list = seed_list, mc.cores = ncores, non_synon_muts_list,
                               sampling_method = sampling_method, include_gene_synons = include_gene_synons, rounding_method)
  
  results_list <- lapply(results_list_all, function(x) x$results)
  
  proportion_sampled_per_patient <- lapply(results_list_all, function(x) x$proportion_sampled)
  
 
  res_df <- do.call(rbind, results_list)
  
  res_df_lm <- res_df %>%
    filter(method == "RVdriver_lm")
  
  proportion_sampled_per_patient <- do.call(rbind, proportion_sampled_per_patient)
  proportion_sampled_per_patient <- unique(proportion_sampled_per_patient)
 
  # do the same for weighted results
  res_df_weighted_lm <- res_df %>%
    filter(method == "weighted_RVdriver_lm")
  
  # make summary df
  weighted_df <- res_df_lm %>%
    group_by(func) %>%
    summarise(
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
    ) %>%
    # mutate(total_proportion_sampled = proportion_sampled_from_patient + proportion_sampled_from_mutated_patients +
    #         proportion_sampled_from_global_synon + proportion_sampled_from_same_mut_type) %>%
    # ungroup() %>%
    select(method, gene, func, median_tval, median_df, pval, pval_sd, num_signif, num_mutations, num_mutated_samples, number_synonymous_filtered, num_muts_greater_than_7_depth, proportion_greater_than_7_depth,
           num_samples_w_more_than_required, num_samples_sampling_from_mutated_samples, num_samples_sampling_from_global_synon, num_samples_sampling_from_mut_type_within_patient, num_samples_sampling_incomplete) 
  
  weighted_df_lm <- res_df_weighted_lm %>%
    group_by(func) %>%
    summarise(
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
    ) %>%
    # mutate(total_proportion_sampled = proportion_sampled_from_patient + proportion_sampled_from_mutated_patients +
    #         proportion_sampled_from_global_synon + proportion_sampled_from_same_mut_type) %>%
    # ungroup() %>%
    select(method, gene, func, median_tval, median_df, pval, pval_sd, num_signif, num_mutations, num_mutated_samples, number_synonymous_filtered, num_muts_greater_than_7_depth, proportion_greater_than_7_depth,
           num_samples_w_more_than_required, num_samples_sampling_from_mutated_samples, num_samples_sampling_from_global_synon, num_samples_sampling_from_mut_type_within_patient, num_samples_sampling_incomplete) 
  
  #    if(all(weighted_df_lm$total_proportion_sampled != 1)){
  #     print(paste0("WARNING: total proportion sampled for ", g, " is not 1"))
  #    }
  
  # res_df <- rbind(weighted_df, non_weighted_df, weighted_df_lm, non_weighted_df_lm)
  res_df <- rbind(weighted_df, weighted_df_lm)
  return(list(res_df = res_df, proportion_sampled_per_patient = proportion_sampled_per_patient))
  
  
}

sample_and_update <- function(patient_ids, threshold, source_df_name, source_name, missing) {
    for (patient_id_to_check in patient_ids) {
        while (missing > 0 && nrow(get(source_df_name)) > 0) {
            sampled_priority <- get(source_df_name) %>%
                sample_n(min(missing, nrow(get(source_df_name))), replace = FALSE) %>%
                select(colnames(syn_df)) %>%
                mutate(pseudo_sample = patient_id_to_check)

            final_synonymous_background <<- rbind(final_synonymous_background, sampled_priority)
            updated_df <- get(source_df_name) %>% anti_join(sampled_priority %>% select(-pseudo_sample), by = colnames(get(source_df_name)))
            assign(source_df_name, updated_df, envir = .GlobalEnv)

            missing <- threshold - nrow(final_synonymous_background %>% filter(pseudo_sample == patient_id_to_check))
            sampling_summary_df[[source_name]] <<- sampling_summary_df[[source_name]] + 1
        }
    }
    return(missing)
}




create_synonymous_background_including_mut_type <- function(syn_df, gene_df, synon_threshold, synonymous_background_filtered, seed, non_synon_muts_list, sampling_method, include_gene_synons = FALSE) {
    set.seed(seed)

    syn_df_summary <- syn_df %>%
        group_by(patient_id) %>%
        summarise(synon_count = n()) %>%
        ungroup() %>%
        complete(patient_id = gene_df$patient_id, fill = list(synon_count = 0))

    patients_above <- syn_df_summary %>%
        filter(synon_count >= synon_threshold) %>%
        pull(patient_id)
    patients_below <- syn_df_summary %>%
        filter(synon_count > 0 & synon_count < synon_threshold) %>%
        pull(patient_id)
    patients_with_0 <- syn_df_summary %>%
        filter(synon_count == 0) %>%
        pull(patient_id)

    final_synonymous_background <- data.frame(matrix(ncol = ncol(syn_df) + 1, nrow = 0))
    colnames(final_synonymous_background) <- c(colnames(syn_df), "pseudo_sample") 

  # Match the class of columns in final_synonymous_background to those in syn_df
  for (col in colnames(syn_df)) {
    class(final_synonymous_background[[col]]) <- class(syn_df[[col]])
  }
  
  # Set the pseudo_sample column as character class
  final_synonymous_background$pseudo_sample <- character()
  
    sampling_summary_df <- data.frame(
        num_samples_w_more_than_required = length(patients_above),
        num_samples_sampling_from_mutated_samples = 0,
        num_samples_sampling_from_global_synon = 0,
        num_samples_sampling_from_mut_type_within_patient = 0,
        num_samples_sampling_incomplete = 0
    )

    synonymous_background_filtered <- synonymous_background_filtered %>%
        select(colnames(syn_df)) %>%
        mutate(func = "synonymous") %>%
        filter(!patient_id %in% unique(gene_df$patient_id))


    if(include_gene_synons){
        # # first sample any synonymous mutations from our gene of interest and remove them from the pool
        gene_synons <- syn_df %>% filter(gene == unique(gene_df$gene)) %>% 
            mutate(pseudo_sample = patient_id)

        gene_synons_all <- synonymous_background_filtered %>% filter(gene == unique(gene_df$gene)) %>% 
            mutate(pseudo_sample = patient_id)

        all_gene_synons <- rbind(gene_synons, gene_synons_all)

        remaining_synonymous <- syn_df %>% anti_join(all_gene_synons, by = colnames(syn_df))
        synonymous_background_filtered <- synonymous_background_filtered %>% anti_join(all_gene_synons, by = colnames(syn_df))

        final_synonymous_background <- rbind(final_synonymous_background, all_gene_synons)
    }

    if (length(patients_above) > 0) {
        # Handle patients_above separately
        for (patient_id_to_check in patients_above) {
            patient_syn_df <- syn_df %>% filter(patient_id == patient_id_to_check)
            sampled_synonymous <- patient_syn_df %>%
                sample_n(synon_threshold, replace = FALSE) %>%
                mutate(pseudo_sample = patient_id_to_check)
            final_synonymous_background <- rbind(final_synonymous_background, sampled_synonymous)
        }

        # Remaining synonymous mutations after sampling patients_above
        remaining_synonymous <- syn_df %>% anti_join(final_synonymous_background, by = colnames(syn_df))
    }

    # Handle patients_below separately initially
    for (patient_id_to_check in patients_below) {
        patient_syn_df <- syn_df %>% filter(patient_id == patient_id_to_check)
        sampled_synonymous <- patient_syn_df %>%
            mutate(pseudo_sample = patient_id_to_check)
        final_synonymous_background <- rbind(final_synonymous_background, sampled_synonymous)
    }

    # Remaining synonymous mutations after initial handling of patients_below
    remaining_synonymous <- syn_df %>% anti_join(final_synonymous_background, by = colnames(syn_df))

    # Define the sample_and_update function without using global environment
    sample_and_update <- function(patient_ids, threshold, source_df, source_name, missing) {
        for (patient_id_to_check in patient_ids) {
            while (missing > 0 && nrow(source_df) > 0) {
                source_df_before <- nrow(source_df)
                sampled_priority <- source_df %>%
                    sample_n(min(missing, nrow(source_df)), replace = FALSE) %>%
                    select(colnames(syn_df)) %>%
                    mutate(pseudo_sample = patient_id_to_check)

                final_synonymous_background <<- rbind(final_synonymous_background, sampled_priority)
                source_df <- source_df %>% anti_join(sampled_priority %>% select(-pseudo_sample), by = colnames(source_df))

                source_df_after <- nrow(source_df)
                if (source_df_after >= source_df_before) {
                    stop(paste("Error: Rows not removed from", source_name, "for patient", patient_id_to_check))
                }

                missing <- threshold - nrow(final_synonymous_background %>% filter(pseudo_sample == patient_id_to_check))
            }
        }
        return(list(missing = missing, source_df = source_df))
    }

    # Update sample and handle patients_below
    for (patient_id_to_check in patients_below) {
        num_sampled <- nrow(final_synonymous_background %>% filter(pseudo_sample == patient_id_to_check))
        missing <- synon_threshold - num_sampled

        if (sampling_method == "strict") {
            missing <- 0

            next
        }

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_mutated_samples <- sampling_summary_df$num_samples_sampling_from_mutated_samples + 1

            next
        }

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, remaining_synonymous, "remaining_synonymous", missing)
        missing <- result$missing
        remaining_synonymous <- result$source_df

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_mutated_samples <- sampling_summary_df$num_samples_sampling_from_mutated_samples + 1

            next
        }

        if (sampling_method == "strict") {
            missing <- 0

            next
        }

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, synonymous_background_filtered, "synonymous_background_filtered", missing)
        missing <- result$missing
        synonymous_background_filtered <- result$source_df

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_global_synon <- sampling_summary_df$num_samples_sampling_from_global_synon + 1

            next
        }

        mut_type_avail <- non_synon_muts_list[[patient_id_to_check]]$patient_df %>%
            filter(gene != unique(gene_df$gene), func %in% unique(gene_df$func[gene_df$patient_id == patient_id_to_check])) %>%
            filter(RNA_depth >= depth_filt) %>%
            select(colnames(syn_df))

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, mut_type_avail, "mut_type_avail", missing)
        missing <- result$missing
        mut_type_avail <- result$source_df

        if(missing == 0) {

            sampling_summary_df$num_samples_sampling_from_mut_type_within_patient <- sampling_summary_df$num_samples_sampling_from_mut_type_within_patient + 1

            next
        }
        if (missing > 0) {
            # if there are still missing mutations at this point - add 1 to the not enough mutations counter and finish here
            sampling_summary_df$num_samples_sampling_incomplete <- sampling_summary_df$num_samples_sampling_incomplete + 1
        }
    }

    for (patient_id_to_check in patients_with_0) {
        num_sampled <- nrow(final_synonymous_background %>% filter(pseudo_sample == patient_id_to_check))
        missing <- synon_threshold - num_sampled

        if (sampling_method == "strict") {
            missing <- 0

            next
        }

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_mutated_samples <- sampling_summary_df$num_samples_sampling_from_mutated_samples + 1

            next
        }

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, remaining_synonymous, "remaining_synonymous", missing)
        missing <- result$missing
        remaining_synonymous <- result$source_df

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_mutated_samples <- sampling_summary_df$num_samples_sampling_from_mutated_samples + 1

            next
        }

        if (sampling_method == "strict") {
            missing <- 0

            next
        }

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, synonymous_background_filtered, "synonymous_background_filtered", missing)
        missing <- result$missing
        synonymous_background_filtered <- result$source_df

        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_global_synon <- sampling_summary_df$num_samples_sampling_from_global_synon + 1

            next
        }

        mut_type_avail <- non_synon_muts_list[[patient_id_to_check]]$patient_df %>%
            filter(gene != unique(gene_df$gene), func %in% unique(gene_df$func[gene_df$patient_id == patient_id_to_check])) %>%
            filter(RNA_depth >= depth_filt) %>%
            select(colnames(syn_df))

        result <- sample_and_update(c(patient_id_to_check), synon_threshold, mut_type_avail, "mut_type_avail", missing)
        missing <- result$missing
        mut_type_avail <- result$source_df
        
        if(missing == 0) {
            sampling_summary_df$num_samples_sampling_from_mut_type_within_patient <- sampling_summary_df$num_samples_sampling_from_mut_type_within_patient + 1

            next
        }
        
        if (missing > 0) {
            # if there are still missing mutations at this point - add 1 to the not enough mutations counter and finish here
            sampling_summary_df$num_samples_sampling_incomplete <- sampling_summary_df$num_samples_sampling_incomplete + 1
        }
    }

    # get the proportion of synon mutations per sample that have been sampled
    proportion_sampled_per_patient <- syn_df %>%
        group_by(patient_id) %>%
        summarise(
            total_synon_mutations = n(),
            sampled_synon_mutations = sum(paste(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func) %in% 
                paste(final_synonymous_background$patient_id, final_synonymous_background$RNA_VAF, final_synonymous_background$gene, final_synonymous_background$RNA_depth, final_synonymous_background$rna_vaf_max_weight, final_synonymous_background$func))
        ) %>%
        mutate(proportion_sampled = sampled_synon_mutations / total_synon_mutations,
               gene_tested = unique(gene_df$gene),
               synon_threshold = synon_threshold)

    return(list(final_synonymous_background = final_synonymous_background %>% select(-pseudo_sample), 
                sampling_summary_df = sampling_summary_df, 
                proportion_sampled_per_patient = proportion_sampled_per_patient))
}