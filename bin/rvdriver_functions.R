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

rvdriver <- function(comp_df, synonymous_background_filt, seed_list, synon_threshold, num_tests, gene){
  all_iterations <- list()

  for(iteration in 1:num_tests){
    set.seed(seed_list$seed[iteration])

    # sample 10 synon from each region
    synonymous_background_filt_threshold <- c()
    for(i in unique(synonymous_background_filt$patient_id)){
        region_synon <- synonymous_background_filt %>%
            filter(patient_id == i)
        region_synon <- region_synon[sample(1:nrow(region_synon), replace = T, size = synon_threshold),]
        synonymous_background_filt_threshold <- rbind(synonymous_background_filt_threshold, region_synon)
    }
    
    number_of_samples_in_synon <- length(unique(synonymous_background_filt_threshold$patient_id))
    numbers_with_synon <-length(unique(synonymous_background_filt$patient_id))

  comp_df <- rbind(comp_df, synonymous_background_filt_threshold)

  comp_df$func <- factor(comp_df$func, levels = c("synonymous", "non_synonymous"))

  gene_res_df <- get_results_df(comp_df, gene)

  gene_res_df <- gene_res_df %>%
    mutate(num_w_synon = numbers_with_synon) %>%
    mutate(num_w_synon_test = number_of_samples_in_synon) %>%
    mutate(canc_type = i) %>%
    mutate(synon_threshold = synon_threshold) %>%
    mutate(iteration = iteration)

  all_iterations[[iteration]] <- gene_res_df
  }

  all_res <- do.call(rbind, all_iterations)

  geometric_mean_pval <- exp(mean(log(all_res$pval), na.rm = T))


  m <- all_res %>%
    group_by(method) %>%
    summarise(t_value = median(t_value), df = median(df)) %>%
    mutate(pval = pt(t_value, df, lower.tail = FALSE)) %>%
    mutate(geometric_mean_pval) %>%
    mutate(num_mutations = length(unique(comp_df$patient_id))) %>%
    mutate(gene = as.character(gene)) %>%
    mutate(number_synonymous_filtered = unique(synon_threshold)) %>%
    mutate(method = as.character(method)) %>% 
    select(method, gene, t_value,
           pval, geometric_mean_pval, num_mutations,
           number_synonymous_filtered)

  return(m) 
}

get_results_df <- function(comp_df, g){
  rna_vaf_res <- myTryCatch(get_rna_vaf_res(comp_df, g))
    if(is.null(rna_vaf_res$value)) {
    rna_vaf_res <- data.frame(method = "RVdriver", gene = g, df = NA, effect_size = NA, pval = NA)
  } else {
    rna_vaf_res <- rna_vaf_res$value
  }
  res_df <- do.call(rbind, list(rna_vaf_res))
}

get_rna_vaf_res <- function(comp_df, g){
  p <- summary(lmer(RNA_VAF~func+(1|patient_id), data = comp_df)) %>%
                `[[`("coefficients") %>%
                as.data.frame() %>% rownames_to_column()

              tvalue <- p$`t value`[2]
              df <- p$df[2]
              pval <- pt(tvalue, df, lower.tail = FALSE)
              
              res <- data.frame(method = "RVdriver",
                                         gene = g,
                                         df, 
                                         t_value = tvalue,
                                         pval = pval)
  return(res)
}

