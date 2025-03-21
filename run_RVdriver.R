cat("Loading required libraries\n")

library(tidyverse)
library(lmerTest)
library(matrixStats)
library(parallel)
library(argparse)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggforce)
library(gridExtra)
library(patchwork)

############## MAIN #####################

# source RVdriver functions
source("bin/rvdriver_functions.R")

parser <- ArgumentParser(description = "Run RVdriver")

parser$add_argument("-i", "--input", help = "input mutations file", required = TRUE)
parser$add_argument("-c", "--cancer_type", help = "Cancer type to run RVdriver on")
parser$add_argument("-d", "--depth_filt", help = "Minimum depth to filter non-synonymous mutations (default = 8)", default = 8)
parser$add_argument("-s", "--synon_depth_filt", help = "Minimum depth to filter synonymous mutations on, default = 8", default = 8)
parser$add_argument("-n", "--num_muts_filt", help = "Minimum number of non-synonymous mutations required to include a gene", default = 1)
parser$add_argument("-p", "--ncores", help = "Number of cores to use")
parser$add_argument("-o", "--outdir", help = "Output directory")
parser$add_argument("-l", "--silent_scaling", help = "Whether to scale silent mutations", default = TRUE)
parser$add_argument("-sm", "--sampling_method", help = "Sampling method to use", default = "relaxed")
parser$add_argument("-g", "--include_gene_synons", help = "include Gene synonymous mutations?", default = TRUE)

args <- parser$parse_args()

if(!file.exists(args$input)) {
  stop("Input mutation table does not exist")
}

# if the suffix is rds - read it in as an rds file, otherwise read it in as a tsv

if(grepl("rds$", args$input)) {
  mutation_table <- readRDS(args$input)
} else {
  mutation_table <- read_tsv(args$input)
}
cancer_type <- args$cancer_type
depth_filt <- as.numeric(args$depth_filt)
synon_depth_filt <- as.numeric(args$synon_depth_filt)
num_muts_filt <- as.numeric(args$num_muts_filt)
ncores <- as.numeric(args$ncores)
outdir <- args$outdir
silent_scaling <- as.logical(args$silent_scaling)
sampling_method <- args$sampling_method
include_gene_synons <- as.logical(args$include_gene_synons)

cat("Running RVdriver on ", cancer_type, " with the following parameters: \n")
cat("Depth filter: ", depth_filt, "\n")
cat("Synon depth filter: ", synon_depth_filt, "\n")
cat("Number of mutations filter: ", num_muts_filt, "\n")
cat("Number of cores: ", ncores, "\n")
cat("Output directory: ", outdir, "\n")
cat("Silent scaling: ", silent_scaling, "\n")
cat("Sampling method: ", sampling_method, "\n")
cat("Include gene synonymous mutations: ", include_gene_synons, "\n")

# Default parameters
max_coding_muts_per_sample <- 3000
max_muts_per_gene_per_sample <- 3

outdir <- paste0(outdir, "/", cancer_type, "/RVdriver_results")
set.seed(999)

# make the output directory if it doesn't exist
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

pdf(paste0(outdir, "/", cancer_type, "_RV_plots.pdf"), width = 15, height = 10)

cat("Reading in CGC list ... \n")
# add cgc list
cgc_list <- read_csv("./assets/cgc_list.csv", col_types = cols()) %>%
  mutate(is_cgc = TRUE)

mutation_table <- mutation_table %>%
  left_join(cgc_list)

mutation_table$is_cgc[is.na(mutation_table$is_cgc)] <- FALSE

# Non-synonymous mutation key
non_synonymous_key <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")

# STEP 1 - Preprocessing...
# remove samples with no expressed mutations
# exclude samples with > 3000 mutations - as in dNdScv
# remove dinucleotide and trinucleotide mutations
# exclude samples with > 3 mutations in a single gene - as in dNdScv

# filter out cases where no mutations are expressed 
prop_expressed <- mutation_table %>%
  group_by(patient_id) %>%
  summarise(prop_expressed = sum(RNA_alt_count > 2) / n(),
            num_expressed = sum(RNA_alt_count > 2),
            num = n()) %>%
  filter(num_expressed > 0)

mutation_table <- mutation_table %>%
  filter(patient_id %in% prop_expressed$patient_id)

genes <- mutation_table %>%
  filter(func %in% non_synonymous_key) %>%
  select(patient_id, gene, RNA_depth) %>% unique() %>% # dont count more than one mutation in a single sample here
  group_by(gene) %>%
  filter(RNA_depth >= depth_filt) %>%
  select(patient_id, gene) %>% unique() %>%
  dplyr::summarise(num_muts = n()) %>%
  filter(num_muts > num_muts_filt) %>%
  select(gene, num_muts_above_depth_filter = num_muts)

genes_of_interest <- genes$gene

# Excluding samples exceeding the limit of mutations/sample - as in dNdScv
nsampl = sort(table(mutation_table$patient_id))
exclsamples = NULL
if (any(nsampl>max_coding_muts_per_sample)) {
  message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample (see the max_coding_muts_per_sample argument in dndscv). %0.0f samples left after filtering.',sum(nsampl>max_coding_muts_per_sample),sum(nsampl<=max_coding_muts_per_sample)))
  exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
  mutation_table = mutation_table[!(mutation_table$patient_id %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
}

mutation_table$key <- paste(mutation_table$patient_id, mutation_table$gene, mutation_table$pos, sep = "_")

# Detect dinucleotide and trinucleotide substitutions and create a unique identifier for them
multi_muts <- mutation_table %>%
  group_by(patient_id, gene) %>%
  mutate(num_muts = n()) %>%
  filter(num_muts > 1) %>%
  # Convert pos to numeric and arrange by position and chromosome
  mutate(pos = as.numeric(pos)) %>%
  arrange(patient_id, gene, pos, chrom) %>%
  # Create a group when the difference between consecutive positions is exactly 1
  mutate(
    diff_pos = c(Inf, diff(pos)),
    group_start = cumsum(diff_pos > 1)
  ) %>%
  # Only keep groups where mutations are truly adjacent
  filter(diff_pos == 1 | lead(diff_pos) == 1) %>%
  # Create the group_key using the patient_id, gene, and first position in the group
  group_by(patient_id, gene, group_start) %>%
  mutate(group_key = paste(patient_id, gene, min(pos), sep = "_"),
         key = paste(patient_id, gene, pos, sep = "_")) %>%
  ungroup()

# filter out the multi muts from the main table
mutation_table <- mutation_table %>%
  filter(!(key %in% multi_muts$key))

# Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) - as in dNdScv
mutrank = ave(mutation_table$pos, paste(mutation_table$patient_id,mutation_table$gene), FUN = function(x) rank(x))
exclmuts = NULL
if (any(mutrank>max_muts_per_gene_per_sample)) {
  message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per geneper sample (see the max_muts_per_gene_per_sample argument)',sum(mutrank>max_muts_per_gene_per_sample)))
  exclmuts = mutation_table[mutrank>max_muts_per_gene_per_sample,]
  mutation_table = mutation_table[mutrank<=max_muts_per_gene_per_sample,]
}


# STEP 2: NMD scaling

# NMD scaling of individual nonsense mutations may have caused RNA VAF to exceed 1, so we need to cap it at 1
mutation_table$RNA_VAF[mutation_table$RNA_VAF > 1] <- 1

# cancer type RNA VAF plots
# 1. non cgc RNA VAF across mutation categories
vaf_missense_comp_plot <- mutation_table %>%
  filter((func %in% non_synonymous_key & gene %in% genes_of_interest) | func == "Silent") %>%
  filter(func %in% c("Missense_Mutation", "Silent")) %>%
  # remove any silent mutations with RNA depth < synon_depth_filt
  filter(!(func == "Silent" & RNA_depth < synon_depth_filt)) %>%
  ggplot(aes(x = func, y = RNA_VAF, fill = func)) +
  geom_boxplot() +
  theme_minimal() +
  stat_compare_means() +
  facet_wrap(~is_cgc) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Mutation Category", y = "RNA VAF", title = paste0(cancer_type, " RNA VAF across mutation categories"))

# get the scaling factor for nonsense mutations - we use the difference between the median RNA VAF of missense and nonsense mutations
nmd_scaling_factor <- mutation_table %>%
  filter(func %in% c("Nonsense_Mutation", "Missense_Mutation")) %>%
  # filter(RNA_depth > 7) %>% # filter out low coverage mutations
  group_by(func) %>%
  dplyr::summarise(median_RNA_VAF = mean(RNA_VAF)) %>%
  spread(func, median_RNA_VAF) %>%
  mutate(scaling_factor = Missense_Mutation / Nonsense_Mutation) %>%
  select(scaling_factor) %>%
  as.numeric() 

# we don't want to scale nonsense mutations down - only up to account for NMD
if(nmd_scaling_factor < 1) {
  nmd_scaling_factor <- 1
}

mutation_table$RNA_VAF[mutation_table$func == "Nonsense_Mutation"] <- mutation_table$RNA_VAF[mutation_table$func == "Nonsense_Mutation"] * nmd_scaling_factor

# This again might lead to cases where RNA VAF exceeds 1
mutation_table$RNA_VAF[mutation_table$RNA_VAF > 1] <- 1

#redo the mutation vaf plot after scaling the nonsense mutations
vaf_nonsense_comp_plot_post_scaling <- mutation_table %>%
  filter((func %in% non_synonymous_key & gene %in% genes_of_interest) | func == "Silent") %>%
  filter(func %in% c("Nonsense_Mutation", "Silent")) %>%
  filter(!(func == "Silent" & RNA_depth < synon_depth_filt)) %>%
  ggplot(aes(x = func, y = RNA_VAF, fill = func)) +
  geom_boxplot() +
  theme_minimal() +
  stat_compare_means() +
  facet_wrap(~is_cgc) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Mutation Category", y = "RNA VAF", title = paste0(cancer_type, " RNA VAF across mutation categories"))

grid.arrange(vaf_missense_comp_plot, vaf_nonsense_comp_plot_post_scaling, ncol = 2)


# STEP 3 - generate weights per mutation - based on RNA depth at the mutation site

# get the upper quartile of RNA depth and cap it at the 90th percentile
upper_quartile <- quantile(mutation_table$RNA_depth, 0.9)
mutation_table$RNA_depth_capped <- mutation_table$RNA_depth
mutation_table$RNA_depth_capped[mutation_table$RNA_depth_capped >= upper_quartile] <- upper_quartile

# plot distribution of rna depths for all mutations
mutation_table %>%
  ggplot(aes(x = RNA_depth)) +
  geom_histogram(binwidth = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "RNA depth", y = "Number of mutations", title = paste0(cancer_type, " RNA depth distribution")) +
  scale_y_log10() + 
  scale_x_log10()

# get the weights for RNA VAF - based on coverage at the mutation site
mutation_table <- mutation_table %>%
  ungroup() %>%
  mutate(log_cov = log2((RNA_depth_capped + 1))) %>%
  mutate(rna_vaf_weight = get_rna_weights(log_cov)) %>%
  mutate(rna_vaf_max_weight = get_rna_max_weights(log_cov))

rna_weight_at_1 <- log2(2)/max(mutation_table$log_cov)

# replace any 0 weights with the weight at 1 - otherwise they will be excluded from the analysis
mutation_table$rna_vaf_max_weight[mutation_table$rna_vaf_max_weight == 0] <- rna_weight_at_1


# STEP 4 - create the required input objects for RVdriver

# get the synonymous mutations
synonymous_background <- mutation_table %>%
  filter(func == "Silent") %>%
  filter(RNA_depth >= synon_depth_filt) %>%
  select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func, RNA_ref_count, RNA_alt_count)

# create a list of patient dataframes for non-synonymous mutations
patient_non_synon_list <- lapply(unique(mutation_table$patient_id), function(p) {
  
  patient_df <- mutation_table %>%
    filter(patient_id == p, func %in% non_synonymous_key) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func, RNA_ref_count, RNA_alt_count)
  
  list(patient_df = patient_df)
  
})

names(patient_non_synon_list) <- unique(mutation_table$patient_id)


genes <- mutation_table %>%
  filter(func %in% non_synonymous_key) %>%
  select(patient_id, gene, RNA_depth) %>% unique() %>% # dont count more than one mutation in a single sample here
  group_by(gene) %>%
  filter(RNA_depth >= depth_filt) %>%
  select(patient_id, gene) %>% unique() %>%
  dplyr::summarise(num_muts = n()) %>%
  filter(num_muts > num_muts_filt) %>%
  select(gene, num_muts_above_depth_filter = num_muts)


# STEP 5 - Perform synonymous mutation scaling
if(silent_scaling) {
  
  
  synon_scaling <- synonymous_background %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF),
              median_RNA_VAF = median(RNA_VAF),
              num_muts = n()) 
  
  synon_scaling_all <- synonymous_background %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF),
              median_RNA_VAF = median(RNA_VAF),
              num_muts = n())
  
  gene_non_synon <- mutation_table %>%
    filter(func %in% non_synonymous_key) %>%
    filter(gene %in% genes$gene) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene, patient_id) %>%
    summarise(mean_rna_vaf = (RNA_VAF)) %>%
    ungroup() %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf), 
              median_RNA_VAF = median(mean_rna_vaf),
              num_muts = n()) 
  
  # comparing agains missense mutations in non-cgc genes
  gene_non_synon_all <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene) %>%
    summarise(mean_rna_vaf = mean(RNA_VAF)) %>%
    ungroup() %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf), 
              median_RNA_VAF = median(mean_rna_vaf),
              num_muts = n()) 
  
  # Combine the data frames
  combined_data <- bind_rows(
    synon_scaling %>% mutate(mutation_type = "Synonymous"),
    gene_non_synon %>% mutate(mutation_type = "Non-Synonymous")
  )
  
  # Plot the boxplot with lines and perform paired Wilcoxon test
  comparison_plot_patient <- combined_data %>%
    # filter(num_muts > 5) %>%
    group_by(patient_id) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    ggplot(aes(x = mutation_type, y = mean_RNA_VAF)) +
    geom_line(aes(group = patient_id), color = "gray", linetype = "dashed", alpha = 0.6) +
    geom_boxplot(aes(fill = mutation_type), alpha = 0.5) + 
    geom_sina(alpha = 0.7, aes(size = num_muts)) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "RNA VAF per Patient for Synonymous vs Non-Synonymous Mutations",
         x = "Mutation Type", y = "RNA VAF") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, 
                       label = "p.signif", label.y = max(combined_data$mean_RNA_VAF) + 0.05)
  
  gene_non_synon_w_cgc <-  mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    # filter(gene %in% cgc_list$gene) %>%
    group_by(gene) %>%
    summarise(mean_rna_vaf = mean(RNA_VAF)) %>%
    ungroup() %>%
    # group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(mean_rna_vaf), 
              median_RNA_VAF = median(mean_rna_vaf),
              num_muts = n())
  
  gene_non_synon_all <- rbind(gene_non_synon_all, gene_non_synon_w_cgc)
  gene_non_synon_all$include_cgc <- c(FALSE, TRUE)
  gene_non_synon_all$cancer_type <- cancer_type
  
  write.csv(gene_non_synon_all, paste0(outdir, "/", cancer_type, "_gene_non_synon_w_cgc_comp.csv"), quote = FALSE, row.names = FALSE)
  
  
  # Define the target mean and median from gene_non_synon_all
  target_mean <- gene_non_synon_all$mean_RNA_VAF[1] # Assuming we use the first row for target
  target_median <- gene_non_synon_all$median_RNA_VAF[1] # Assuming we use the first row for target
  
  # Define the original mean and median from synon_scaling
  original_mean <- synon_scaling_all$mean_RNA_VAF[1]
  original_median <- synon_scaling_all$median_RNA_VAF[1]
  
  find_scaling_factor <- function(k) {
    scaled_mean <- original_mean * k
    scaled_median <- original_median * k
    
    mean_diff <- abs(scaled_mean - target_mean)
    median_diff <- abs(scaled_median - target_median)
    
    # Minimize the maximum of the differences to balance both mean and median scaling
    return(max(mean_diff, median_diff))
  }
  
  result <- optimize(f = find_scaling_factor, interval = c(0, 2), maximum = FALSE)
  
  # Get the optimal scaling factor
  optimal_k <- result$minimum
  
  # Calculate the scaled mean and median
  scaled_mean <- original_mean * optimal_k
  scaled_median <- original_median * optimal_k
  
  # Output the scaling factor and the adjusted means and medians
  cat("Optimal scaling factor:", optimal_k, "\n")
  cat("Scaled mean:", scaled_mean, "\n")
  cat("Scaled median:", scaled_median, "\n")
  cat("original mean:", original_mean, "\n")
  cat("original median:", original_median, "\n")
  cat("Target mean:", target_mean, "\n")
  cat("Target median:", target_median, "\n")
  
  mean_scaling_factor <- gene_non_synon_all$mean_RNA_VAF[1] / synon_scaling_all$mean_RNA_VAF
  median_scaling_factor <- gene_non_synon_all$median_RNA_VAF[1] / synon_scaling_all$median_RNA_VAF
  
  # Calculate the combined scaling factor as an average of mean and median scaling factors
  combined_scaling_factor <- (mean_scaling_factor + median_scaling_factor) / 2
  
  non_synon_in_genes <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    mutate(func = "non_synonymous") %>%
    select(gene, func, RNA_VAF) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene, func) %>%
    summarise(RNA_VAF = mean(RNA_VAF))
  
  silent <- synonymous_background %>%
    filter(RNA_depth >= synon_depth_filt) %>%
    filter(func == "Silent")
  
  silent <- silent %>%
    select(gene, func, RNA_VAF) 
  
  non_synon_in_genes <- rbind(non_synon_in_genes, silent)
  
  medians_per_func <- non_synon_in_genes %>%
    group_by(func) %>%
    summarise(median_RNA_VAF = median(RNA_VAF)) %>%
    ungroup()
  
  comparison_plot_pre <- non_synon_in_genes %>%
    mutate(func = ifelse(func == "Silent", "synonymous", "non_synonymous")) %>%
    ggplot(aes(x = func, y = RNA_VAF)) +
    geom_sina(alpha = 0.3) + 
    geom_boxplot() +
    theme_bw() +
    stat_compare_means() + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    # Label the median value for the boxplot
    stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., digits = 2)), size = 6, vjust = 1) +
    # Label the mean value for the boxplot
    stat_summary(fun.y = mean, geom = "text", aes(label = round(..y.., digits = 2)), size = 6, vjust = -1, color = "red") +
    # Plot the mean as a point
    stat_summary(fun.y = mean, geom = "point", size = 3, color = "red", shape = 18) +
    # Add median scaling factor to the title and the mean scaling factor to the subtitle
    labs(title = paste0(cancer_type, " RNA VAF across mutation categories"), 
         subtitle = paste0("Median scaling factor = ", round(median_scaling_factor, digits = 3),
                           "\nMean scaling factor = ", round(mean_scaling_factor, digits = 3),
                           "\nCombined scaling factor = ", round(combined_scaling_factor, digits = 3),
                           "\nOptimised scaling factor = ", round(optimal_k, digits = 3)))

  scaling_to_choose <- optimal_k

  message(sprintf("    Note: Scaling Factor being applied to synonymous mutations is %s", optimal_k))  
  
  # apply mean scaling factor to synonymous mutations
  for(sample_to_check in unique(synonymous_background$patient_id)){
    synonymous_background$RNA_VAF[synonymous_background$patient_id == sample_to_check] <- synonymous_background$RNA_VAF[synonymous_background$patient_id == sample_to_check] * scaling_to_choose
  }
  
  # cap RNA VAF at 1
  synonymous_background$RNA_VAF[synonymous_background$RNA_VAF > 1] <- 1
  
  # REMAKE THE ABOVE PLOTS AND PUT SIDE BY SIDE
  synon_scaling_after <- synonymous_background %>%
    group_by(patient_id) %>%
    summarise(mean_RNA_VAF = mean(RNA_VAF),
              median_RNA_VAF = median(RNA_VAF),
              num_muts = n()) 
  
  # Step 1: Combine the data frames
  combined_data_alt <- bind_rows(
    synon_scaling_after %>% mutate(mutation_type = "Synonymous"),
    gene_non_synon %>% mutate(mutation_type = "Non-Synonymous")
  )
  
  # Step 2: Plot the boxplot with lines and perform paired Wilcoxon test
  comparison_plot_patient_post <- combined_data_alt %>%
    # filter(num_muts > 5) %>%
    group_by(patient_id) %>%
    filter(n() == 2) %>%
    ungroup() %>%
    ggplot(aes(x = mutation_type, y = mean_RNA_VAF)) +
    geom_line(aes(group = patient_id), color = "gray", linetype = "dashed", alpha = 0.6) +
    geom_boxplot(aes(fill = mutation_type), alpha = 0.5) + 
    geom_sina(alpha = 0.7, aes(size = num_muts)) +
    theme_bw() +
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "RNA VAF per Patient for Synonymous vs Non-Synonymous Mutations",
         x = "Mutation Type", y = "RNA VAF") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, 
                       label = "p.signif", label.y = max(combined_data_alt$mean_RNA_VAF) + 0.05)
  
  # remake the plot with the scaling factor applied
  to_plot <- comparison_plot_patient + comparison_plot_patient_post
  print(to_plot)
  
  silent <- synonymous_background %>%
    filter(RNA_depth >= synon_depth_filt) %>%
    filter(func == "Silent")
  
  silent <- silent %>%
    select(gene, func, RNA_VAF) 
  
  non_synon_in_genes <- mutation_table %>%
    filter(func == "Missense_Mutation") %>%
    filter(gene %in% genes$gene) %>%
    mutate(func = "non_synonymous") %>%
    select(gene, func, RNA_VAF) %>%
    filter(!gene %in% cgc_list$gene) %>%
    group_by(gene, func) %>%
    summarise(RNA_VAF = mean(RNA_VAF))
  
  non_synon_in_genes <- rbind(non_synon_in_genes, silent)
  

  # create the comparison plot after scaling
  comparison_plot <- non_synon_in_genes %>%
    mutate(func = ifelse(func == "Silent", "synonymous", "non_synonymous")) %>%
    ggplot(aes(x = func, y = RNA_VAF)) +
    geom_sina(alpha = 0.3) + 
    geom_boxplot() +
    theme_bw() +
    stat_compare_means() + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    # Label the median value for the boxplot
    stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., digits = 2)), size = 6, vjust = 1) +
    # Label the mean value for the boxplot
    stat_summary(fun.y = mean, geom = "text", aes(label = round(..y.., digits = 2)), size = 6, vjust = -1, color = "red") +
    # Plot the mean as a point
    stat_summary(fun.y = mean, geom = "point", size = 3, color = "red", shape = 18) +
    # Add median scaling factor to the title and the mean scaling factor to the subtitle
    labs(title = paste0(cancer_type, " RNA VAF across mutation categories"), 
         subtitle = paste0("Median scaling factor = ", round(median_scaling_factor, digits = 3),
                           "\nMean scaling factor = ", round(mean_scaling_factor, digits = 3),
                           "\nCombined scaling factor = ", round(combined_scaling_factor, digits = 3),
                           "\nOptimised scaling factor = ", round(optimal_k, digits = 3)))
  
  comp_plot_all <- comparison_plot_pre + comparison_plot
  
  comp_plot_all
  
}



# List to store the tables for each gene
gene_tables_list <- lapply(unique(genes$gene), function(g) {
  
  gene_df <- mutation_table %>%
    filter(gene == g, func %in% non_synonymous_key) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func)
  
  synonymous_background_filt <- synonymous_background %>%
    filter(patient_id %in% gene_df$patient_id) %>%
    filter(RNA_depth >= synon_depth_filt) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth,rna_vaf_max_weight, func) %>%
    mutate(func = "synonymous")
  
  list(gene_df = gene_df, synonymous_background_filt = synonymous_background_filt)
})

names(gene_tables_list) <- unique(genes$gene)

seed_list <- read_tsv("assets/seed_list.txt", col_types = cols()) %>%
  pull(seed)


# Function to calculate the rate based on the total cohort size
calculate_rate <- function(total_cohort, synonymous_sampled = 10, max_rate = 0.202) {
  five_percent_cohort <- 0.05 * total_cohort
  
  # Calculate the rate to ensure that at 5% of the cohort, we get the desired_result (10)
  rate_general = log(synonymous_sampled) / five_percent_cohort
  
  # Return the minimum of the two rates to satisfy both conditions
  rate <- min(rate_general, 0.202)
  
  return(rate)
}


# Function to calculate the exponential growth using the calculated rate
exponential_growth_limited_2x <- function(n_samples, base_num = 1, rate) {
  # Calculate the result using the given rate
  result <- round(base_num * exp(rate * n_samples))
  
  # Limit the result to be no more than 2x the number of samples
  result <- pmin(result, 2 * n_samples)
  
  # Ensure that if the result is less than 1, it is set to 1
  result[result < 1] <- 1
  
  return(result)
}


total_cohort <- length(unique(mutation_table$patient_id))
rate <- calculate_rate(total_cohort = total_cohort)


results_table <- lapply(1:length(gene_tables_list),
                        FUN = function(g, gene_tables_list, seed_list, patient_non_synon_list, sampling_method, ncores, include_gene_synons) {
                          lm_summary <- data.frame()
                          msg <- paste0("Testing gene ", g, " out of ", nrow(genes), "\n")
                          cat(msg)
                          g <- names(gene_tables_list)[g]
                          synonymous_background_filt <- gene_tables_list[[g]]$synonymous_background_filt
                          
                          if(nrow(synonymous_background_filt) == 0 && sampling_method == "strict"){
                            return(NULL)
                          }

                          gene_df <- gene_tables_list[[g]]$gene_df
                          
                          comp_df <- rbind(gene_df, synonymous_background_filt)
                          
                          num_samples <- length(unique(comp_df$patient_id))
                          
                          num_to_sample <- exponential_growth_limited_2x(num_samples, rate = rate)
                          
                          res_expo <- rvdriver(comp_df, seed_list, synon_threshold = num_to_sample, 25,
                                                g, ncores, patient_non_synon_list, sampling_method,
                                                include_gene_synons, "expo")
                          
                          proportion_sampled_per_patient_expo <- res_expo$proportion_sampled_per_patient %>%
                            mutate(sampling_func = "expo")
                          
                          res_expo <- res_expo$res_df %>%
                            mutate(sampling_func = "expo")
                          
                          proportion_sampled_per_patient <- proportion_sampled_per_patient_expo
                          
                          res <- res_expo
                          # print res to terminal
                          print(as.data.frame(res %>% select(gene, method, pval, number_synonymous_filtered, num_mutations)))
                          return(list(res = res, proportion_sampled_per_patient = proportion_sampled_per_patient))
                        },
                        gene_tables_list = gene_tables_list,
                        seed_list = seed_list,
                        patient_non_synon_list = patient_non_synon_list,
                        sampling_method = sampling_method,
                        ncores = ncores,
                        include_gene_synons = include_gene_synons)


# Combine the results into a single data frame
res_list <- list()
proportion_sampled_per_patient_list <- list()

# Loop through each element in results_table
for (i in seq_along(results_table)) {
  res_list[[i]] <- results_table[[i]]$res
  proportion_sampled_per_patient_list[[i]] <- results_table[[i]]$proportion_sampled_per_patient
}

# Combine the data frames 
res_df <- bind_rows(res_list)
proportion_sampled_per_patient_df <- bind_rows(proportion_sampled_per_patient_list)

write.csv(proportion_sampled_per_patient_df, paste0(outdir, "/", cancer_type, "_proportion_sampled_per_patient.csv"), row.names = FALSE)

# add cgc info
res_df <- res_df %>%
  left_join(cgc_list, by = "gene")

res_df$is_cgc[is.na(res_df$is_cgc)] <- FALSE

# write the results to a csv file
write.csv(res_df, paste0(outdir, "/", cancer_type, "_RV_results_all_mutations.csv"), row.names = FALSE)



# filter for only the lowest p-value for each gene
res_df <- res_df %>%
  group_by(method, gene) %>%
  # filter(num_mutations > 1) %>%
  # only keep the best p-value for each gene
  filter(pval == min(pval)) %>%
  ungroup() %>%
  mutate(canc_type = cancer_type)

# sometimes there can still be ties - take the category with the most mutations
filtered_results <- res_df %>%
  group_by(gene, method) %>%
  filter(n() > 1) %>%  # Filter groups with more than one row
  slice_max(num_mutations, with_ties = FALSE) %>%  # Pick the row with the largest num_mutations
  ungroup() %>%
  mutate(method_gene = paste0(method, gene)) 

res_df <- res_df %>%
  mutate(method_gene = paste0(method, gene)) %>%
  filter(!method_gene %in% filtered_results$method_gene) %>%
  bind_rows(filtered_results) %>%
  select(-method_gene) %>%
  group_by(method) %>%
  mutate(p.adj = p.adjust(pval, "BH")) %>%
  arrange(p.adj) 

write.csv(res_df, paste0(outdir, "/", cancer_type, "_RV_results_most_sig_func.csv"), row.names = FALSE)

plot_qq <- function(data, method_name) {
  # Cap p-values and remove NA values
  data <- data %>%
    filter(method == method_name) %>%
    filter(!is.na(pval) & !is.na(p.adj))
  
  # Calculate expected p-values
  n <- nrow(data)
  
  # Calculate lambda
  
  data <- data %>%
    filter(method == method_name) %>%
    mutate(pval = ifelse(pval < 1e-6, 1e-6, pval)) %>%
    filter(!is.na(pval) & !is.na(p.adj)) %>%
    arrange(pval) %>%
    mutate(expected = -log10(ppoints(n)),
           observed = -log10(pval))
  lambda <- median(data$observed) / -log10(0.5)
  
  
  # Plot
  p <- ggplot(data, aes(x = expected, y = observed)) +
    geom_point(aes(color = ifelse(p.adj < 0.25 & is_cgc, "CGC Gene",
                                  ifelse(p.adj < 0.25, "Significant", "Not Significant"))), size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_text(data = subset(data, p.adj < 0.25), 
              aes(label = gene), 
              vjust = -1, 
              hjust = -0.5, 
              size = 3, 
              check_overlap = TRUE) +
    labs(title = paste("QQ Plot for", method_name),
         x = "-log10(Expected P-value)",
         y = "-log10(Observed P-value)",
         subtitle = paste("Lambda =", round(lambda, 2))) +
    scale_color_manual(values = c("CGC Gene" = "red", "Significant" = "blue", "Not Significant" = "gray")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
}


# Plot QQ plots for each method
methods <- unique(res_df$method)
for (method in methods) {
  data_subset <- res_df %>% filter(method == method)
  plot_qq(data_subset, method)
}

# Close the PDF file
dev.off()

sig <- res_df %>%
  filter(p.adj < 0.25) %>%
  arrange(p.adj)

table(sig$is_cgc, sig$method)