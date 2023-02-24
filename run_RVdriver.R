library(tidyverse)
library(lmerTest)
library(argparse)

########################################## START ##########################################

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--mut_path',  nargs=1,
                    help='Path to mutation file',
                    required=TRUE)
parser$add_argument('--canc_type',  nargs=1,
                    help='cancer type of interest',
                    required=TRUE)
parser$add_argument('--mutation_filter',  nargs=1,
                    help='number of required mutations for a gene to be tested',
                    required=FALSE, default=3, type = "double")
parser$add_argument('--out_dir',  nargs=1,
                    help='cancer type of interest',
                    required=FALSE, default = "RVdriver_results_table.csv")

args <- parser$parse_args()

# parse args
input_mutations_path <- args$mut_path
cancer_type <- args$canc_type
num_muts_filt <- as.numeric(args$mutation_filter)
outdir <- args$out_dir
depth_filt <- 8
seed_list <- read_tsv("assets/seed_list.txt")


print(args)
# set seed
set.seed(999)

# source functions required to run RVdriver
source("./bin/rvdriver_functions.R")

cat("Reading in mutation table ... ")
# load mutation table
rna_muts_table <- read_tsv(input_mutations_path) %>%
    filter(canc_type == cancer_type)

cat("Reading in CGC list ... ")
# add cgc list
cgc_list <- read_csv("./assets/cgc_list.csv") %>%
  mutate(is_cgc = TRUE)

rna_muts_table <- rna_muts_table %>% left_join(cgc_list)

rna_muts_table$is_cgc[is.na(rna_muts_table$is_cgc)] <- FALSE

# calculate rna depth
rna_muts_table$RNA_depth <- rna_muts_table$RNA_alt_count + rna_muts_table$RNA_ref_count

# Specific to TCGA - change so input is just synonymous vs nonsynonymous
non_synonymous_key <- c("Nonsense_Mutation", "Nonstop_Mutation",
                                        "Splice_Site", "Missense_Mutation")

# determine which genes to test
passed_genes <- rna_muts_table %>%
    filter(func %in% non_synonymous_key) %>%
    select(patient_id, gene, RNA_depth) %>% unique() %>% # dont count more than one mutation in a single sample here
    group_by(gene) %>%
    filter(RNA_depth >= depth_filt) %>%
    select(patient_id, gene) %>% unique() %>%
    summarise(num_muts = n()) %>%
    filter(num_muts > num_muts_filt) %>%
    select(gene, num_muts_above_depth_filter = num_muts) 
        
# Synonymous mutations
synonymous_background <- rna_muts_table %>% filter(func == "Silent")

# run RVdriver
results_table <- lapply(unique(passed_genes$gene), synonymous_background, FUN = function(g, synonymous_background){
  cat("Testing", g,"...\n")
    # get gene df
    gene_df <- rna_muts_table %>% filter(gene == g, func %in% non_synonymous_key) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth) %>%
      mutate(func = "non_synonymous") 

    # synonymous mutations of samples with mutations in GOI - remove low coverage mutations 
    synonymous_background_filt <- synonymous_background %>%
      filter(patient_id %in% gene_df$patient_id) %>%
      filter(RNA_depth > 7) %>%
      select(patient_id, RNA_VAF, gene, RNA_depth) %>%
      mutate(func = "synonymous") 

  # perform test
  res <- rvdriver(gene_df, synonymous_background_filt, seed_list, synon_threshold = 10, 50, g) 

    return(res)
        })

        # merge cgc list with results table and add gene info
        results_table <- do.call(rbind, results_table)

        results_table <- results_table %>% left_join(cgc_list)

        results_table$is_cgc[is.na(results_table$is_cgc)] <- FALSE

        results_table <- results_table %>%
          left_join(passed_genes)

        results_table <- results_table %>%
          mutate(depth_filter = depth_filt) %>%
          mutate(cancer_type) %>%
          mutate(prop_cohort_filt = num_muts_above_depth_filter / length(unique(rna_muts_table$patient_id))) %>%
          mutate(p.adj = p.adjust(pval, "BH")) %>%
          arrange(p.adj)

# make out dir if it doesnt exist
if(!dir.exists(outdir)){
  dir.create(outdir)
}

# write results 

write.csv(results_table, file = paste0(outdir,"/", cancer_type,"_rvdriver_results.csv"),        
          quote = F, row.names = F)