# run_RVdriver_example.R
# ------------------------------------------------------------------------------
# This script demonstrates how to run the complete RVdriver pipeline.
#
# It performs the following steps:
#   1. Reads the input mutation table, CGC list, and seed list.
#   2. Preprocesses the mutation table.
#   3. Performs NMD scaling on nonsense mutations.
#   4. Computes RNA depth weights.
#   5. Performs synonymous (silent) mutation scaling.
#   6. Builds gene-level RVdriver input objects.
#   7. Runs gene-wise RVdriver analysis over multiple iterations.
#   8. Writes output CSV files and saves output plots to the specified output directory.
#
# Make sure your package (and all functions such as preprocess_mutation_table,
# plot_nmd_scaling_comparison, compute_rna_depth_weights, perform_synonymous_scaling,
# exponential_growth_limited_2x, calculate_rate, rvdriver, and plot_qq) are installed.
# ------------------------------------------------------------------------------
library(RVdriver)




# --- PARAMETERS ---
cancer_type       <- "UCS"                              # Change to your cancer type label
depth_filt        <- 8                                  # Minimum RNA depth for non-synonymous mutations
synon_depth_filt  <- 8                                  # Minimum RNA depth for synonymous mutations
num_muts_filt     <- 1                                  # Minimum number of mutations per gene
ncores            <- 1                                  # Number of cores for parallel processing
outdir            <- paste0("example_data/",
                            cancer_type, "_results/")   # Output directory (ensure it exists or will be created)
silent_scaling    <- TRUE                               # Whether to perform silent (synonymous) scaling
sampling_method   <- "relaxed"                          # Sampling method ("relaxed" or "strict")
include_gene_synons <- TRUE                             # Whether to include gene-specific synonymous mutations

# File paths for inputs 
mutation_path <- system.file("example_data", "UCS_mutation_table.rds", package = "RVdriver")
cgc_path      <- system.file("assets", "cgc_list.csv", package = "RVdriver")                # Shipped with the package
seed_path     <- system.file("assets", "seed_list.txt", package = "RVdriver")                # Shipped with the package


# Print parameters
cat("Running RVdriver on", cancer_type, "with the following parameters:\n")
cat("Depth filter:", depth_filt, "\n")
cat("Synon depth filter:", synon_depth_filt, "\n")
cat("Number of mutations filter:", num_muts_filt, "\n")
cat("Number of cores:", ncores, "\n")
cat("Output directory:", outdir, "\n")
cat("Silent scaling:", silent_scaling, "\n")
cat("Sampling method:", sampling_method, "\n")
cat("Include gene synonymous mutations:", include_gene_synons, "\n")

# Set seed
set.seed(999)

### Read input tables
if(grepl("rds$", mutation_path)) {
  mutation_table <- readRDS(mutation_path)
} else {
  mutation_table <- read_tsv(mutation_path)
}

if("is_cgc" %in% colnames(mutation_table)){
  mutation_table <- mutation_table %>%
    select(-is_cgc)
}

# Read CGC list
cgc_list <- read.csv(cgc_path) %>%
  mutate(is_cgc = TRUE)

# Merge CGC info into the mutation table
mutation_table <- mutation_table %>%
  left_join(cgc_list, by = "gene")

mutation_table$is_cgc[is.na(mutation_table$is_cgc)] <- FALSE

# Define non-synonymous mutation keys
non_synonymous_key <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")

### Preprocessing
processed_mutations <- preprocess_mutation_table(mutation_table = mutation_table,
                                                 non_synonymous_key = non_synonymous_key,
                                                 depth_filt = depth_filt,
                                                 num_muts_filt = num_muts_filt,
                                                 max_coding_muts_per_sample = 3000,
                                                 max_muts_per_gene_per_sample = 3)

mutation_table <- processed_mutations$mutation_table
# Assume that preprocess_mutation_table() returns genes_of_interest in its summary
genes_of_interest <- processed_mutations$genes_of_interest

### NMD Scaling
scaled_vafs <- plot_nmd_scaling_comparison(mutation_table = mutation_table,
                                           cancer_type = cancer_type,
                                           outdir = outdir)

mutation_table <- scaled_vafs$mutation_table
# If a new column RNA_VAF_scaled was created, override RNA_VAF with it.
if("RNA_VAF_scaled" %in% names(mutation_table)){
  mutation_table$RNA_VAF <- mutation_table$RNA_VAF_scaled
}

### Compute mutation weights
mutation_table <- compute_rna_depth_weights(mutation_table, cancer_type, outdir)

### RVdriver input objects
synonymous_background <- mutation_table %>%
  filter(func == "Silent", RNA_depth >= synon_depth_filt) %>%
  select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func,
         RNA_ref_count, RNA_alt_count)

patient_non_synon_list <- lapply(unique(mutation_table$patient_id), function(p) {
  patient_df <- mutation_table %>%
    filter(patient_id == p, func %in% non_synonymous_key) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func,
           RNA_ref_count, RNA_alt_count)
  list(patient_df = patient_df)
})
names(patient_non_synon_list) <- unique(mutation_table$patient_id)

### Genes RVdriver will test
genes <- mutation_table %>%
  filter(func %in% non_synonymous_key) %>%
  select(patient_id, gene, RNA_depth) %>%
  distinct() %>%
  group_by(gene) %>%
  filter(RNA_depth >= depth_filt) %>%
  distinct(patient_id, gene) %>%
  summarise(num_muts = n(), .groups = "drop") %>%
  filter(num_muts > num_muts_filt) %>%
  select(gene, num_muts_above_depth_filter = num_muts)

### Synonymous mutation scaling
scaled_synonymous <- perform_synonymous_scaling(synonymous_background = synonymous_background,
                                                mutation_table = mutation_table,
                                                non_synonymous_key = non_synonymous_key,
                                                genes = genes,
                                                cgc_list = cgc_list,
                                                synon_depth_filt = synon_depth_filt,
                                                cancer_type = cancer_type,
                                                outdir = outdir)

synonymous_background <- scaled_synonymous$synonymous_background

### Another RVdriver input object
gene_tables_list <- lapply(unique(genes$gene), function(g) {
  gene_df <- mutation_table %>%
    filter(gene == g, func %in% non_synonymous_key) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func)
  synonymous_background_filt <- synonymous_background %>%
    filter(patient_id %in% gene_df$patient_id) %>%
    filter(RNA_depth >= synon_depth_filt) %>%
    select(patient_id, RNA_VAF, gene, RNA_depth, rna_vaf_max_weight, func) %>%
    mutate(func = "synonymous")
  list(gene_df = gene_df, synonymous_background_filt = synonymous_background_filt)
})
names(gene_tables_list) <- unique(genes$gene)

### Seed list
seed_list <- read.table(seed_path, header = T) %>%
  pull(seed)

### Calculate the sampling rqate based on the size of the cohort
total_cohort <- length(unique(mutation_table$patient_id))
rate <- calculate_rate(total_cohort = total_cohort)

# Run RV
results_table <- lapply(
  seq_along(gene_tables_list),
  FUN = function(g_index,
                 gene_tables_list,
                 seed_list,
                 patient_non_synon_list,
                 sampling_method,
                 ncores,
                 include_gene_synons) {
    cat("Testing gene", g_index, "out of", nrow(genes), "\n")
    gene_name <- names(gene_tables_list)[g_index]
    synonymous_background_filt <- gene_tables_list[[gene_name]]$synonymous_background_filt
    if (nrow(synonymous_background_filt) == 0 && sampling_method == "strict") {
      return(NULL)
    }
    gene_df <- gene_tables_list[[gene_name]]$gene_df
    comp_df <- bind_rows(gene_df, synonymous_background_filt)
    num_samples <- length(unique(comp_df$patient_id))
    num_to_sample <- exponential_growth_limited_2x(num_samples, rate = rate)
    res_expo <- rvdriver(comp_df,
                         seed_list,
                         synon_threshold = num_to_sample,
                         num_iter = 25,
                         g = gene_name,
                         ncores = ncores,
                         non_synon_muts_list = patient_non_synon_list,
                         sampling_method = sampling_method,
                         include_gene_synons = include_gene_synons,
                         rounding_method = "expo")
    proportion_sampled_per_patient_expo <- res_expo$proportion_sampled_per_patient %>%
      mutate(sampling_func = "expo")
    res_expo <- res_expo$res_df %>% mutate(sampling_func = "expo")

    # Paste output to terminal per gene
    print(res_expo %>% select(gene, method, pval, number_synonymous_filtered, num_mutations))
    return(list(res = res_expo,
                proportion_sampled_per_patient = proportion_sampled_per_patient_expo))
  },
  gene_tables_list = gene_tables_list,
  seed_list = seed_list,
  patient_non_synon_list = patient_non_synon_list,
  sampling_method = sampling_method,
  ncores = ncores,
  include_gene_synons = include_gene_synons
)

# Combine gene-level results
res_list <- lapply(results_table, function(x) x$res)
proportion_sampled_per_patient_list <- lapply(results_table, function(x) x$proportion_sampled_per_patient)

res_df <- bind_rows(res_list)
proportion_sampled_per_patient_df <- bind_rows(proportion_sampled_per_patient_list) %>% distinct()

### Save output
if(!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
write.csv(proportion_sampled_per_patient_df,
          file = file.path(outdir, paste0(cancer_type, "_proportion_sampled_per_patient.csv")),
          row.names = FALSE)
res_df <- res_df %>%
  left_join(cgc_list, by = "gene")
res_df$is_cgc[is.na(res_df$is_cgc)] <- FALSE
write.csv(res_df,
          file = file.path(outdir, paste0(cancer_type, "_RV_results_all_mutations.csv")),
          row.names = FALSE)

# Process results table for genes with multiple tests - we take the lowest p-value for such cases
res_df <- res_df %>%
  group_by(method, gene) %>%
  filter(pval == min(pval)) %>%
  ungroup() %>%
  mutate(canc_type = cancer_type)
filtered_results <- res_df %>%
  group_by(gene, method) %>%
  filter(n() > 1) %>%
  slice_max(num_mutations, with_ties = FALSE) %>%
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
write.csv(res_df,
          file = file.path(outdir, paste0(cancer_type, "_RV_results_most_sig_func.csv")),
          row.names = FALSE)

cat("RVdriver analysis complete.\n")
plot_qq(res_df, method_name = "weighted_RVdriver", outdir = outdir)
plot_qq(res_df, method_name = "non_weighted_RVdriver", outdir = outdir)
