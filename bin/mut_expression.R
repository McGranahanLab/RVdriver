library(deepSNV)
library(argparse)
library(tidyverse)

########################################## START ##########################################

parser <- ArgumentParser()

# Patient specific parameters

parser$add_argument('--mut_path',  nargs=1,
                    help='Path to mutation file',
                    required=TRUE)
parser$add_argument('--rna_bam',  nargs=1,
                    help='Path to RNA bam file',
                    required=TRUE)
parser$add_argument('--sample_id',  nargs=1,
                    help='RNA sample ID',
                    required=TRUE)
parser$add_argument('--patient_id',  nargs=1,
                    help='Patient ID',
                    required=TRUE)
parser$add_argument('--out_dir',  nargs=1,
                    help='output directory',
                    required=TRUE)

args <- parser$parse_args()

mut_path <- args$mut_path
rna_bam <- args$rna_bam
sample_id <- args$sample_id
pat_id <- args$patient_id
out_dir <- args$out_dir

print(args)

# get sample mutation file
muts_df <- read.csv(mut_path) %>%
    filter(patient_id == pat_id)

# add rna sample id
muts_df$rna_sample_id = as.character(sample_id)

# add bam path
muts_df$rna_bam_path = as.character(rna_bam)

# get RNA VAF
for(x_idx in 1:nrow(muts_df)){
    if(!is.na(rna_bam) & !is.na(muts_df[x_idx,"pos"])){
        print(x_idx)
        rna_reads_output <- as.data.frame(bam2R(file = rna_bam,
                                                chr = muts_df[x_idx,"chr"],
                                                start = muts_df[x_idx,"pos"],
                                                stop = muts_df[x_idx,"pos"],
                                                mask = 1280)) # mask removes dup reads
        rna_reads_output <- rna_reads_output[,c("A", "C", "G", "T", "a", "c", "g", "t")]
        depth_no_dup <- sum(rna_reads_output[1,muts_df[x_idx,"alt"]] + rna_reads_output[1,tolower(muts_df[x_idx,"alt"])] + rna_reads_output[1,muts_df[x_idx,"ref"]] + rna_reads_output[1,tolower(muts_df[x_idx,"ref"])])
        alt_count_no_dup <- sum(rna_reads_output[1,muts_df[x_idx,"alt"]] + rna_reads_output[1,tolower(muts_df[x_idx,"alt"])])   

        RNA_VAF_no_dup <- alt_count_no_dup / depth_no_dup
        muts_df$RNA_pos_depth[x_idx] = as.numeric(depth_no_dup)
        muts_df$RNA_alt_count[x_idx] = as.numeric(alt_count_no_dup)
        muts_df$RNA_ref_count[x_idx] = as.numeric(depth_no_dup) - as.numeric(alt_count_no_dup)
        muts_df$RNA_VAF[x_idx] = as.numeric(alt_count_no_dup) / as.numeric(depth_no_dup)

    }
}

#create output table
muts_df <- muts_df %>%
    mutate(CSRA = paste0(chr,":",pos,":",ref,":",alt)) %>%
    select(patient_id, sample_id, gene, CSRA, RNA_ref_count, RNA_alt_count, RNA_VAF, func, canc_type)

# write output
write.csv(muts_df, file = paste0(out_dir, "/", sample_id, "_RNA_VAF_mutable.csv"), quote = F, row.names = F)
