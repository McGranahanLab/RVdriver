# RVdriver
A statistical approach to detect cancer driver genes utilising the RNA Variant Allele Frequency of somatic mutations

## Inputs

RVdriver only requires a simple tab deliminted or .RDS file mutation table in the following format

patient_id | gene | CSRA | RNA_ref_count | RNA_alt_count | RNA_VAF | func | canc_type 
----|----|------|-----|------|------|------|------

an example mutation table is given [here](./assets/UVM_mutation_table.rds)

within the ```func``` column, the following mutation types are accepted:
```
Missense_Mutation
Nonsense_Mutation
Nonstop_Mutation
Splice_Site
Silent
```

### The following steps are needed to run RVdriver on the example dataset:
1. clone this repo
```
git clone https://github.com/McGranahanLab/RVdriver.git
cd RVdriver
```
if necessary, create the conda environment containing the required R packages (all that is required to run RVdriver is argparse, tidyverse and lmerTest)
```
conda create --name RVdriver r-tidyverse r-argparse r-lmertest
```
Activate the conda environment
```
conda activate RVdriver
```
OR using a singularity container
```
mkdir singularity_images
cd singularity_images/
singularity pull --arch amd64 library://tpjones15/default/rvdriver:latest
cd ../
```
2. run RVdriver
```
# for conda
Rscript run_RVdriver.R -i "assets/UVM_mutation_table.rds" -c UVM -d 8 -s 8 -n 1 -p 2 -o test -l TRUE -sm "relaxed" -g TRUE
```
```
# for singularity 
singularity exec -B ${PWD}:${PWD} singularity_images/rvdriver_latest.sif Rscript run_RVdriver.R "assets/UVM_mutation_table.rds" -c UVM -d 8 -s 8 -n 1 -p 2 -o test -l TRUE -sm "relaxed" -g TRUE
```
This will create a results table in the ```test/``` directory.  
These results should be the same as those presented [here](./test_data_results/UVM)

#### The following steps were taken to calculate the RNA VAFs for somatic mutations across the pan cancer TCGA dataset

1. Download the following files and store them in the assets directory. These files are required for GATK preprocessing (MarkDuplicates and BQSR):
    - [TCGA reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) 
    - [dbsnp_138_vcf](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf)
    - [dbsnp_138_vcf_idx](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx)

2. Pull the relevent singularity container (this singularity image can also be used to run RVdriver as above) 
```
mkdir singularity_images
cd singularity_images/
singularity pull --arch amd64 library://tpjones15/default/rvdriver_rnavaf:latest
```    
The mutation table can be a single table with all samples or multiple sample specific tables. It must have the following columns:

patient_id | gene | chr | pos | ref | alt | func | canc_type 
----|----|------|-----|-----|-----|------|-----

3. Run the following bash script to get the RNA VAFs for all SNVs within a single sample
```
bash get_RNA_VAF.sh -b ${rna_bam} -p ${patient_id} \
        -m ${mutation_table} -o ${out_dir} -w ${work_dir} \
        -s ${rna_sample_id}
```
