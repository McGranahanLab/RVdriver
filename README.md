# RVdriver
A statistical approach to detect cancer driver genes utilising the RNA Variant Allele Frequency of somatic mutations

## Inputs

RVdriver only requires a simple tab deliminted mutation table in the following format

patient_id | gene | CSRA | RNA_ref_count | RNA_alt_count | RNA_VAF | func | canc_type 
----|----|------|-----|------|------|------|------

an example mutation table is given [here](./assets/UVM_example_data.txt)

### The following steps are needed to run RVdriver on the example dataset:
1. clone this repo
if necessary, create the required conda environment (all that is required to run RVdriver is argparse, tidyverse and lmerTest)
```
conda env create --name RVdriver --file=RVdriver.yml
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
```
4. run RVdriver
```
Rscript run_RVdriver.R --mut_path assets/UVM_example_data.txt --canc_type UVM --out_dir test
```
The results should be the same as those presented [here](./test_data_results/UVM_rvdriver_results.csv)

#### The following steps were taken to get the RNA VAFs for somatic mutations across the pan cancer TCGA dataset

1. Download the following files and store them in the assets directory. These files are required for GATK preprocessing (MarkDuplicates and BQSR):
    - [TCGA reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) 
    - [dbsnp_138_vcf](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf)
    - [dbsnp_138_vcf_idx](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx)

2. Pull the RVdriver singularity container
```
mkdir singularity_images
cd singularity_images/
singularity pull --arch amd64 library://tpjones15/default/rvdriver:latest
```    
The mutation table can be a single table with all samples or multiple sample specific tables. It must have the following columns:

patient_id | gene | chr | pos | ref | alt | func | canc_type 
----|----|------|-----|-----|-----|------|-----

3. Run the following bash script to get the RNA VAF for a single sample
```
bash get_RNA_VAF.sh -b ${rna_bam_path} -p ${patient_id} \
        -m ${mut_path} -o ${out_dir} -w ${work_dir} \
        -s ${rna_sample_id}
```
