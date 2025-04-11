# RVdriver
A statistical approach to detect cancer driver genes utilising the RNA Variant Allele Frequency of somatic mutations

## Inputs

RVdriver only requires a simple tab deliminted or .RDS file mutation table in the following format

patient_id | gene | CSRA | RNA_ref_count | RNA_alt_count | RNA_VAF | func | canc_type 
----|----|------|-----|------|------|------|------

an example mutation table is given [here](./inst/example_data/UCS_mutation_table.rds)

within the ```func``` column, the following mutation types are accepted:
```
Missense_Mutation
Nonsense_Mutation
Nonstop_Mutation
Splice_Site
Silent
```

### The following steps are needed to run RVdriver on the example dataset:
1. Install RVdriver and dependencies
```
devtools::install_github("https://github.com/McGranahanLab/RVdriver.git")
```
2. Run the [example script](./inst/scripts/run_RVdriver_example.R) for either UCS or UVM
This will create a results table wherever specified which should match up the results found [here for UCS](.inst/example_data/UCS_results) or [here for UVM](.inst/example_data/UVM_results)

The same script can be used for your own paired genomic-transcriptomic mutation data!

#### The following steps were taken to calculate the RNA VAFs for somatic mutations across the pan cancer TCGA dataset

1. Download the following files and store them in the assets directory. These files are required for GATK preprocessing (MarkDuplicates and BQSR):
    - [TCGA reference genome](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files) 
    - [dbsnp_138_vcf](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf)
    - [dbsnp_138_vcf_idx](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx)

2. Pull the relevent singularity container 
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
