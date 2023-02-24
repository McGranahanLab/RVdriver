# RVdriver
A statistical approach to detect cancer driver genes utilising the RNA Variant Allele Frequency of somatic mutations

## Inputs

RVdriver only requires a simple tab deliminted mutation table in the following input

patient_id | gene | CSRA | RNA_ref_count | RNA_alt_count | RNA_VAF | func | canc_type 
----|----|------|-----|------|------|------|------

an example mutation table is geiven [here](./assets/UVM_example_data.txt)

### The following steps are needed to run RVdriver on the example dataset:
1. clone this repo
2. if necessary, create the required conda environment (all that is required to run RVdriver is argparse, tidyverse and lmerTest)
```
conda env create --name RVdriver --file=RVdriver.yml
```
4. Activate the conda environment
```
conda activate RVdriver
```
5. run RVdriver
```
Rscript run_RVdriver.R --mut_path assets/UVM_example_data.txt --canc_type UVM --out_dir test
```
The results should be the same as those presented [here](./test_data_results/UVM_rvdriver_results.csv)
