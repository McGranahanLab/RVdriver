# RVdriver
A statistical approach to detect cancer driver genes utilising the RNA Variant Allele Frequency of somatic mutations

## Inputs

RVdriver only requires a simple tab deliminted mutation table in the following input

patient_id | gene | CSRA | RNA_ref_count | RNA_alt_count | RNA_VAF | func | canc_type 
----|----|------|-----|------|------|------|------

an example mutation table is geiven [here](./assets/UVM_example_data.txt)
