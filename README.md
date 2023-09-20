# Description

The repository contains the data analysis pipeline for UKBiobank data, focusing on the analysis of heterozygous carriers of Pathogenic/Likely Pathogenic variaints in recessive genes. 

The pipeline scripts are under `scripts`, reusable code is under `ukbb_recessive`. 

There are 2 main pipelines: 

1. Collection of singleton high qualuty LoF variants in highly constrained (s-het >=0.15) genes from the UK Biobank. 

2. Collection of PLPs in recessive genes from the UK Biobank. 

Data from both of the pipelines are analyzed using generalized linear models. 

# Pipelines

The reason for these pipelines to be difefrent is that historically recessive genes pipeline was developed first. Therefore in order to add information from other genes, the new approach was developed independently. 

## Collection of singleton high qualuty LoF variants in highly constrained (s-het >=0.15) genes from the UK Biobank.

The pipeline is under `scripts/lof_carriers_selection` and consists of the following steps:

1. Creation of the list of highly constrained genes with genomic locations (`1_create_gene_list.ipynb`)

2. Variant collection procedure, executed on RAP (`2_data_collection_paper_rap_450k.ipynb`)

3. Variant selection procedure -- singleton high-quality LoF (`3_variants_vep_annotation.ipynb`)

## Collection of PLPs in recessive genes from the UK Biobank.

The pipeline is under `scripts/pipeline` and consists of the following steps:

0. Converting all the conserning genes into gencode v.34 notation (`0_convert_data_into_gencode_v34.ipynb`)

1. Creation of the list of all unrelated samples (`1_generate_unrelaed_samples.ipynb`)

2. Variant collection procedure for 1929 recessive genes executed in RAP (`2_collect_variants_rap.ipynb`)

3. Variant and genotype normalization procedure from multi-sample VCF format to single-sample VCF format (`3_normalize_data.ipynb`)

4. PLP selection procedure (`5_plps_selection.ipynb`)

5. Filtration of the cohort data to contain PLPs-only variants (`6_filter_plps_in_cohort_data.ipynb`)

# Regressions

Regressions analysis combines previous two pipelines outputs and therefore described separatly. The pipeline is under `scripts/pipeline/7_regressions`. 