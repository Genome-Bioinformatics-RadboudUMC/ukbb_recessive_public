# Description

The repository contains the data analysis pipeline for UKBiobank data, focusing on the analysis of heterozygous carriers of Pathogenic/Likely Pathogenic variaints in recessive genes. 

The pipeline scripts are under `scripts`, reusable code is under `ukbb_recessive`. 

There are 2 different variant selection procedures: 

1. Collection of singleton high qualuty LoF variants in highly constrained (s-het >=0.15) genes from the UK Biobank. 

2. Collection of PLPs in recessive genes from the UK Biobank. 

Data from both of the sources are analyzed using generalized linear models. 

# Pipeline

## Collection of singleton high qualuty LoF variants in highly constrained (s-het >=0.15) genes from the UK Biobank.

The pipeline is under `scripts/1_data_preparation/lof_carriers_selection` and consists of the following steps:

1. Creation of the list of highly constrained genes with genomic locations (`1_create_gene_list.ipynb`)

2. Variant collection procedure, executed on RAP (`2_data_collection_paper_rap_450k.ipynb`)

3. Variant selection procedure -- singleton high-quality LoF (`3_variants_vep_annotation.ipynb`)

## Collection of PLPs in recessive genes from the UK Biobank.

The pipeline is under `scripts/pipeline/1_data_preparation` and consists of the following steps:

0. Converting all the conserning genes into gencode v.34 notation (`0_convert_data_into_gencode_v34.ipynb`)

1. Creation of the list of all unrelated samples (`1_generate_unrelaed_samples.ipynb`)

2. Variant collection procedure for 1929 recessive genes executed in RAP (`2_collect_variants_rap.ipynb`)

3. Variant and genotype normalization procedure from multi-sample VCF format to single-sample VCF format (`3_normalize_data.ipynb`)

4. PLP selection procedure (`5_plps_selection.ipynb`)

5. Filtration of the cohort data to contain PLPs-only variants (`6_filter_plps_in_cohort_data.ipynb`)

# Regressions

Regressions analysis combines previous two pipelines outputs and therefore described separatly. The pipeline is under `scripts/pipeline/2_regressions`. It consusts of the different steps and analyses:

0. Creation of the dataset for all of consequent regression analyses (`0_dataset_creation`)

1. Regressions for all phenotypes for PLPs in recessive genes and singleton LoFs in non-recessive highly constrained genes (`1_basic_regressions`)

2. Regressions for different disorder groups (`2_disorder_groups`)

3. Analyses based on sampling  (`3_sampling_analysis`):
    a. Sampling dataset size (`dataset`)
    b. Samling non-ID recessive genes to resemble the same s-het distribution as ID genes (`genes`)
    c. Sampling synonymous variants as a negative control (`synonymous`)

4. Sex-specific analyses (`4_sex_specific_analyses`)

5. Consanguinity ratio (`5_cr_analysis`)

6. Correction for deprivation scores (`6_deprivation_abalysis`)

# Figures

All main and supplementary/extended figures generation is under `scripts/pipeline/3_generate_figures`.
