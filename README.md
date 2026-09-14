# Code Repository for Analysis of B/T Cell Receptor Sequencing From a Cohort of Patients With Immunodeficiency

This repository contains code associated with the analysis of B-cell and T-cell receptor sequencing data to replicate the results in the the manuscript "_Multi-omics approach to map a clinically diverse large cohort with homozygous founder RAG1 p.C176F variant among Old Order Mennonites_" ([Published in Cell](http://dx.doi.org/10.2139/ssrn.7096245)).

The code in this repository is largely based on the BCR/TCR sequencing and clone detection pipeline by Dr. Peter Blazso ([AIRRMINE](https://github.com/blazsop/airrmine)). 
Additional R code was contributed by Oscar Ospina and Marta Toth. Additional Unix code was contributed by Oscar Ospina. The Slurm job scheduler was used to analyze this data set on the Johns Hopkins University DISCOVERY HPC. Whenever possible and/or relevant, intermediate or result files are included in this repository. 

## TCR/BCR sequence processing
The raw .fastq files can be obtained from the Sequence Read Archive ([SRA accession PRJNA1526658](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1526658). Mapping and identifcation of CD3 sequences was performed with the [IMGT/HighV-QUEST tool](https://www.imgt.org/HighV-QUEST/). The pipeline generated compressed outputs (.txz files) that were used as input for the code in this repository.

## Modfied AIRRMINE pipeline (Dr. Peter Blazso in [Csomos et al. (2022), Nat Immunology](https://www.nature.com/articles/s41590-022-01271-6)) 
The modified AIRRMINE pipeline is located within the folder `code/airrmine_MOD`. Only minor changes were made to account for differences in our AIRR database formats and those from the manuscript originally presenting the pipeline. We thank Dr. Blazso for his help to re-implement the pipeline.

The AIRRMINE pipeline includes reference BCR sequences downloaded from the [IMGT database]((https://www.imgt.org/vquest/refseqh.html#refdir)). We also downloaded reference TCR sequences to be used by the pipeline.

## SLURM scripts to execute the modified AIRRMINE pipeline
A series of SLURM scripts with named `step*` are available within the `code` folder. Executed in order produces the clone identification outputs, including the clone treemaps.

## Manuscript figures/tables
The code in `code/BCRseq_Gene_usage_SHM_calculations.R` and `code/TCRseq_gene_usage_PCA.R` was used to generate plots for the quantitation of clone diversity present in the manuscript.