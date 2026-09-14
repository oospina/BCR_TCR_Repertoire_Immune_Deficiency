#!/bin/bash
#SBATCH --job-name=gene_usage
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00
#SBATCH -a 1-36%36

##
# Calculate gene (V,D,J) usages
#

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/dest

fname="../data/coll/${id}_IGH_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/geneusage.R \
  ${fname} \
  "../data/dest/${id}_IGH_genes.xlsx"
  
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/geneusage_m.R \
  ${fname} \
  "../data/dest/${id}_IGH_genes.tsv"

fname="../data/coll/${id}_TRB_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/geneusage.R \
  ${fname} \
  "../data/dest/${id}_TRB_genes.xlsx"
  
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/geneusage_m.R \
  ${fname} \
  "../data/dest/${id}_TRB_genes.tsv"


