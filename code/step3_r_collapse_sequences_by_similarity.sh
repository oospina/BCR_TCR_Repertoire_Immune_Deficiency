#!/bin/bash
#SBATCH --job-name=collapse_seqs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --time=10:00:00
#SBATCH -a 1-36%36

##
# Collapse sequences based on similarity and copy abundance
#

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/coll

fname="../data/prep/${id}_IGH_prepared.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/collapserep.R \
  ${fname} \
  "../data/coll/${id}_IGH_collapsed.tsv"

fname="../data/prep/${id}_TRB_prepared.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/collapserep.R \
  ${fname} \
  "../data/coll/${id}_TRB_collapsed.tsv"


