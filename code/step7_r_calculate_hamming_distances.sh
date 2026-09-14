#!/bin/bash
#SBATCH --job-name=hamming_dist
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00
#SBATCH -a 1-36%36

##
# Calculate Hamming distances and define clone split thresholds
#

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/dest
mkdir -p ../results/hamming_dist

fname="../data/coll/${id}_IGH_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/findthreshold.R \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../data/dest/${id}_IGH_dist_thrs.data" \
  "../results/hamming_dist/${id}_IGH_threshold_histo.pdf"

fname="../data/coll/${id}_TRB_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/findthreshold.R \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../data/dest/${id}_TRB_dist_thrs.data" \
  "../results/hamming_dist/${id}_TRB_threshold_histo.pdf"


