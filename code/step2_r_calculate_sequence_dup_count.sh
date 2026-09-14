#!/bin/bash
#SBATCH --job-name=calc_dups
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Calculate duplicate_count column within AIRR file
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/prep

fname="../data/tsv/${id}_IGH_db-pass.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./utils/calculate_dup_counts.R \
  ${fname} \
  "../data/prep/${id}_IGH_prepared.tsv"

fname="../data/tsv/${id}_TRB_db-pass.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./utils/calculate_dup_counts.R \
  ${fname} \
  "../data/prep/${id}_TRB_prepared.tsv"

echo "PROCESS FINISHED: $(date)"


