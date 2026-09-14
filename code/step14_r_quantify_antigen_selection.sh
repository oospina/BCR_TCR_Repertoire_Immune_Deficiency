#!/bin/bash
#SBATCH --job-name=select_baseline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Calculate mutation frequencies
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/clone

fname="../data/clone/${id}_IGH_withgl.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/agselection.R \
  ${fname} \
  "../data/clone/${id}_IGH_agsel.tsv"

fname="../data/clone/${id}_TRB_withgl.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/agselection.R \
  ${fname} \
  "../data/clone/${id}_TRB_agsel.tsv"
  
echo "PROCESS FINISHED: $(date)"


