#!/bin/bash
#SBATCH --job-name=count_mutation
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
mkdir -p ../results/clone_stats

fname="../data/clone/${id}_IGH_withgl.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ../code/airrmine_MOD/airrnat/scripts/R_S_mutations.R \
  ${fname} \
  "../results/clone_stats/${id}_IGH_rsmut.xlsx"

fname="../data/clone/${id}_TRB_withgl.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ../code/airrmine_MOD/airrnat/scripts/R_S_mutations.R \
  ${fname} \
  "../results/clone_stats/${id}_TRB_rsmut.xlsx"
  
echo "PROCESS FINISHED: $(date)"


