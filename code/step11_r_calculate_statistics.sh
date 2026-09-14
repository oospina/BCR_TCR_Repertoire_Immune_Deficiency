#!/bin/bash
#SBATCH --job-name=diversity_stats
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Calculate diversity statistic
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../results/clone_stats

fname="../data/clone/${id}_IGH_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/cstats.R \
  ${fname} \
  "../results/clone_stats/${id}_IGH_cstats.tsv"

fname="../data/clone/${id}_TRB_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/cstats.R \
  ${fname} \
  "../results/clone_stats/${id}_TRB_cstats.tsv"
  
echo "PROCESS FINISHED: $(date)"


