#!/bin/bash
#SBATCH --job-name=airr_net
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Generate network with germlines from AIRR
#

echo "PROCESS STARTED: $(date)"

module load shared
module load anaconda3

# Use environment containing phylip
source activate btcr_tools

module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/graph

fname="../data/clone/${id}_IGH_agsel.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/generate_nd.R \
  ${fname} \
  "../data/graph/${id}_IGH.gml"

fname="../data/clone/${id}_TRB_agsel.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/generate_nd.R \
  ${fname} \
  "../data/graph/${id}_TRB.gml"
  
echo "PROCESS FINISHED: $(date)"


