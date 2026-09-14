#!/bin/bash
#SBATCH --job-name=viz_airr_net
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Plot AIRR germline network
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../results/airr_network

fname="../data/graph/${id}_IGH.gml"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/visualize_nd.R \
  ${fname} \
  "../results/airr_network/${id}_IGH_AIRRND.pdf" \
  "../results/airr_network/${id}_IGH_AIRRND_agsel.pdf" \
  "../results/airr_network/${id}_IGH_AIRRND_clone.pdf" \
  "../results/airr_network/${id}_IGH_AIRRND_mutations.pdf" \
  "../data/dest/IGH_v_gene_colorcodes.tsv"

fname="../data/graph/${id}_TRB.gml"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/visualize_nd.R \
  ${fname} \
  "../results/airr_network/${id}_TRB_AIRRND.pdf" \
  "../results/airr_network/${id}_TRB_AIRRND_agsel.pdf" \
  "../results/airr_network/${id}_TRB_AIRRND_clone.pdf" \
  "../results/airr_network/${id}_TRB_AIRRND_mutations.pdf" \
  "../data/dest/TRB_v_gene_colorcodes.tsv"
  
echo "PROCESS FINISHED: $(date)"


