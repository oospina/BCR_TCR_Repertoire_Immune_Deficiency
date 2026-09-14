#!/bin/bash
#SBATCH --job-name=vgene_clone_tmap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00
#SBATCH -a 1-36%36

##
# Generate V-gene clone abundances treemap
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../results/vgene_abund_treemap

fname="../data/clone/${id}_IGH_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/treemap_clonal.R \
  "../data/dest/IGH_v_gene_colorcodes.tsv" \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/${id}_IGH_clone_tmap.pdf"

fname="../data/clone/${id}_TRB_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/treemap_clonal.R \
  "../data/dest/TRB_v_gene_colorcodes.tsv" \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/${id}_TRB_clone_tmap.pdf"
  
echo "PROCESS FINISHED: $(date)"


