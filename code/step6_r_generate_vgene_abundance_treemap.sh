#!/bin/bash
#SBATCH --job-name=vgene_abund_tmap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00
#SBATCH -a 1-36%36

##
# Generate V gene abundance treemap
#

module load shared
module load R/4.4.2+Bioconductor

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../results/vgene_abund_treemap

fname="../data/coll/${id}_IGH_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/treemap.R \
  "../data/dest/IGH_v_gene_colorcodes.tsv" \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/${id}_IGH_tmap.pdf"

fname="../data/coll/${id}_TRB_collapsed.tsv"
Rscript --vanilla ./airrmine_MOD/airrnat/scripts/treemap.R \
  "../data/dest/TRB_v_gene_colorcodes.tsv" \
  ${fname} \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/${id}_TRB_tmap.pdf"
  

