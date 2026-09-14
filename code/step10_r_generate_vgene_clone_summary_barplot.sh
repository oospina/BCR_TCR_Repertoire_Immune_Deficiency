#!/bin/bash
#SBATCH --job-name=vgene_clone_bplot
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00

##
# Generate V-gene clone summary bar plot
#

echo "PROCESS STARTED: $(date)"

module load shared
module load R/4.4.2+Bioconductor

# Make output folder
mkdir -p ../results/vgene_abund_treemap

fname=$( find ../data/clone/ | grep "_IGH_allclone.tsv" )
Rscript --vanilla ./utils/barplot_clonal.R \
  "../data/dest/IGH_v_gene_colorcodes.tsv" \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/IGH_clone_summary_bplot.pdf" \
  ${fname}

fname=$( find ../data/clone/ | grep "_TRB_allclone.tsv" )
Rscript --vanilla ./utils/barplot_clonal.R \
  "../data/dest/TRB_v_gene_colorcodes.tsv" \
  ../data/raw_data/tcrbcr_index_sample_key.txt \
  "../results/vgene_abund_treemap/TRB_clone_summary_bplot.pdf" \
  ${fname}
  
echo "PROCESS FINISHED: $(date)"


