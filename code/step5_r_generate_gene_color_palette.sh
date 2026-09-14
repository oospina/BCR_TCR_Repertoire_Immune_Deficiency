#!/bin/bash
#SBATCH --job-name=colorpals
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=10:00:00

##
# Generate color palettes 
#

module load shared
module load R/4.4.2+Bioconductor

# Get gene usage files
igh=$( find ../data/dest/ | grep "_IGH_genes\.xlsx$" )
trb=$( find ../data/dest/ | grep "_TRB_genes\.xlsx$" )

# Make output folder
mkdir -p ../data/dest

Rscript --vanilla ./airrmine_MOD/airrnat/scripts/vgenecolors.R \
  "../data/dest/IGH_v_gene_colorcodes.tsv" \
  ${igh}

Rscript --vanilla ./airrmine_MOD/airrnat/scripts/vgenecolors.R \
  "../data/dest/TRB_v_gene_colorcodes.tsv" \
  ${trb}


