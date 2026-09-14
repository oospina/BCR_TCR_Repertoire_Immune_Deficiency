#!/bin/bash
#SBATCH --job-name=def_germlines
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Define germlines
#

echo "PROCESS STARTED: $(date)"

module load shared
module load anaconda3

# Use environment containing Change-O adn load R
source activate btcr_tools

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/clone/

fname="../data/clone/${id}_IGH_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
CreateGermlines.py \
  -d "${fname}" \
  -g dmask \
  --cloned \
  -r ./airrmine_MOD/samples/igdata/imgt/human/IGHV.fasta \
     ./airrmine_MOD/samples/igdata/imgt/human/IGHD.fasta \
     ./airrmine_MOD/samples/igdata/imgt/human/IGHJ.fasta \
  -o "../data/clone/${id}_IGH_withgl.tsv"

fname="../data/clone/${id}_TRB_allclone.tsv"
echo "PROCESSING FILE: ${fname}"
CreateGermlines.py \
  -d "${fname}" \
  -g dmask \
  --cloned \
  -r ../data/IMGT_refs/trdata/imgt/human/TRBV.fasta \
     ../data/IMGT_refs/trdata/imgt/human/TRBD.fasta \
     ../data/IMGT_refs/trdata/imgt/human/TRBJ.fasta \
  -o "../data/clone/${id}_TRB_withgl.tsv"
  
echo "PROCESS FINISHED: $(date)"


