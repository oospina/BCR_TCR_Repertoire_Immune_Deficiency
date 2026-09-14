#!/bin/bash
#SBATCH --job-name=define_clones
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=1:00:00
#SBATCH -a 1-36%36

##
# Define clones using Change-O and calculated Hamming distance split threshold
#

echo "PROCESS STARTED: $(date)"

module load shared
module load anaconda3

# Use environment containing Change-O adn load R
source activate btcr_tools

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/clone

fname="../data/coll/${id}_IGH_collapsed.tsv"
hamd=$( cat "../data/dest/${id}_IGH_dist_thrs_MOD.data" )
echo "PROCESSING FILE: ${fname}"
DefineClones.py \
  -d ${fname} \
  --nproc 32 \
  --mode gene \
  --act set \
  --norm len \
  --model ham \
  --dist ${hamd} \
  -o "../data/clone/${id}_IGH_allclone.tsv"

fname="../data/coll/${id}_TRB_collapsed.tsv"
hamd=$( cat "../data/dest/${id}_TRB_dist_thrs_MOD.data" )
echo "PROCESSING FILE: ${fname}"
DefineClones.py \
  -d ${fname} \
  --nproc 32 \
  --mode gene \
  --act set \
  --norm len \
  --model ham \
  --dist ${hamd} \
  -o "../data/clone/${id}_TRB_allclone.tsv"

echo "PROCESS FINISHED: $(date)"


