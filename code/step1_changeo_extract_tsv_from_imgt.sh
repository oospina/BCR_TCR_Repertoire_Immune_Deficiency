#!/bin/bash
#SBATCH --job-name=extract_tsv
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=256GB
#SBATCH --time=2:00:00
#SBATCH -a 1-36%36

##
# Extract TSV database from IMGT output
#

echo "PROCESS STARTED: $(date)"

module load shared
module load anaconda3

# Use environment containing Change-O and load R
source activate btcr_tools

# Get task ID
id=$( head -n ${SLURM_ARRAY_TASK_ID} ../data/raw_data/tcrbcr_index_sample_key.txt | tail -n1 | cut -f1 )

# Make output folder
mkdir -p ../data/tsv

fname="../data/raw_data/imgt/${id}_IGH.txz"
echo "PROCESSING FILE: ${fname}"
MakeDb.py imgt \
  -i ${fname} \
  --extended \
  --outdir ../data/tsv/ \
  --outname "${id}_IGH"

fname="../data/raw_data/imgt/${id}_TRB.txz"
echo "PROCESSING FILE: ${fname}"
MakeDb.py imgt \
  -i ${fname} \
  --extended \
  --outdir ../data/tsv/ \
  --outname "${id}_TRB"

echo "PROCESS FINISHED: $(date)"


