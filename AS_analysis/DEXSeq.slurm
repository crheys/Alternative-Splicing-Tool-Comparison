#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

#Define correct paths
AS_DIR=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis
DEXSEQ_DIR=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/DEXSeq

# Activate conda environment
source ~/.bashrc
conda activate DEXSeq_env

#Make necessary directories if they don't already exist
mkdir -p "$AS_DIR" "$DEXSEQ_DIR"

#Run DEXSeq in R
Rscript DEXSeq.R
