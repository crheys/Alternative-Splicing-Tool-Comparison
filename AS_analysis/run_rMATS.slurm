#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

# Make sure rMATS is installed 
# Change paths to correct directories: line 13-22, and change dataset specific info (read length and paired vs. single)
# Change sample names in the defining sample groups section to fit your dataset: line 44-

# Set paths - change them according to your working environments/ directories/ dataset info
CONDA_ENV_PATH="/mainfs/scratch/ch5g20/.conda/MISO_env"
BAM_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/bams"
GTF_FILE="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf"
AS_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/AS_analysis"
RMATS_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/AS_analysis/rMATS"
OUTPUT_DIR="$RMATS_DIR/output"
SIG_DIR_1="$OUTPUT_DIR/significant_padj_0.05_dPSI_0.1"
SIG_DIR_2="$OUTPUT_DIR/significant_padj_0.05_dPSI_0.05"
READ_LENGTH=101
TYPE=paired

# Activate environment
source ~/.bashrc
conda activate "$CONDA_ENV_PATH"

#Make required directories 
mkdir -p "$AS_DIR"
mkdir -p "$RMATS_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$SIG_DIR_1"
mkdir -p "$SIG_DIR_2"
mkdir -p "$SIG_DIR_3"

# Define sample groups
CON_1_SAMPLES=(
    "SRR22776675_Aligned.sortedByCoord.out.bam"
    "SRR22776676_Aligned.sortedByCoord.out.bam"
    "SRR22776677_Aligned.sortedByCoord.out.bam"
    "SRR22776686_Aligned.sortedByCoord.out.bam"
    "SRR22776687_Aligned.sortedByCoord.out.bam"
    "SRR22776688_Aligned.sortedByCoord.out.bam"
    "SRR22776689_Aligned.sortedByCoord.out.bam"
    "SRR22776690_Aligned.sortedByCoord.out.bam"
    "SRR22776691_Aligned.sortedByCoord.out.bam"
    "SRR22776692_Aligned.sortedByCoord.out.bam"
    "SRR22776693_Aligned.sortedByCoord.out.bam"
)

CON_2_SAMPLES=(
    "SRR22776660_Aligned.sortedByCoord.out.bam"
    "SRR22776678_Aligned.sortedByCoord.out.bam"
    "SRR22776679_Aligned.sortedByCoord.out.bam"
    "SRR22776680_Aligned.sortedByCoord.out.bam"
    "SRR22776681_Aligned.sortedByCoord.out.bam"
    "SRR22776682_Aligned.sortedByCoord.out.bam"
    "SRR22776683_Aligned.sortedByCoord.out.bam"
    "SRR22776684_Aligned.sortedByCoord.out.bam"
    "SRR22776685_Aligned.sortedByCoord.out.bam"
)

# Create BAM list files
CON1_FILE="Con1.txt"
CON2_FILE="Con2.txt"

# Join BAM paths with commas
(IFS=,; echo "${CON_1_SAMPLES[*]/#/$BAM_DIR/}") > "$CON1_FILE"
(IFS=,; echo "${CON_2_SAMPLES[*]/#/$BAM_DIR/}") > "$CON2_FILE"

# Run rMATS
rmats.py -t "$TYPE" \
  --b1 "$CON1_FILE" \
  --b2 "$CON2_FILE" \
  --gtf "$GTF_FILE" \
  --od "$OUTPUT_DIR" \
  --readLength "$READ_LENGTH"

#Filter output for padj < 0.05 and -0.1 > dPSI > 0.1
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' "$OUTPUT_DIR/SE.MATS.JCEC.txt" > "$SIG_DIR_1/SE.filtered.txt"
awk 'NR==1 || ($22 < 0.05 && ($25 > 0.1 || $25 < -0.1))' "$OUTPUT_DIR/MXE.MATS.JCEC.txt" > "$SIG_DIR_1/MXE.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' "$OUTPUT_DIR/RI.MATS.JCEC.txt" > "$SIG_DIR_1/RI.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' "$OUTPUT_DIR/A3SS.MATS.JCEC.txt" > "$SIG_DIR_1/A3SS.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.1 || $23 < -0.1))' "$OUTPUT_DIR/A5SS.MATS.JCEC.txt" > "$SIG_DIR_1/A5SS.filtered.txt"

#Filter output for padj	< 0.05 and -0.05 > dPSI > 0.05
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.05 || $23 < -0.05))' "$OUTPUT_DIR/SE.MATS.JCEC.txt" > "$SIG_DIR_2/SE.filtered.txt"
awk 'NR==1 || ($22 < 0.05 && ($25 > 0.05 || $25 < -0.05))' "$OUTPUT_DIR/MXE.MATS.JCEC.txt" > "$SIG_DIR_2/MXE.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.05 || $23 < -0.05))' "$OUTPUT_DIR/RI.MATS.JCEC.txt" > "$SIG_DIR_2/RI.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.05 || $23 < -0.05))' "$OUTPUT_DIR/A3SS.MATS.JCEC.txt" > "$SIG_DIR_2/A3SS.filtered.txt"
awk 'NR==1 || ($20 < 0.05 && ($23 > 0.05 || $23 < -0.05))' "$OUTPUT_DIR/A5SS.MATS.JCEC.txt" > "$SIG_DIR_2/A5SS.filtered.txt"
