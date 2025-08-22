#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

# CONFIGURATION
SUPPA2_ENV_PATH="/mainfs/scratch/ch5g20/SUPPA2_env"
SALMON_ENV_PATH="/mainfs/scratch/ch5g20/.conda/salmon_env"
SUPPA_PATH="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/plant/SUPPA2/SUPPA_v2.2/SUPPA-2.2"
GTF_FILE="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf"
TRANSCRIPT_INDEX="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/FTD/AS_analysis/SUPPA2/transcript_index"
AS_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis"
SUPPA_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/SUPPA2"
EVENT_OUTPUT_DIR="$SUPPA_DIR/Events"
QUANT_DIR="$SUPPA_DIR/Quant"
QUANT_OUT="$QUANT_DIR"
TPM_MATRIX="$SUPPA_DIR/Quant/iso_tpm.txt"
PSI_DIR="$SUPPA_DIR/PSI"
SIGNIFICANT_DIR="$SUPPA_DIR/PSI/Significant"

# Sample FASTQ pairs (Update these based on your dataset)
declare -A SAMPLES=(
  ["SRR3306823"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306823_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306823_R2_paired.fastq.gz"
  ["SRR3306824"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306824_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306824_R2_paired.fastq.gz"
  ["SRR3306825"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306825_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306825_R2_paired.fastq.gz"
  ["SRR3306826"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306826_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306826_R2_paired.fastq.gz"
  ["SRR3306827"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306827_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306827_R2_paired.fastq.gz"
  ["SRR3306828"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306828_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306828_R2_paired.fastq.gz"
  ["SRR3306829"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306829_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306829_R2_paired.fastq.gz"
  ["SRR3306830"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306830_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306830_R2_paired.fastq.gz"
  ["SRR3306831"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306831_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306831_R2_paired.fastq.gz"
  ["SRR3306832"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306832_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306832_R2_paired.fastq.gz"
  ["SRR3306833"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306833_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306833_R2_paired.fastq.gz"
  ["SRR3306834"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306834_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306834_R2_paired.fastq.gz"
  ["SRR3306835"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306835_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306835_R2_paired.fastq.gz"
  ["SRR3306836"]="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306836_R1_paired.fastq.gz /mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/filtering/SRR3306836_R2_paired.fastq.gz"
  # Add more samples here
)

# Assign groupings based on metadata (make sure it's even)
GROUP1=("SRR3306823" "SRR3306824" "SRR3306825" "SRR3306826" "SRR3306827" "SRR3306828" "SRR3306829")
GROUP2=("SRR3306830" "SRR3306831" "SRR3306832" "SRR3306833" "SRR3306834" "SRR3306835" "SRR3306836")

# Activate environment
source ~/.bashrc
conda activate "$SUPPA2_ENV_PATH"

# Create output dirs
mkdir -p "$AS_DIR" "$SUPPA_DIR" "$EVENT_OUTPUT_DIR" "$QUANT_DIR" "$PSI_DIR" "$SIGNIFICANT_DIR"

# STEP 1: Generate splicing events
python "$SUPPA_PATH/suppa.py" generateEvents \
  -i "$GTF_FILE" \
  -o "$EVENT_OUTPUT_DIR/GRCh38" \
  -e SE SS MX RI FL

#Deactivate SUPPA2 environment and activate salmon environment
conda deactivate
conda activate "$SALMON_ENV_PATH"

# STEP 2: Quantification with Salmon
for SAMPLE in "${!SAMPLES[@]}"; do
  read -r R1 R2 <<< "${SAMPLES[$SAMPLE]}"
  salmon quant -i "$TRANSCRIPT_INDEX" -p 6 --libType A \
    -1 "$R1" -2 "$R2" \
    -o "$QUANT_OUT/${SAMPLE}.quant"
done

#Deactivate salmon environment and activate SUPPA2 environment 
conda deactivate
conda activate "$SUPPA2_ENV_PATH"

# STEP 3: Generate TPM matrix
QUANT_FILES=""
for SAMPLE in "${!SAMPLES[@]}"; do
  QUANT_FILES+="$QUANT_OUT/${SAMPLE}.quant/quant.sf "
done

python "$SUPPA_PATH/multipleFieldSelection.py" -i $QUANT_FILES -k 1 -f 4 -o "$TPM_MATRIX"

# STEP 4: Calculate PSI per event
EVENT_TYPES=(SE RI MX AL AF A3 A5)
for EVENT in "${EVENT_TYPES[@]}"; do
  python "$SUPPA_PATH/suppa.py" psiPerEvent \
    -i "$EVENT_OUTPUT_DIR/GRCh38_${EVENT}_strict.ioe" \
    -e "$TPM_MATRIX" \
    -o "$PSI_DIR/${EVENT}_events"
done
