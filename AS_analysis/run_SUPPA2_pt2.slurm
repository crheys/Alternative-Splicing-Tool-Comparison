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
SUPPA_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/SUPPA2"
EVENT_OUTPUT_DIR="$SUPPA_DIR/Events"
QUANT_DIR="$SUPPA_DIR/Quant"
QUANT_OUT="$QUANT_DIR"
TPM_MATRIX="$SUPPA_DIR/Quant/iso_tpm.txt"
PSI_DIR="$SUPPA_DIR/PSI"
DIFFSPLICE_DIR="$SUPPA_DIR/DiffSplice"
SIGNIFICANT_DIR_1="$SUPPA_DIR/DiffSplice/Significant_dPSI_0.1"
SIGNIFICANT_DIR_2="$SUPPA_DIR/DiffSplice/Significant_dPSI_0.05"

# Activate environment
source ~/.bashrc
conda activate "$SUPPA2_ENV_PATH"

mkdir -p "$DIFFSPLICE_DIR" "$SIGNIFICANT_DIR"

# Group 1 and 2 TPM matrices - change these values based on the arrangment of the iso_tpm.txt file (column 1 must always be added as its the gene ids, then add the columns for samples in each group to the two condition files)
cut -f 1,2,3,4,5,6,7,8 "$TPM_MATRIX" > "$PSI_DIR/con1.tmp"
cut -f 1,9,10,11,12,13,14,15 "$TPM_MATRIX" > "$PSI_DIR/con2.tmp"

# Group-specific PSI files for each event type
for EVENT in A3 A5 AF AL MX RI SE; do
  cut -f 1,2,3,4,5,6,7,8 "$PSI_DIR/${EVENT}_events.psi" > "$PSI_DIR/${EVENT}_con1.psi"
  cut -f 1,9,10,11,12,13,14,15 "$PSI_DIR/${EVENT}_events.psi" > "$PSI_DIR/${EVENT}_con2.psi"
done

# STEP 6: Run diffSplice for each event type
for EVENT in A3 A5 AF AL MX RI SE; do
  python "$SUPPA_PATH/suppa.py" diffSplice \
    -m empirical -gc \
    -i "$EVENT_OUTPUT_DIR/GRCh38_${EVENT}_strict.ioe" \
    --save_tpm_events \
    -p "$PSI_DIR/${EVENT}_con1.psi" "$PSI_DIR/${EVENT}_con2.psi" \
    -e "$PSI_DIR/con1.tmp" "$PSI_DIR/con2.tmp" \
    -o "$DIFFSPLICE_DIR/diffsplice_${EVENT}"
done

# STEP 6: Filtering significant events
for EVENT in A3 A5 AF AL MX RI SE; do
  echo "Filtering significant events for $EVENT..."
  awk 'NR==1 || ($3 <= 0.05 && ($2 >= 0.1 || $2 <= -0.1))' "$DIFFSPLICE_DIR/diffsplice_${EVENT}.dpsi" > "$SIGNIFICANT_DIR_1/significant_${EVENT}_events.dpsi"
done

for EVENT in A3 A5 AF AL MX RI SE; do
  echo "Filtering significant events for $EVENT..."
  awk 'NR==1 || ($3 <= 0.05 && ($2 >= 0.05 || $2 <= -0.05))' "$DIFFSPLICE_DIR/diffsplice_${EVENT}.dpsi" > "$SIGNIFICANT_DIR_2/significant_${EVENT}_events.dpsi"
done
