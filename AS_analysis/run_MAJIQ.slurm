#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

#Lines 12-21 change paths to correct directories, and lines 61-72, change SRR codes to correct codes for each sample group relevant to your dataset

#Define paths
SIF_PATH=/lyceum/ch5g20/majiq_latest.sif
AS_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis"
MAJIQ_DIR=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/MAJIQ
REF_DIR=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14
ANNOTATION=${REF_DIR}/Homo_sapiens.GRCh38.113.gff3
BUILD_DIR=${MAJIQ_DIR}/build
BAMS_DIR=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/bams
PSI_DIR=${MAJIQ_DIR}/psi
DELTA_DIR=${MAJIQ_DIR}/deltapsi
VOILA_DIR=${MAJIQ_DIR}/voila
CONFIG_FILE=/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/majiq_config.txt

#Make required directories
mkdir -p "$AS_DIR"
mkdir -p "$MAJIQ_DIR"
mkdir -p "$BUILD_DIR"
mkdir -p "$PSI_DIR"
mkdir -p "$DELTA_DIR"
mkdir -p "$VOILA_DIR"
mkdir -p logs

#Load required module
module load apptainer

#Run MAJIQ build
apptainer run --unsquash -H "$PWD" \
    --bind "$MAJIQ_DIR":"$MAJIQ_DIR" \
    --bind "$REF_DIR":"$REF_DIR" \
    --bind "$BAMS_DIR":"$BAMS_DIR" \
    "$SIF_PATH" \
    majiq build "$ANNOTATION" -c "$CONFIG_FILE" -o "$BUILD_DIR" -j 16 --min-intronic-cov 1 --simplify 

#Define variables to be input into MAJIQ psi
ALL_SAMPLES=($(ls "$BUILD_DIR"/*.majiq | xargs -n 1 basename | cut -d_ -f1))

#Run MAJIQ psi per sample
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    apptainer run --unsquash -H "$PWD" \
        --bind "$BUILD_DIR":"$BUILD_DIR" \
        --bind "$PSI_DIR":"$PSI_DIR" \
        "$SIF_PATH" \
        majiq psi "$BUILD_DIR/${SAMPLE}_Aligned.sortedByCoord.out.majiq" \
        -j 4 -o "$PSI_DIR" -n "$SAMPLE"
done

#Define full paths for deltapsi
CON_1=(
    $BUILD_DIR/SRR3306823_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306824_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306825_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306826_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306827_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306828_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306829_Aligned.sortedByCoord.out.majiq
)

CON_2=(
    $BUILD_DIR/SRR3306830_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306831_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306832_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306833_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306834_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306835_Aligned.sortedByCoord.out.majiq
    $BUILD_DIR/SRR3306836_Aligned.sortedByCoord.out.majiq
)

#Run MAJIQ deltapsi
apptainer run --unsquash -H "$PWD" \
    --bind "$BUILD_DIR":"$BUILD_DIR" \
    --bind "$DELTA_DIR":"$DELTA_DIR" \
    "$SIF_PATH" \
    majiq deltapsi \
    -grp1 "${CON_1[@]}" \
    -grp2 "${CON_2[@]}" \
    -j 1 -o "$DELTA_DIR" -n Condition_1 Condition_2

#Run MAJIQ voila
apptainer run --unsquash -H "$PWD" \
    --bind "$BUILD_DIR":"$BUILD_DIR" \
    --bind "$DELTA_DIR":"$DELTA_DIR" \
    --bind "$VOILA_DIR":"$VOILA_DIR" \
    "$SIF_PATH" \
    voila tsv -j 1 \
    --threshold 0.1 \
    --probability-threshold 0.95 \
    "$BUILD_DIR/splicegraph.sql" \
    "$DELTA_DIR/Condition_1-Condition_2.deltapsi.voila" \
    -f "$VOILA_DIR/Condition_1-Condition_2.deltapsi_padj_0.05_dPSI_0.1_probability-threshold.tsv"

apptainer run --unsquash -H "$PWD" \
    --bind "$BUILD_DIR":"$BUILD_DIR" \
    --bind "$DELTA_DIR":"$DELTA_DIR" \
    --bind "$VOILA_DIR":"$VOILA_DIR" \
    "$SIF_PATH" \
    voila tsv -j 1 \
    --threshold 0.05 \
    --probability-threshold 0.95 \
    "$BUILD_DIR/splicegraph.sql" \
    "$DELTA_DIR/Condition_1-Condition_2.deltapsi.voila" \
    -f "$VOILA_DIR/Condition_1-Condition_2.deltapsi_padj_0.05_dPSI_0.05_probability-threshold.tsv"

apptainer run --unsquash -H "$PWD" \
    --bind "$BUILD_DIR":"$BUILD_DIR" \
    --bind "$DELTA_DIR":"$DELTA_DIR" \
    --bind "$VOILA_DIR":"$VOILA_DIR" \
    "$SIF_PATH" \
    voila tsv -j 1 \
    --threshold 0.1 \
    --changing-pvalue-threshold 0.05 \
    "$BUILD_DIR/splicegraph.sql" \
    "$DELTA_DIR/Condition_1-Condition_2.deltapsi.voila" \
    -f "$VOILA_DIR/Condition_1-Condition_2.deltapsi_padj_0.05_dPSI_0.1_changing-pvalue-threshold.tsv"

apptainer run --unsquash -H "$PWD" \
    --bind "$BUILD_DIR":"$BUILD_DIR" \
    --bind "$DELTA_DIR":"$DELTA_DIR" \
    --bind "$VOILA_DIR":"$VOILA_DIR" \
    "$SIF_PATH" \
    voila tsv -j 1 \
    --threshold 0.05 \
    --changing-pvalue-threshold 0.05 \
    "$BUILD_DIR/splicegraph.sql" \
    "$DELTA_DIR/Condition_1-Condition_2.deltapsi.voila" \
    -f "$VOILA_DIR/Condition_1-Condition_2.deltapsi_padj_0.05_dPSI_0.05_changing-pvalue-threshold.tsv"
echo "MAJIQ pipeline complete."
