#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

# Define base paths
BAM_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/bams"
SCRIPT_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/FTD-ALS_data/leafcutter/scripts"
CLUSTER_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/FTD-ALS_data/leafcutter/clustering"
AS_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis"
LEAFCUTTER_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/leafcutter"
JUNC_DIR="$LEAFCUTTER_DIR/bam_to_junc"
INTR_CLUST_DIR="$LEAFCUTTER_DIR/intron_clustering"
DIFF_SPLICING_DIR="$LEAFCUTTER_DIR/diff_splicing"
JUNC_LIST="$LEAFCUTTER_DIR/intron_clustering/junc_file.txt"
EXONS="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/leafcutter_exons.txt"
GROUP_FILE="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/leafcutter/groups.txt"


#Needed for differential splicing input, tells leafcutter the sample number of the group with the least samples
SAMPLE_GROUP_NUM=7

# Create necessary output directories
mkdir -p "$AS_DIR"
mkdir -p "$LEAFCUTTER_DIR"
mkdir -p "$JUNC_DIR"
mkdir -p "$INTR_CLUST_DIR"
mkdir -p "$DIFF_SPLICING_DIR"
mkdir -p tmp

# List of sample IDs
SAMPLES=(
    SRR3306823
    SRR3306824
    SRR3306825
    SRR3306826
    SRR3306827
    SRR3306828
    SRR3306829
    SRR3306830
    SRR3306831
    SRR3306832
    SRR3306833
    SRR3306834
    SRR3306835
    SRR3306836
)

# Empty the junc list file
> "$JUNC_LIST"

# Activate conda environment
source ~/.bashrc
conda activate /mainfs/scratch/ch5g20/.conda/MISO_env

# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."

    TMP_DIR="$JUNC_DIR/${SAMPLE}"
    mkdir -p "$TMP_DIR"

    BAM_FILE="$BAM_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    BED_FILE="$TMP_DIR/${SAMPLE}.bed"
    JUNC_FILE="$JUNC_DIR/${SAMPLE}.junc"

    # Convert BAM to BED
    samtools view "$BAM_FILE" | \
    "$SCRIPT_DIR/filter_cs.py" | \
    "$SCRIPT_DIR/sam2bed.pl" --use-RNA-strand - "$BED_FILE"

    # Convert BED to JUNC
    "$SCRIPT_DIR/bed2junc.pl" "$BED_FILE" "$JUNC_FILE"

    # Append to junc list
    echo "$JUNC_FILE" >> "$JUNC_LIST"

    # Remove temporary directory
    rm -r "$TMP_DIR"
done

# Activate conda environment
source ~/.bashrc
conda activate /mainfs/scratch/ch5g20/leafcutter_env

# Run LeafCutter clustering
python "$CLUSTER_DIR/leafcutter_cluster.py" \
-j "$JUNC_LIST" \
-r "$INTR_CLUST_DIR"

echo "LeafCutter clustering complete."



#Run python script to generate groups.txt file based on format of intron clustering output file, needed for input in diff splicing
leafcutter_generate_groups.py

#Run differential splicing analysis with leafcutter
echo "Running differential splicing analysis with LeafCutter..."

Rscript "$SCRIPT_DIR/leafcutter_ds.R" \
    "$INTR_CLUST_DIR/leafcutter_perind_numers.counts.gz" \
    "$GROUP_FILE" \
    -e "$EXONS" \
    --output_prefix "$DIFF_SPLICING_DIR/leafcutter_results" \
    -i "$SAMPLE_GROUP_NUM"

# -------------------------
# Filtering results
# -------------------------
echo "Filtering significant clusters..."

cd "$DIFF_SPLICING_DIR"
mkdir -p significant 

# 1. Filter for adjusted p-value (FDR < 0.05)
awk -F'\t' 'NR == 1 || $5 < 0.05' leafcutter_results_cluster_significance.txt > significant/significant_results.txt

# 2. Filter for absolute effect size (delta PSI > 0.1 or < -0.1)
awk -F'\t' 'NR == 1 || $5 > 0.1 || $5 < -0.1' leafcutter_results_effect_sizes.txt > significant/significant_effect_sizes.txt

cd significant 

# 3. Combine both files to get final significant events
awk 'NR==FNR && FNR > 1 {genes[$1]=$8; next}  
     FNR == 1 {print $0 "\tgene"; next}  
     $1 in genes {print $0 "\t" genes[$1]}' \
     significant_results.txt \
     significant_effect_sizes.txt \
     > final_significant_events.txt

echo "Differential splicing and filtering complete."
