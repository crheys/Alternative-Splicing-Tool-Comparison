#!/bin/bash
#SBATCH --job-name=intron_exon_batch
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=%j-intron_exon_batch.out
#SBATCH --error=%j-intron_exon_batch.err

# === Load environment ===
source ~/.bashrc
conda activate MISO_env

# === Paths ===
BAM_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/bams"
INTRON_BED="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/FTD/bams/rseqc/introns_trimmed.bed"
EXON_BED="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/FTD/bams/rseqc/exons_trimmed.bed"
OUTFILE="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/bams/intron_exon_summary.tsv"

# === Check input ===
if [[ ! -f "$INTRON_BED" || ! -f "$EXON_BED" ]]; then
    echo "ERROR: introns.bed or exons.bed file missing!"
    exit 1
fi

# === Header for summary file ===
echo -e "Sample\tIntronic_Reads\tExonic_Reads\tTotal_Reads\tFraction_Intronic" > "$OUTFILE"

# === Loop through BAM files ===
for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .bam)

    echo "Processing $SAMPLE..."

    INTRON_READS=$(bedtools coverage -a "$INTRON_BED" -b "$BAM" | awk '{sum+=$7} END {print sum+0}')
    EXON_READS=$(bedtools coverage -a "$EXON_BED" -b "$BAM" | awk '{sum+=$7} END {print sum+0}')
    TOTAL_READS=$((INTRON_READS + EXON_READS))

    FRACTION_INTRON=$(awk -v i=$INTRON_READS -v t=$TOTAL_READS 'BEGIN {if(t==0) print 0; else print i/t}')

    echo -e "$SAMPLE\t$INTRON_READS\t$EXON_READS\t$TOTAL_READS\t$FRACTION_INTRON" >> "$OUTFILE"
done

echo "All samples processed. Results saved to $OUTFILE"
