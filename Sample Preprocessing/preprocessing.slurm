#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=2
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12 
#SBATCH --output=%j-%x.out
#SBATCH --error=%j-%x.err

# Make sure trimmomatic, FASTQC, STAR, and samtools are installed
# Change directory paths to correct paths for your analysis
# Make sure you have a reference index build by genomegenerate with STAR, using your chosen .gtf and .fa files (GENOME_DIR) 

# Set correct paths
RAW_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/raw_data"
FILTERED_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/filtering"
TRIMMOMATIC_ADAPTERS="/mainfs/lyceum/ch5g20/trimmomatic_adapters/TruSeq3-PE.fa"
GENOME_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/GRCh38_index"
FASTA="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/Homo_sapiens.GRCh38.113.gtf"
FASTQC_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/FASTQC"
BAM_DIR="/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/DS/bams"

# Load required modules/ activate required environments for trimmomatic and STAR
source ~/.bashrc
conda activate MISO_env

# Make required directories if they don't already exist
mkdir -p "$FILTERED_DIR"
mkdir -p "$BAM_DIR"
mkdir -p "$FASTQC_DIR"


# Loop through all R1 files
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"

    echo "Processing $SAMPLE"

    # Trimmomatic
    trimmomatic PE -threads 4 -phred33 \
        "$R1" "$R2" \
        "$FILTERED_DIR/${SAMPLE}_R1_paired.fastq.gz" "$FILTERED_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
        "$FILTERED_DIR/${SAMPLE}_R2_paired.fastq.gz" "$FILTERED_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Remove unpaired
    rm "$FILTERED_DIR/${SAMPLE}_R1_unpaired.fastq.gz" "$FILTERED_DIR/${SAMPLE}_R2_unpaired.fastq.gz"

    # FASTQC
    fastqc -o "$FASTQC_DIR" "$FILTERED_DIR/${SAMPLE}_R1_paired.fastq.gz" "$FILTERED_DIR/${SAMPLE}_R2_paired.fastq.gz"

    # STAR alignment
    STAR --runThreadN 8 \
        --genomeDir "$GENOME_DIR" \
        --readFilesCommand zcat \
        --readFilesIn "$FILTERED_DIR/${SAMPLE}_R1_paired.fastq.gz" "$FILTERED_DIR/${SAMPLE}_R2_paired.fastq.gz" \
        --outFileNamePrefix "$BAM_DIR/${SAMPLE}_" \
        --outSAMtype BAM SortedByCoordinate

    # Index BAM
    samtools index "$BAM_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam"
done
