#Direct R to correct library directory
.libPaths("/mainfs/scratch/ch5g20/RESEARCH_PROJECT/R")

#Load required packages 
packages <- c("DEXSeq", "Rsubread", "tibble", "dplyr", "EnhancedVolcano", "edgeR", "GenomicFeatures", "GenomicAlignments", "Rsamtools", "BiocParallel")
invisible(lapply(packages, library, character.only = TRUE))

#Specifiy correct directory paths
gtf_file <- "/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/reference/GRCh38.p14/Homo_sapiens.GRCh38.113.gtf"
output_dir <- "/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/AS_analysis/DEXSeq"
bamDir <- "/mainfs/scratch/ch5g20/RESEARCH_PROJECT/COMPARISON/HD/bams"

#List names of input bam files in order (where all samples from the first sample group are listed together followed by all from the second sample group)
sampleOrder <- c("SRR3306823_Aligned.sortedByCoord.out.bam", "SRR3306824_Aligned.sortedByCoord.out.bam", "SRR3306825_Aligned.sortedByCoord.out.bam", "SRR3306826_Aligned.sortedByCoord.out.bam", "SRR3306827_Aligned.sortedByCoord.out.bam", "SRR3306828_Aligned.sortedByCoord.out.bam", "SRR3306829_Aligned.sortedByCoord.out.bam", "SRR3306830_Aligned.sortedByCoord.out.bam", "SRR3306831_Aligned.sortedByCoord.out.bam", "SRR3306832_Aligned.sortedByCoord.out.bam", "SRR3306833_Aligned.sortedByCoord.out.bam", "SRR3306834_Aligned.sortedByCoord.out.bam", "SRR3306835_Aligned.sortedByCoord.out.bam", "SRR3306836_Aligned.sortedByCoord.out.bam")

#Create a transcript database from the GTF
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

#Flatten into exonic parts (counting bins)
flattenedAnnotation <- exonicParts(txdb, linked.to.single.gene.only = TRUE)

#Set proper naming convention for exon IDs
names(flattenedAnnotation) <- sprintf(
  "%s:E%03d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part
)


#List bam files in specfied order
bamFiles <- file.path(bamDir, sampleOrder)

#Wrap into BamFileList (this enables parallel processing too)
bamList <- BamFileList(bamFiles, yieldSize=2000000)

#Set seqlevels to match (e.g., "UCSC" or "Ensembl", depending on STAR output)
seqlevelsStyle(flattenedAnnotation) <- "Ensembl"  # or "UCSC" if needed

#Adjust parameters based on your experiment (strand-specific, PE/SE)
se <- summarizeOverlaps(
  features = flattenedAnnotation,
  reads = bamList,
  mode = "Union",
  singleEnd = FALSE,            # set to TRUE if single-end
  fragments = TRUE,             # important for PE data
  ignore.strand = TRUE,         # change if your protocol is strand-specific
  BPPARAM = MulticoreParam(6)   # use your available cores
)

#Extract sample names from BAM filenames (remove file extension)
sampleNames <- tools::file_path_sans_ext(basename(bamFiles))

#Add to colData
colnames(se) <- sampleNames
colData(se)$sample <- factor(sampleNames)

#Construct your sample metadata table
colData(se)$condition <- factor(c("Con1", "Con1", "Con1", "Con1", "Con1", "Con1", "Con1",
                                  "Con2", "Con2", "Con2", "Con2", "Con2", "Con2", "Con2"))

colData(se)$libType <- factor(rep("paired-end", length(sampleNames)))


#Save intermediate files for potential troubleshooting
#Save count data 
count_matrix <- assay(se)
write.csv(as.data.frame(count_matrix), file = file.path(output_dir, "DEXSeq_raw_counts.csv"))
# Save metadata
metadata_table <- as.data.frame(colData(se))
write.csv(metadata_table, file = file.path(output_dir, "DEXSeq_metadata.csv"), row.names = TRUE)


#Run differental exon usage analysis with DEXSeq
dxd <- DEXSeqDataSetFromSE(se, design = ~ sample + exon + condition:exon)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)

dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
dxr <- DEXSeqResults(dxd)

#Save result object
save(dxd, dxr, file = file.path(output_dir, "DEXSeq_results.RData"))


#Extract results
results <- as.data.frame(dxr)

#Drop the 'transcripts' column (which prevents saving as a table)
results_clean <- results[, names(results) != "transcripts"]

#Write to TSV
write.table(results_clean, file = file.path(output_dir, "DEXSeq_exon_results_unfiltered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

#Filter results
filtered_results <- subset(results_clean, 
                           padj < 0.05 & 
                           (log2fold_Con2_Con1 < -0.1 | log2fold_Con2_Con1 > 0.1))

#Write to TSV
write.table(filtered_results, file = file.path(output_dir, "DEXSeq_exon_results_filtered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
