#Load files in
DEXSeq <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/DEXSeq_exon_results_with_names_0.05.tsv", header = TRUE)
leafcutter <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/leafcutter_with_gene_ids_all_samples.txt", header = TRUE)
MAJIQ <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/MAJIQ_dPSI_0.05_filtered_cleaned_changing-pvalue-threshold.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "")
whippet <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/whippet_with_genes_dPSI_0.05.diff", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
rMATS_A3SS <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/rMATS.A3SS.cleaned.txt", header = TRUE)
rMATS_A5SS <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/rMATS.A5SS.cleaned.txt", header = TRUE)
rMATS_MXE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/rMATS.MXE.cleaned.txt", header = TRUE)
rMATS_SE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/rMATS.SE.cleaned.txt", header = TRUE)
rMATS_RI <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/rMATS.RI.cleaned.txt", header = TRUE)
SUPPA_A3 <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_A3.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_A5 <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_A5.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_AF <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_AF.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_AL <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_AL.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_MX <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_MX.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_RI <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_RI.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_SE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/SUPPA2_SE.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)

#Specify output directory
output_dir <- ("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/ASD_padj_0.05_dPSI_0.05/for_tool_comparison/graphs")

#Combine all gene IDs and gene names from the different rMATS event types into one file 
#Extract just columns 2 and 3 from each data set
A3SS_genes <- rMATS_A3SS[, 2:3]
A5SS_genes <- rMATS_A5SS[, 2:3]
SE_genes   <- rMATS_SE[,   2:3]
RI_genes   <- rMATS_RI[,   2:3]
MXE_genes  <- rMATS_MXE[,  2:3]

#Combine all into one dataframe
rMATS_all_genes <- rbind(A3SS_genes, A5SS_genes, SE_genes, RI_genes, MXE_genes)

SUPPA_A3_genes <- SUPPA_A3[, 1:2]
SUPPA_A5_genes <- SUPPA_A5[, 1:2]
SUPPA_AF_genes <- SUPPA_AF[, 1:2]
SUPPA_AL_genes <- SUPPA_AL[, 1:2]
SUPPA_MX_genes <- SUPPA_MX[, 1:2]
SUPPA_RI_genes <- SUPPA_RI[, 1:2]
SUPPA_SE_genes <- SUPPA_SE[, 1:2]

#combine all SUPPA genes
SUPPA_all_genes <- do.call(rbind, lapply(list(SUPPA_A3_genes, SUPPA_A5_genes,
                                              SUPPA_AF_genes, SUPPA_AL_genes,
                                              SUPPA_MX_genes, SUPPA_RI_genes,
                                              SUPPA_SE_genes), as.matrix))
SUPPA_all_genes <- as.data.frame(SUPPA_all_genes)


#Remove unneeded files
remove(A3SS_genes)
remove(A5SS_genes)
remove(SE_genes)
remove(RI_genes)
remove(MXE_genes)
remove(SUPPA_A3_genes)
remove(SUPPA_A5_genes)
remove(SUPPA_AF_genes)
remove(SUPPA_AL_genes)
remove(SUPPA_MX_genes)
remove(SUPPA_RI_genes)
remove(SUPPA_SE_genes)
remove(rMATS_A3SS) 
remove(rMATS_A5SS)
remove(rMATS_MXE) 
remove(rMATS_SE) 
remove(rMATS_RI) 
remove(SUPPA_A3) 
remove(SUPPA_A5)
remove(SUPPA_AF)
remove(SUPPA_AL)
remove(SUPPA_MX)
remove(SUPPA_RI)
remove(SUPPA_SE)

#Extract gene IDs and gene names from other tool files 
DEXSeq_genes <- DEXSeq[, 2]
leafcutter_genes <- leafcutter[, 1]
MAJIQ_genes <- MAJIQ[, 2]
whippet_genes <- whippet[, 1]
rMATS_genes <- rMATS_all_genes[, 1]
SUPPA_genes <- SUPPA_all_genes[, 2]
#Load required library
library(VennDiagram)
library(RColorBrewer)


###ggvenn
#Load required package for generating 6 group venn diagrams
library(ggvenn)
# Convert to character vectors (gene IDs)
DEXSeq_genes     <- as.character(na.omit(DEXSeq_genes))
leafcutter_genes <- as.character(na.omit(leafcutter_genes))
MAJIQ_genes      <- as.character(na.omit(MAJIQ_genes))
whippet_genes    <- as.character(na.omit(whippet_genes))
rMATS_genes      <- as.character(na.omit(rMATS_genes))
SUPPA_genes  <- as.character(na.omit(SUPPA_genes))
#List gene sets from each tool
gene_lists <- list(rMATS = rMATS_genes, leafcutter = leafcutter_genes, MAJIQ = MAJIQ_genes, whippet = whippet_genes, DEXSeq = DEXSeq_genes, SUPPA2 = SUPPA_genes)
#Create venn diagram
png(filename = file.path(output_dir, "ggvenn.png"), width = 1100, height = 1000, res = 150)
ggvenn(gene_lists)
dev.off()

###eulerr
library(eulerr)
# Deduplicate each gene set
rMATS_genes_unique     <- unique(rMATS_genes)
leafcutter_genes_unique <- unique(leafcutter_genes)
MAJIQ_genes_unique     <- unique(MAJIQ_genes)
whippet_genes_unique   <- unique(whippet_genes)
DEXSeq_genes_unique    <- unique(DEXSeq_genes)
SUPPA_genes_unique     <- unique(SUPPA_genes)
#Prep gene lists
fit <- euler(list(rMATS = rMATS_genes_unique, leafcutter = leafcutter_genes_unique, MAJIQ = MAJIQ_genes_unique, whippet = whippet_genes_unique, DEXSeq = DEXSeq_genes_unique, SUPPA2 = SUPPA_genes_unique))
#Plot Venn diagram 
png(filename = file.path(output_dir, "eulerr.png"), width = 1100, height = 1000, res = 150)
plot(fit, quantities = TRUE)
dev.off()

plot(fit, quantities = TRUE, labels = FALSE)  # Hides set labels
plot(fit, 
     quantities = list(cex = 0.9,   # Smaller font size for counts
     labels = list(cex = 0.1))       # Smaller font size for set names

svg(filename = file.path(output_dir, "eulerr.svg"), width = 7, height = 7)
plot(fit, quantities = list(cex = 0.7), labels = list(cex = 0.8))
dev.off()

#remove intermediate files
remove(rMATS_genes_unique)
remove(whippet_genes_unique)
remove(SUPPA_genes_unique)
remove(DEXSeq_genes_unique)
remove(MAJIQ_genes_unique)
remove(leafcutter_genes_unique)
remove(fit)




###UpSetR
#Load libraries
library(UpSetR)
#Prep gene lists
gene_lists <- list(rMATS = rMATS_genes, leafcutter = leafcutter_genes, MAJIQ = MAJIQ_genes, whippet = whippet_genes, DEXSeq = DEXSeq_genes, SUPPA2 = SUPPA_genes)
upset_input <- fromList(gene_lists)
#Generate Upset plot
png(filename = file.path(output_dir,"UpSetR_diff_colours.png"), width = 1800, height = 800, res = 150)
upset(
  upset_input,
  nsets = 6,
  nintersects = 100,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = "#db9925",       # Change set bar color
  main.bar.color = "#7f2054",           # Change intersection bar color
  matrix.color = "black",           # Change dot color
  point.size = 3,                       # Enlarge dots
  line.size = 1                         # Thicker lines
)
dev.off()
#colours used for MD dataset #db9925 and #7f205

png(filename = file.path(output_dir,"UpSetR_diff_larger_text.png"), width = 1800, height = 800, res = 150)
upset(
  upset_input,
  nsets = 6,
  nintersects = 100,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = "#db9925",       # Change set bar color
  main.bar.color = "#7f2054",       # Change intersection bar color
  matrix.color = "black",           # Change dot color
  point.size = 5,                   # Enlarge dots
  line.size = 1,                    # Thicker lines
  text.scale = c(2, 2, 1.8, 1.8, 1.5, 1.5) # Adjust text sizes
)
dev.off()

#Look at overlapping genes 
#Generate lists of overlapping genes
#Look at which genes are detected by all events
common_genes_all <- Reduce(intersect, list(DEXSeq_genes, leafcutter_genes, MAJIQ_genes, rMATS_genes, whippet_genes, SUPPA_genes))
#Look at genes detected by 5 of the tools
common_genes_not_whippet <- Reduce(intersect, list(DEXSeq_genes, leafcutter_genes, MAJIQ_genes, rMATS_genes, SUPPA_genes))
common_genes_not_rMATS <- Reduce(intersect, list(DEXSeq_genes, leafcutter_genes, MAJIQ_genes, whippet_genes, SUPPA_genes))
common_genes_not_MAJIQ <- Reduce(intersect, list(DEXSeq_genes, leafcutter_genes, rMATS_genes, whippet_genes, SUPPA_genes))
common_genes_not_leafcutter <- Reduce(intersect, list(DEXSeq_genes, MAJIQ_genes, rMATS_genes, whippet_genes, SUPPA_genes))
common_genes_not_DEXSeq <- Reduce(intersect, list(leafcutter_genes, MAJIQ_genes, rMATS_genes, whippet_genes, SUPPA_genes))
common_genes_not_SUPPA <- Reduce(intersect, list(leafcutter_genes, MAJIQ_genes, rMATS_genes, whippet_genes, DEXSeq_genes))

common_rMATS_SUPPA2 <- Reduce(intersect, list(rMATS_genes, SUPPA_genes))
common_rMATS_MAJIQ <- Reduce(intersect, list(rMATS_genes, MAJIQ_genes))
common_rMATS_leafcutter <- Reduce(intersect, list(rMATS_genes, leafcutter_genes))
common_rMATS_DEXSeq <- Reduce(intersect, list(rMATS_genes, DEXSeq_genes))
common_rMATS_whippet <- Reduce(intersect, list(rMATS_genes, whippet_genes))
common_SUPPA2_MAJIQ <- Reduce(intersect, list(SUPPA_genes, MAJIQ_genes))
common_SUPPA2_leafcutter <- Reduce(intersect, list(SUPPA_genes, leafcutter_genes))
common_SUPPA2_DEXSeq <- Reduce(intersect, list(SUPPA_genes, DEXSeq_genes))
common_SUPPA2_whippet <- Reduce(intersect, list(SUPPA_genes, whippet_genes))
common_MAJIQ_leafcutter <- Reduce(intersect, list(MAJIQ_genes, leafcutter_genes))
common_MAJIQ_DEXSeq <- Reduce(intersect, list(MAJIQ_genes, DEXSeq_genes))
common_MAJIQ_whippet <- Reduce(intersect, list(MAJIQ_genes, whippet_genes))
common_leafcutter_DEXSeq <- Reduce(intersect, list(leafcutter_genes, DEXSeq_genes))
common_leafcutter_whippet <- Reduce(intersect, list(leafcutter_genes, whippet_genes))
common_whippet_DEXSeq <- Reduce(intersect, list(whippet_genes, DEXSeq_genes))

#Heatmap to look at overlap of tools (genes identified by at least 4 of them)
# Assume these are vectors of gene names
gene_lists <- list(
  DEXSeq = unique(DEXSeq_genes),
  leafcutter = unique(leafcutter_genes),
  MAJIQ = unique(MAJIQ_genes),
  whippet = unique(whippet_genes),
  rMATS = unique(rMATS_genes),
  SUPPA2 = unique(SUPPA_genes)
)
# Get all unique genes across all lists
all_genes <- unique(unlist(gene_lists))
# Create a binary matrix: rows = genes, columns = lists
presence_matrix <- sapply(gene_lists, function(g) all_genes %in% g)
rownames(presence_matrix) <- all_genes
# Filter for genes present in 4 or more lists
filtered_matrix <- presence_matrix[rowSums(presence_matrix) >= 1, ]
# Convert logical matrix to numeric (1 = present, 0 = absent)
filtered_numeric <- apply(filtered_matrix, 2, as.numeric)
rownames(filtered_numeric) <- rownames(filtered_matrix)

# Install if needed
library(pheatmap)
library(viridis)
# Define your custom color gradient
my_colours <- colorRampPalette(c("#0D0887FF", "black", "#E16462FF"))(100)
# Plot with the custom colors
png(filename = file.path(output_dir, "heatmap_gene_overlap_with_geneids_1_opt2.png"), width = 1100, height = 2500, res = 150)
pheatmap(filtered_numeric,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Genes Detected by ≥1 Tools",
         color = my_colours,
         fontsize = 10,
         fontsize_col = 35,
         treeheight_col = 200)
dev.off()

png(filename = file.path(output_dir, "heatmap_gene_overlap_with_geneids_1.png"), width = 1100, height = 3000, res = 150)
pheatmap(filtered_numeric,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Genes Detected by ≥1 Tools",
         color = my_colours,
         fontsize = 10,
         fontsize_col = 35,
         treeheight_col = 200)
dev.off()


#remove intermediate files
remove(upset_input)
remove(all_genes)
remove(filtered_matrix)
remove(filtered_numeric)
remove(presence_matrix)
remove(gene_lists)
remove(common_genes_all)
remove(common_genes_not_DEXSeq)
remove(common_genes_not_MAJIQ)
remove(common_genes_not_leafcutter)
remove(common_genes_not_rMATS)
remove(common_genes_not_SUPPA)
remove(common_genes_not_whippet)
remove(my_colours)
