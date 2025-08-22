#Run comparison of gene overlap script first
#Load in files
DEXSeq <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/DEXSeq_filtered_padj_0.05_dPSI_0.05.tsv", header = TRUE)
leafcutter <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/leafcutter_with_gene_ids.txt", header = TRUE)
MAJIQ <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/MAJIQ_filtered_cleaned_changing-pvalue-threshold.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "")
whippet <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/whippet_with_genes.diff", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
rMATS_A3SS <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/rMATS.A3SS.cleaned.txt", header = TRUE)
rMATS_A5SS <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/rMATS.A5SS.cleaned.txt", header = TRUE)
rMATS_MXE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/rMATS.MXE.cleaned.txt", header = TRUE)
rMATS_SE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/FTD_padj_0.05_dPSI_0.05/rMATS.SE.cleaned.txt", header = TRUE)
rMATS_RI <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/FTD_padj_0.05_dPSI_0.05/rMATS.RI.cleaned.txt", header = TRUE)
SUPPA_A3 <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_A3_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_A5 <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_A5_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_AF <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_AF_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_AL <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_AL_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_MX <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_MX_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_RI <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_RI_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)
SUPPA_SE <- read.table("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/significant_SE_events_named.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "", row.names = NULL)

#Specify output directory
output_dir <- ("~/Documents/Big_Data_Biology/Research_Project/AS_tool_results/MD_padj_0.05_dPSI_0.05/graphs")

#manipulate dataframes to add a column with event location in the format "4:17604588-17604667"

##DEXSeq
#Combine relevant columns to create event location  
DEXSeq$genomic_coord <- paste0(DEXSeq$genomicData.seqnames, ":", DEXSeq$genomicData.start, "-", DEXSeq$genomicData.end)
#Preview
head(DEXSeq$genomic_coord)
##leafcutter
#Split the status column into parts
parts <- strsplit(leafcutter$status, ":")
#Extract chr, start, and end from each part
leafcutter$genomic_coord <- sapply(parts, function(x) {
  paste0(x[1], ":", x[2], "-", x[3])
})
#Remove "chr" prefix
leafcutter$genomic_coord <- gsub("chr", "", leafcutter$genomic_coord)  # removes "chr" prefix and cluster info
#Preview
head(leafcutter$genomic_coord)
#Remove intermmediate files
remove(parts)

##MAJIQ
#Create genomic coordinates in the format "chr#:start-end"
MAJIQ$genomic_coord <- paste0("chr", MAJIQ$seqid, ":",
                              sub(".*:(\\d+-\\d+)$", "\\1", MAJIQ$lsv_id))
#Remove "chr" prefix
MAJIQ$genomic_coord <- gsub("chr", "", MAJIQ$genomic_coord)  #removes "chr" prefix and cluster info
#Preview
head(MAJIQ$genomic_coord)

##whippet
#Alteady formatted
#View
head(whippet$Node)

##rMATS
#Combine all into one dataframe
A3SS <- rMATS_A3SS[, c(2, 3, 4, 6, 7, 20, 23)]
colnames(A3SS) <- c("GeneID", "GeneSymbol", "Chromosome", "ExonStart", "ExonEnd", "FDR", "incleveldiff")
A5SS <- rMATS_A5SS[, c(2, 3, 4, 6, 7, 20, 23)]
colnames(A5SS) <- c("GeneID", "GeneSymbol", "Chromosome", "ExonStart", "ExonEnd", "FDR", "incleveldiff")
SE   <- rMATS_SE[,   c(2, 3, 4, 6, 7, 20, 23)]
colnames(SE) <- c("GeneID", "GeneSymbol", "Chromosome", "ExonStart", "ExonEnd", "FDR", "incleveldiff")
RI   <- rMATS_RI[,   c(2, 3, 4, 6, 7, 20, 23)]
colnames(RI) <- c("GeneID", "GeneSymbol", "Chromosome", "ExonStart", "ExonEnd", "FDR", "incleveldiff")
MXE  <- rMATS_MXE[,  c(2, 3, 4, 7, 8, 22, 25)]
colnames(MXE) <- c("GeneID", "GeneSymbol", "Chromosome", "ExonStart", "ExonEnd", "FDR", "incleveldiff")
#Combine all rMATS event types
rMATS_all <- rbind(A3SS, A5SS, MXE, RI, SE)
#Remove "chr" prefix
rMATS_all$Chromosome <- gsub("chr", "", rMATS_all$Chromosome)
#Now manipulate exon start column to account for the ExonStart_0base output from rMATS
# Adjust ExonStart to 1-based
rMATS_all$ExonStart <- rMATS_all$ExonStart + 1
#Combine relevant columns to create event location  
rMATS_all$genomic_coord <- paste0(
  rMATS_all$Chromosome, ":", rMATS_all$ExonStart, "-", rMATS_all$ExonEnd)
#Preview result
head(rMATS_all$genomic_coord)
#Remove intermmediate files
remove(A3SS)
remove(A5SS)
remove(SE)
remove(MXE)
remove(RI)
remove(rMATS_A3SS)
remove(rMATS_A5SS)
remove(rMATS_MXE)
remove(rMATS_SE)
remove(rMATS_RI)


##SUPPA2
#change column names to match format 
colnames(SUPPA_A3) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_A5) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_AF) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_AL) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_RI) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_MX) <- c("row.names", "X", "gene_id", "dPSI", "pval")
colnames(SUPPA_SE) <- c("row.names", "X", "gene_id", "dPSI", "pval")
#Combine all rMATS event types
SUPPA2_all <- rbind(SUPPA_A3, SUPPA_A5 , SUPPA_AF, SUPPA_AL, SUPPA_MX, SUPPA_RI, SUPPA_SE)
#Create gene ID column
# Extract chromosome, start, and end from EventID
SUPPA_A3$genomic_coord <- sub(".*A3:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_A3$gene_id)

SUPPA_A5$genomic_coord <- sub(".*A5:([^:]+):([0-9]+)-([0-9]+):.*", 
                                "chr\\1:\\2-\\3", 
                                SUPPA_A5$gene_id)

SUPPA_AF$genomic_coord <- sub(".*AF:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_AF$gene_id)

SUPPA_AL$genomic_coord <- sub(".*AL:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_AL$gene_id)

SUPPA_SE$genomic_coord <- sub(".*SE:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_SE$gene_id)

SUPPA_MX$genomic_coord <- sub(".*MX:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_MX$gene_id)

SUPPA_RI$genomic_coord <- sub(".*RI:([^:]+):([0-9]+)-([0-9]+):.*", 
                              "chr\\1:\\2-\\3", 
                              SUPPA_RI$gene_id)

SUPPA2_all <- rbind(SUPPA_A3, SUPPA_A5 , SUPPA_AF, SUPPA_AL, SUPPA_MX, SUPPA_RI, SUPPA_SE)

#Remove "chr" prefix
SUPPA2_all$genomic_coord <- gsub("chr", "", SUPPA2_all$genomic_coord)
# Preview result
head(SUPPA2_all$genomic_coord)
#Remove intermediate files
remove(SUPPA_A3)
remove(SUPPA_A5)
remove(SUPPA_AF)
remove(SUPPA_AL)
remove(SUPPA_MX)
remove(SUPPA_SE)
remove(SUPPA_RI)



###extract event locations from each file 
DEXSeq_locations <- DEXSeq[, 41]
leafcutter_locations <- leafcutter[, 9]
whippet_locations <- whippet[, 3]  
rMATS_locations <- rMATS_all[, 8]  
MAJIQ_locations <- MAJIQ[, 19]  
SUPPA2_locations <- SUPPA2_all[, 6]



#Compare the overlap 
library(UpSetR)

locations_lists <- list(rMATS = rMATS_locations, leafcutter = leafcutter_locations, MAJIQ = MAJIQ_locations, whippet = whippet_locations, DEXSeq = DEXSeq_locations, SUPPA2 = SUPPA2_locations)
upset_input <- fromList(locations_lists)

png(filename = file.path(output_dir, "coordinate_overlap_upsetplot.png"), width = 1500, height = 800, res = 150)
upset(
  upset_input,
  nsets = 6,
  nintersects = 35,       # adjust depending on how many intersections you want shown
  order.by = "freq",      # order intersections by number of genes
  keep.order = TRUE,       # keep your set order
  sets.bar.color = "#db9925",           # Change set bar color
  main.bar.color = "#062379",           # Change intersection bar color
  matrix.color = "black",           # Change dot color
  point.size = 4,                       # Enlarge dots
  line.size = 1                         # Thicker lines
)
dev.off()

png(filename = file.path(output_dir, "coordinate_overlap_upsetplot_opt3.png"), width = 1000, height = 800, res = 150)
upset(
  upset_input,
  nsets = 6,
  nintersects = 40,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = "#db9925",       # Change set bar color
  main.bar.color = "#062379",       # Change intersection bar color
  matrix.color = "black",           # Change dot color
  point.size = 5,                   # Enlarge dots
  line.size = 1,                    # Thicker lines
  text.scale = c(2, 2, 1.8, 1.8, 1.5, 1.5) # Adjust text sizes
)
dev.off()

#From upset plot, look at which of the event locations overlap between certain tools
DEXSeq_whippet_overlap <- Reduce(intersect, list(DEXSeq_locations, whippet_locations))
DEXSeq_MAJIQ_overlap <- Reduce(intersect, list(DEXSeq_locations, MAJIQ_locations))
DEXSeq_rMATS_overlap <- Reduce(intersect, list(DEXSeq_locations, rMATS_locations))
MAJIQ_rMATS_whippet_overlap <- Reduce(intersect, list(whippet_locations, MAJIQ_locations, rMATS_locations))
MAJIQ_rMATS_overlap <- Reduce(intersect, list(MAJIQ_locations, rMATS_locations))
rMATS_whippet_overlap <- Reduce(intersect, list(rMATS_locations, whippet_locations))
leafcutter_SUPPA2_overlap <- Reduce(intersect, list(leafcutter_locations, SUPPA2_locations))

####Visualise overlapping events
###eulerr
library(eulerr)
# Deduplicate each gene set
rMATS_locations_unique     <- unique(rMATS_locations)
leafcutter_locations_unique <- unique(leafcutter_locations)
MAJIQ_locations_unique     <- unique(MAJIQ_locations)
whippet_locations_unique   <- unique(whippet_locations)
DEXSeq_locations_unique    <- unique(DEXSeq_locations)
SUPPA2_locations_unique     <- unique(SUPPA2_locations)
#Prep gene lists
fit <- euler(list(rMATS = rMATS_locations_unique, 
                  leafcutter = leafcutter_locations_unique, 
                  MAJIQ = MAJIQ_locations_unique, 
                  whippet = whippet_locations_unique, 
                  DEXSeq = DEXSeq_locations_unique, 
                  SUPPA2 = SUPPA2_locations_unique))
#Plot overlapping events 
png(filename = file.path(output_dir, "eulerr_event_overlap.png"), width = 1100, height = 1000, res = 150)
plot(fit, quantities = TRUE)
dev.off()

svg(filename = file.path(output_dir, "eulerr_event_overlap.svg"), width = 7, height = 7)
plot(fit, quantities = list(cex = 0.7), labels = list(cex = 0.8))
dev.off()

#remove intermediate files
remove(fit)
remove(SUPPA2_locations_unique)
remove(DEXSeq_locations_unique)     
remove(MAJIQ_locations_unique)
remove(rMATS_locations_unique)
remove(leafcutter_locations_unique)
remove(whippet_locations_unique)




#Find events unique to DEXSeq
all_other_locations <- Reduce(union, list(whippet_locations, rMATS_locations, MAJIQ_locations, leafcutter_locations, SUPPA2_locations))
unique_DEXSeq_locations <- setdiff(DEXSeq_locations, all_other_locations)
DEXSeq_unique <- DEXSeq[DEXSeq$genomic_coord %in% unique_DEXSeq_locations, ]

#Find events unique to rMATS
all_other_locations <- Reduce(union, list(whippet_locations, DEXSeq_locations, MAJIQ_locations, leafcutter_locations, SUPPA2_locations))
unique_rMATS_locations <- setdiff(rMATS_locations, all_other_locations)
rMATS_unique <- rMATS_all[rMATS_all$genomic_coord %in% unique_rMATS_locations, ]

#Find events unique to whippet
all_other_locations <- Reduce(union, list(rMATS_locations, DEXSeq_locations, MAJIQ_locations, leafcutter_locations, SUPPA2_locations))
unique_whippet_locations <- setdiff(whippet_locations, all_other_locations)
whippet_unique <- whippet[whippet$Node %in% unique_whippet_locations, ]

#Find events unique to SUPPA2
all_other_locations <- Reduce(union, list(whippet_locations, DEXSeq_locations, MAJIQ_locations, leafcutter_locations, rMATS_locations))
unique_SUPPA2_locations <- setdiff(SUPPA2_locations, all_other_locations)
SUPPA2_unique <- SUPPA2_all[SUPPA2_all$genomic_coord %in% unique_SUPPA2_locations, ]

#Find events unique to leafcutter
all_other_locations <- Reduce(union, list(whippet_locations, DEXSeq_locations, MAJIQ_locations, rMATS_locations, SUPPA2_locations))
unique_leafcutter_locations <- setdiff(leafcutter_locations, all_other_locations)
leafcutter_unique <- leafcutter[leafcutter$genomic_coord %in% unique_leafcutter_locations, ]

#Find events unique to MAJIQ
all_other_locations <- Reduce(union, list(whippet_locations, DEXSeq_locations, rMATS_locations, leafcutter_locations, SUPPA2_locations))
unique_MAJIQ_locations <- setdiff(MAJIQ_locations, all_other_locations)
MAJIQ_unique <- MAJIQ[MAJIQ$genomic_coord %in% unique_MAJIQ_locations, ]


#Sort files by highest dPSI
###MAJIQ
# Step 1: Get the maximum absolute dPSI per row
max_abs_dpsi <- apply(MAJIQ_unique["mean_dpsi_per_lsv_junction"], 1, function(x) {
  values <- as.numeric(unlist(strsplit(x, ";")))
  max(abs(values), na.rm = TRUE)
})
# Step 2: Add this as a new column (optional, for inspection)
MAJIQ_unique$max_abs_dpsi <- max_abs_dpsi
# Step 3: Sort the data frame by this value (descending)
MAJIQ_sorted <- MAJIQ_unique[order(-MAJIQ_unique$max_abs_dpsi), ]
MAJIQ_unique_sorted <- MAJIQ_sorted[1:10, c(2, 1, 19, 4)]
###rMATS
# Step 1 (optional): create a column for absolute dPSI
rMATS_unique$abs_incleveldiff <- abs(rMATS_unique$incleveldiff)
# Step 2: Sort descending by that
rMATS_sorted <- rMATS_unique[order(-rMATS_unique$abs_incleveldiff), ]
rMATS_unique_sorted <- rMATS_sorted[1:10, c(1, 2, 7, 6)]

###leafcutter
# Step 1: Calculate absolute dPSI
leafcutter_unique$abs_dpsi <- abs(leafcutter_unique$p.adjust)
# Step 2: Sort by that (descending)
leafcutter_sorted <- leafcutter_unique[order(-leafcutter_unique$abs_dpsi), ]
leafcutter_unique_sorted <- leafcutter_sorted[1:10, c(1, 8, 9, 7)]

###SUPPA2
# Step 1: Calculate absolute dPSI
SUPPA2_unique$abs_dpsi <- abs(SUPPA2_unique$dPSI)
# Step 2: Sort by that (descending)
SUPPA2_sorted <- SUPPA2_unique[order(-SUPPA2_unique$abs_dpsi), ]
SUPPA2_unique_sorted <- SUPPA2_sorted[1:10, c(2, 1, 6, 4)]
# Remove everything before the first colon, including semicolon
SUPPA2_unique_sorted$genomic_coord <- sub(".*?;", "", SUPPA2_unique_sorted$genomic_coord)

###DEXSeq
# Step 1: Calculate absolute dPSI
DEXSeq_unique$abs_dpsi <- abs(DEXSeq_unique$dPSI)
# Step 2: Sort by that (descending)
DEXSeq_sorted <- DEXSeq_unique[order(-DEXSeq_unique$abs_dpsi), ]
DEXSeq_unique_sorted <- DEXSeq_sorted[1:10, c(2, 1, 45, 49)]

###whippet
# Step 1: Calculate absolute dPSI
whippet_unique$abs_dpsi <- abs(whippet_unique$Psi_B)
# Step 2: Sort by that (descending)
whippet_sorted <- whippet_unique[order(-whippet_unique$abs_dpsi), ]
whippet_unique_sorted <- whippet_sorted[1:10, c(1, 13, 3, 5, 8)]
colnames(whippet_unique_sorted) <- c("Gene_ID", "Gene_Name", "Genomic_Coord", "Event_Type", "dPSI")



#Heatmap to look at overlap of tools (genes identified by at least 4 of them)
# Assume these are vectors of gene names
gene_lists <- list(
  DEXSeq = unique(DEXSeq_locations),
  leafcutter = unique(leafcutter_locations),
  MAJIQ = unique(MAJIQ_locations),
  whippet = unique(whippet_locations),
  rMATS = unique(rMATS_locations),
  SUPPA2 = unique(SUPPA2_locations)
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
png(filename = file.path(output_dir, "heatmap_event_overlap_with_geneids.png"), width = 1100, height = 2000, res = 150)
pheatmap(filtered_numeric,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Events Detected by ≥1 Tools",
         color = my_colours)
dev.off()





#Find common events between each pair of tools
common_rMATS_SUPPA2 <- Reduce(intersect, list(rMATS_locations, SUPPA2_locations))
common_rMATS_MAJIQ <- Reduce(intersect, list(rMATS_locations, MAJIQ_locations))
common_rMATS_leafcutter <- Reduce(intersect, list(rMATS_locations, leafcutter_locations))
common_rMATS_DEXSeq <- Reduce(intersect, list(rMATS_locations, DEXSeq_locations))
common_rMATS_whippet <- Reduce(intersect, list(rMATS_locations, whippet_locations))
common_SUPPA2_MAJIQ <- Reduce(intersect, list(SUPPA2_locations, MAJIQ_locations))
common_SUPPA2_leafcutter <- Reduce(intersect, list(SUPPA2_locations, leafcutter_locations))
common_SUPPA2_DEXSeq <- Reduce(intersect, list(SUPPA2_locations, DEXSeq_locations))
common_SUPPA2_whippet <- Reduce(intersect, list(SUPPA2_locations, whippet_locations))
common_MAJIQ_leafcutter <- Reduce(intersect, list(MAJIQ_locations, leafcutter_locations))
common_MAJIQ_DEXSeq <- Reduce(intersect, list(MAJIQ_locations, DEXSeq_locations))
common_MAJIQ_whippet <- Reduce(intersect, list(MAJIQ_locations, whippet_locations))
common_leafcutter_DEXSeq <- Reduce(intersect, list(leafcutter_locations, DEXSeq_locations))
common_leafcutter_whippet <- Reduce(intersect, list(leafcutter_locations, whippet_locations))
common_whippet_DEXSeq <- Reduce(intersect, list(whippet_locations, DEXSeq_locations))
