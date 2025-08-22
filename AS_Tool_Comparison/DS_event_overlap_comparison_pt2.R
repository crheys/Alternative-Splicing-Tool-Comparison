# Define sample names for each group
#DS
control_samples <- c("SRR22776675", "SRR22776676", "SRR22776677", "SRR22776686", "SRR22776687", "SRR22776688", "SRR22776689", "SRR22776690", "SRR22776691", "SRR22776692", "SRR22776693")
condition_samples <- c("SRR22776660", "SRR22776678", "SRR22776679", "SRR22776680", "SRR22776681", "SRR22776682", "SRR22776683", "SRR22776684", "SRR22776685")
#HD2
control_samples <- c("SRR22220826", "SRR22220827", "SRR22220828")
condition_samples <- c("SRR22220832", "SRR22220833", "SRR22220834")
#FTD
control_samples <- c("SRR8571937", "SRR8571941", "SRR8571944", "SRR8571947")
condition_samples <- c("SRR8571938", "SRR8571942", "SRR8571945", "SRR8571948")
#HD
control_samples <- c("SRR3306823", "SRR3306824", "SRR3306825", "SRR3306826", "SRR3306827", "SRR3306828", "SRR3306829")
condition_samples <- c("SRR3306830", "SRR3306831", "SRR3306832", "SRR3306833", "SRR3306834", "SRR3306835", "SRR3306836")
#DM
control_samples <- c("SRR12582120", "SRR12582121", "SRR12582122", "SRR12582123", "SRR12582124", "SRR12582125", "SRR12582126", "SRR12582127", "SRR12582128", "SRR12582129", "SRR12582130", "SRR12582131", "SRR12582132", "SRR12582133", "SRR12582134", "SRR12582135", "SRR12582136", "SRR12582137", "SRR12582138", "SRR12582139", "SRR12582140")
condition_samples <- c("SRR12582145", "SRR12582146", "SRR12582147", "SRR12582148", "SRR12582149", "SRR12582150", "SRR12582151", "SRR12582152")
#ASD
control_samples <- c("SRR5938419", "SRR5938420", "SRR5938421", "SRR5938422", "SRR5938423", "SRR5938424", "SRR5938425", "SRR5938426", "SRR5938427", "SRR5938428", "SRR5938429", "SRR5938430")
condition_samples <- c("SRR5938456", "SRR5938457", "SRR5938458", "SRR5938459", "SRR5938460", "SRR5938461", "SRR5938462", "SRR5938463", "SRR5938464", "SRR5938465", "SRR5938466", "SRR5938467")
# Get only the count columns
count_cols <- grep("^countData\\.", colnames(DEXSeq), value = TRUE)

# Extract exon counts
counts <- DEXSeq[, count_cols]

# Rename for simplicity
colnames(counts) <- gsub("countData\\.|_Aligned.*", "", count_cols)

DEXSeq$avg_control <- rowMeans(counts[, control_samples], na.rm = TRUE)
DEXSeq$avg_condition <- rowMeans(counts[, condition_samples], na.rm = TRUE)

# Add pseudocount to avoid division by zero
pseudocount <- 1

DEXSeq$dPSI <- (DEXSeq$avg_condition + pseudocount) / 
  ((DEXSeq$avg_control + DEXSeq$avg_condition) + 2 * pseudocount) - 
  (DEXSeq$avg_control + pseudocount) / 
  ((DEXSeq$avg_control + DEXSeq$avg_condition) + 2 * pseudocount)


DEXSeq_event_dpsi <- DEXSeq[, c(1, 45, 48)]
SUPPA2_event_dpsi <- SUPPA2_all[, c(2, 6, 4)]
MAJIQ_event_dspi <- MAJIQ[, c(2, 19, 4)]
rMATS_event_dpsi <- rMATS_all[, c(1, 7, 6)]
whippet_event_dpsi <- whippet[, c(1, 3, 8)]
leafcutter_event_dpsi <- leafcutter[, c(1, 9, 7)]

#remove intermmediate files
remove(control_samples)
remove(condition_samples)
remove(count_cols)
remove(counts)
remove(pseudocount)

#rMATS overlapping genes with MAJIQ
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
rMATS_matched_MAJIQ <- rMATS_event_dpsi[rMATS_event_dpsi$genomic_coord %in% overlap_coords, ]
#rMATS overlapping genes with whippet
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
rMATS_matched_whippet <- rMATS_event_dpsi[rMATS_event_dpsi$genomic_coord %in% overlap_coords, ]
#rMATS overlapping genes with MAJIQ and whippet
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
rMATS_matched_whippet_and_MAJIQ <- rMATS_event_dpsi[rMATS_event_dpsi$genomic_coord %in% overlap_coords, ]
#rMATS overlapping genes with DEXSeq
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_rMATS_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
rMATS_matched_DEXSeq <- rMATS_event_dpsi[rMATS_event_dpsi$genomic_coord %in% overlap_coords, ]

#MAJIQ overlapping genes with rMATS
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_rMATS <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]
#MAJIQ overlapping genes with rMATS and whippet
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_rMATS_and_whippet <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_MAJIQ_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_DEXSeq <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]

#MAJIQ overlapping genes with rMATS
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_rMATS <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]
#MAJIQ overlapping genes with rMATS and whippet
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_rMATS_and_whippet <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]
#MAJIQ overlapping genes with DEXSeq
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_MAJIQ_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
MAJIQ_matched_DEXSeq <- MAJIQ_event_dspi[MAJIQ_event_dspi$genomic_coord %in% overlap_coords, ]

#whippet overlapping genes with MAJIQ and rMATS
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(MAJIQ_rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
whippet_matched_rMATS_and_MAJIQ <- whippet_event_dpsi[whippet_event_dpsi$Node %in% overlap_coords, ]
#whippet overlapping genes with DEXSeq
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
whippet_matched_DEXSeq <- whippet_event_dpsi[whippet_event_dpsi$Node %in% overlap_coords, ]
#whippet overlapping genes with rMATS
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(rMATS_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
whippet_matched_rMATS <- whippet_event_dpsi[whippet_event_dpsi$Node %in% overlap_coords, ]

#DEXSeq overlapping genes with whippet
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_whippet_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
DEXSeq_matched_whippet <- DEXSeq_event_dpsi[DEXSeq_event_dpsi$genomic_coord %in% overlap_coords, ]
#DEXSeq overlapping genes with MAJIQ
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_MAJIQ_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
DEXSeq_matched_MAJIQ <- DEXSeq_event_dpsi[DEXSeq_event_dpsi$genomic_coord %in% overlap_coords, ]
#DEXSeq overlapping genes with MAJIQ
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(DEXSeq_rMATS_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
DEXSeq_matched_rMATS <- DEXSeq_event_dpsi[DEXSeq_event_dpsi$genomic_coord %in% overlap_coords, ]

#SUPPA2 overlapping genes with leafcutter
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(leafcutter_SUPPA2_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
SUPPA2_matched_leafcutter <- SUPPA2_event_dpsi[SUPPA2_event_dpsi$genomic_coord %in% overlap_coords, ]
#leafcutter overlapping genes with SUPPA2
# Example: convert the string of coordinates into a proper character vector
overlap_coords <- unlist(strsplit(leafcutter_SUPPA2_overlap, " "))
# Filter rMATS_event_dpsi rows where genomic_coord matches any in overlap_coords
leafcutter_matched_SUPPA2 <- leafcutter_event_dpsi[leafcutter_event_dpsi$genomic_coord %in% overlap_coords, ]




library(viridis)
viridis(6)
#rMATS == #440154FF
#MAJIQ == #414487FF
#whippet == #2A788EFF
#leafcutter == #22A884FF
#DEXSeq == #7AD151FF
#SUPPA2 == #FDE725FF

#For MAJIQ and DEXSeq:
#First create column for MAJIQ which lists the most extreme absolute value (MAJIQ lists multiple dPSI values per event)
library(dplyr)
# Split and compute max absolute value
MAJIQ_matched_DEXSeq <- MAJIQ_matched_DEXSeq %>%
  mutate(
    MAJIQ_dpsi = sapply(strsplit(mean_dpsi_per_lsv_junction, ";"), function(x) {
      x_num <- as.numeric(x)
      x_num[which.max(abs(x_num))]  # max absolute dPSI
    })
  )
# Merge on genomic_coord (can also use gene_id if unique per row)
combined_dpsi <- merge(DEXSeq_matched_MAJIQ, MAJIQ_matched_DEXSeq[, c("gene_id", "genomic_coord", "MAJIQ_dpsi")],
                       by = "genomic_coord", suffixes = c("_DEXSeq", "_MAJIQ"))
library(tidyr)
library(dplyr)

# Rename for clarity
plot_data <- combined_dpsi %>%
  select(genomic_coord, dPSI, MAJIQ_dpsi) %>%
  rename(DEXSeq = dPSI, MAJIQ = MAJIQ_dpsi)

# Convert to long format
plot_data_long <- pivot_longer(
  plot_data,
  cols = c(DEXSeq, MAJIQ),
  names_to = "Tool",
  values_to = "dPSI"
)
library(ggplot2)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

png(filename = file.path(output_dir, "dexseq_and_MAJIQ_overlapping_events_dPSI.png"), width = 1200, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#7AD151FF", "#414487FF")) +
  labs(
    title = "Comparison of dPSI Values: DEXSeq vs MAJIQ",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()





#DEXSeq and whippet
# First, rename the columns for consistency
DEXSeq_matched_whippet_clean <- DEXSeq_matched_whippet %>%
  select(genomic_coord, DEXSeq = dPSI)

whippet_matched_DEXSeq_clean <- whippet_matched_DEXSeq %>%
  select(genomic_coord = Node, Whippet = Psi_B)

# Merge by genomic coordinate
combined_dpsi <- inner_join(DEXSeq_matched_whippet_clean, whippet_matched_DEXSeq_clean, by = "genomic_coord")

library(tidyr)

# Reshape to long format
plot_data_long <- pivot_longer(
  combined_dpsi,
  cols = c("DEXSeq", "Whippet"),
  names_to = "Tool",
  values_to = "dPSI"
)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

library(ggplot2)

png(filename = file.path(output_dir,"dexseq_and_whippet_overlapping_events_dPSI.png"), width = 1200, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#7AD151FF", "#2A788EFF")) +
  labs(
    title = "Comparison of dPSI Values: DEXSeq vs Whippet",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()






#For MAJIQ and rMATS:
#First create column for MAJIQ which lists the most extreme absolute value (MAJIQ lists multiple dPSI values per event)
library(dplyr)
# Split and compute max absolute value
MAJIQ_matched_rMATS <- MAJIQ_matched_rMATS %>%
  mutate(
    MAJIQ_dpsi = sapply(strsplit(mean_dpsi_per_lsv_junction, ";"), function(x) {
      x_num <- as.numeric(x)
      x_num[which.max(abs(x_num))]  # max absolute dPSI
    })
  )
# Merge on genomic_coord (can also use gene_id if unique per row)
combined_dpsi <- merge(rMATS_matched_MAJIQ, MAJIQ_matched_rMATS[, c("gene_id", "genomic_coord", "MAJIQ_dpsi")],
                       by = "genomic_coord", suffixes = c("_rMATS", "_MAJIQ"))
library(tidyr)
library(dplyr)

# Rename for clarity
plot_data <- combined_dpsi %>%
  select(genomic_coord, incleveldiff, MAJIQ_dpsi) %>%
  rename(rMATS = incleveldiff, MAJIQ = MAJIQ_dpsi)

# Convert to long format
plot_data_long <- pivot_longer(
  plot_data,
  cols = c(rMATS, MAJIQ),
  names_to = "Tool",
  values_to = "dPSI"
)
library(ggplot2)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

png(filename = file.path(output_dir,"rMATS_and_MAJIQ_overlapping_events_dPSI.png"), width = 2000, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#87ceeb", "#440192")) +
  labs(
    title = "Comparison of dPSI Values: rMATS vs MAJIQ",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()








#rMATS and whippet
# First, rename the columns for consistency
rMATS_matched_whippet_clean <- rMATS_matched_whippet %>%
  select(genomic_coord, rMATS = incleveldiff)

whippet_matched_rMATS_clean <- whippet_matched_rMATS %>%
  select(genomic_coord = Node, Whippet = Psi_B)

# Merge by genomic coordinate
combined_dpsi <- inner_join(rMATS_matched_whippet_clean, whippet_matched_rMATS_clean, by = "genomic_coord")

library(tidyr)

# Reshape to long format
plot_data_long <- pivot_longer(
  combined_dpsi,
  cols = c("rMATS", "Whippet"),
  names_to = "Tool",
  values_to = "dPSI"
)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

library(ggplot2)

png(filename = file.path(output_dir,"rMATS_and_whippet_overlapping_events_dPSI.png"), width = 1600, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#440154FF", "#2A788EFF")) +
  labs(
    title = "Comparison of dPSI Values: rMATS vs Whippet",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()



#SUPPA2 and leafcutter
# First, rename the columns for consistency
leafcutter_matched_SUPPA2_clean <- leafcutter_matched_SUPPA2 %>%
  select(genomic_coord, leafcutter = p.adjust)

SUPPA2_matched_leafcutter_clean <- SUPPA2_matched_leafcutter %>%
  select(genomic_coord, SUPPA2 = dPSI)

SUPPA2_matched_leafcutter_clean <- SUPPA2_matched_leafcutter_clean[c(2, 3, 4, 6, 7) ,]
# Merge by genomic coordinate
combined_dpsi <- inner_join(leafcutter_matched_SUPPA2_clean, SUPPA2_matched_leafcutter_clean, by = "genomic_coord")

library(tidyr)

# Reshape to long format
plot_data_long <- pivot_longer(
  combined_dpsi,
  cols = c("leafcutter", "SUPPA2"),
  names_to = "Tool",
  values_to = "dPSI"
)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

library(ggplot2)

png(filename = file.path(output_dir,"leafcutter_and_SUPPA2_overlapping_events_dPSI.png"), width = 1600, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#22A884FF", "#FDE725FF")) +
  labs(
    title = "Comparison of dPSI Values: SUPPA2 vs leafcutter",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()






#DEXSeq and rMATS
# First, rename the columns for consistency
DEXSeq_matched_rMATS_clean <- DEXSeq_matched_rMATS %>%
  select(genomic_coord, DEXSeq = dPSI)

rMATS_matched_DEXSeq_clean <- rMATS_matched_DEXSeq %>%
  select(genomic_coord, rMATS = incleveldiff)

# Merge by genomic coordinate
combined_dpsi <- inner_join(DEXSeq_matched_rMATS_clean, rMATS_matched_DEXSeq_clean, by = "genomic_coord")

library(tidyr)

# Reshape to long format
plot_data_long <- pivot_longer(
  combined_dpsi,
  cols = c("DEXSeq", "rMATS"),
  names_to = "Tool",
  values_to = "dPSI"
)

#Convert to absolute dPSI values
plot_data_long$dPSI <- abs(plot_data_long$dPSI)

library(ggplot2)

png(filename = file.path(output_dir,"DEXSeq_and_rMATS_overlapping_events_dPSI.png"), width = 1600, height = 800, res = 150)
ggplot(plot_data_long, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#7AD151FF", "#440154FF")) +
  labs(
    title = "Comparison of dPSI Values: rMATS vs DEXSeq",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()










#SUPPA2 and leafcutter
#Comparisom of dPSI for event detected by rMATS, whippet and MAJIQ
overlap <- data.frame(genomic_coord= c("1:26225703-26227404", "1:26225703-26227404", "1:3060283-3060501", "1:3060283-3060501", "14:102508549-102508641", "14:102508549-102508641", "16:979110-979601", "16:979110-979601", "2:147973524-147975902", "2:147973524-147975902", "9:127896332-127897956", "9:127896332-127897956", "X:48898093-48898492", "X:48898093-48898492"),
                              Tool = c("SUPPA2", "leafcutter"),
                              dPSI = c(0.28991791, 0.07032137, 0.15905283, 0.12867293, 0.17794743, 0.25556320, 0.34326697, 0.14646452, 0.16177138, 0.05442887, 0.08672121, 0.06674656, 0.15589548, 0.24987834)
)

#Convert to absolute dPSI values
overlap$dPSI <- abs(overlap$dPSI)

library(ggplot2)

png(filename = file.path(output_dir,"SUPPA2_leafcutter_overlapping_events_dPSI.png"), width = 1200, height = 800, res = 150)
ggplot(overlap, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#22A884FF", "#FDE725FF")) +
  labs(
    title = "Comparison of dPSI Values: SUPPA2 vs leafcutter",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()




#Comparisom of dPSI for event detected by rMATS, whippet and MAJIQ
ENSG00000060237 <- data.frame(genomic_coord= c("10:95375973-95376056", "10:95375973-95376056", "10:95375973-95376056", "3:33592213-33592239", "3:33592213-33592239", "3:33592213-33592239", "8:104068512-104068612", "8:104068512-104068612", "8:104068512-104068612"),
                              Tool = c("rMATS", "MAJIQ", "whippet"),
                              dPSI = c(0.214, 0.1199, 0.21112, 0.272, 0.3556, 0.20042, 0.108, 0.1649, 0.19885)
                              )

#Convert to absolute dPSI values
ENSG00000060237$dPSI <- abs(ENSG00000060237$dPSI)

library(ggplot2)

png(filename = file.path(output_dir, "dPSI_rMATS_and_whippet_and_MAJIQ_overlapping_events_comp.png"), width = 800, height = 800, res = 150)
ggplot(ENSG00000060237, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#440154FF", "#2A788EFF", "#414487FF")) +
  labs(
    title = "Comparison of dPSI Values: rMATS vs Whippet vs MAJIQ",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()




#Comparisom of dPSI for event detected by rMATS, whippet and MAJIQ
OVERLAP <- data.frame(genomic_coord= c("16:979110-979601", "16:979110-979601"),
                              Tool = c("leafcutter", "SUPPA2"),
                              dPSI = c(-0.1464645, 0.2451645)
)

#Convert to absolute dPSI values
OVERLAP$dPSI <- abs(OVERLAP$dPSI)

library(ggplot2)

png("~/Documents/Big_Data_Biology/Research_Project/HD_padj_0.05_dPSI_0.1/graphs/dPSI_leafcutter_and_SUPPA2_overlapping_events_comp_vir.png", width = 800, height = 800, res = 150)
ggplot(OVERLAP, aes(x = genomic_coord, y = dPSI, fill = Tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#22A884FF", "#FDE725FF")) +
  labs(
    title = "Comparison of dPSI Values: leafcutter vs SUPPA2",
    x = "Splicing Event (Genomic Coordinate)",
    y = "dPSI",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()
