#Convert dPSI columns to names incidating tool
#Convert columns with genomic coordiates to "genomic_coord"
colnames(DEXSeq)[colnames(DEXSeq) == "dPSI"] <- "dPSI_DEXSeq"
colnames(rMATS_all)[colnames(rMATS_all) == "incleveldiff"] <- "dPSI_rMATS"
colnames(whippet)[colnames(whippet) == "Node"] <- "genomic_coord"
colnames(whippet)[colnames(whippet) == "Psi_B"] <- "dPSI_whippet"
colnames(SUPPA2_all)[colnames(SUPPA2_all) == "dPSI"] <- "dPSI_SUPPA2"
colnames(leafcutter)[colnames(leafcutter) == "p.adjust"] <- "dPSI_leafcutter"
colnames(MAJIQ)[colnames(MAJIQ) == "mean_dpsi_per_lsv_junction"] <- "dPSI_MAJIQ"


#Merge dataframes and keep rows with only matching AS events for each AS tool pair
rMATS_SUPPA_common <- merge(rMATS_all, SUPPA2_all, by = "genomic_coord")
rMATS_MAJIQ_common <- merge(rMATS_all, MAJIQ, by = "genomic_coord")
rMATS_leafcutter_common <- merge(rMATS_all, leafcutter, by = "genomic_coord")
rMATS_DEXSeq_common <- merge(rMATS_all, DEXSeq, by = "genomic_coord")
rMATS_whippet_common <- merge(rMATS_all, whippet, by = "genomic_coord")

SUPPA2_MAJIQ_common <- merge(SUPPA2_all, MAJIQ, by = "genomic_coord")
SUPPA2_leafcutter_common <- merge(SUPPA2_all, leafcutter, by = "genomic_coord")
SUPPA2_DEXSeq_common <- merge(SUPPA2_all, DEXSeq, by = "genomic_coord")
SUPPA2_whippet_common <- merge(SUPPA2_all, whippet, by = "genomic_coord")

MAJIQ_leafcutter_common <- merge(MAJIQ, leafcutter, by = "genomic_coord")
MAJIQ_DEXSeq_common <- merge(MAJIQ, DEXSeq, by = "genomic_coord")
MAJIQ_whippet_common <- merge(MAJIQ, whippet, by = "genomic_coord")

leafcutter_DEXSeq_common <- merge(leafcutter, DEXSeq, by = "genomic_coord")
leafcutter_whippet_common <- merge(leafcutter, whippet, by = "genomic_coord")

DEXSeq_whippet_common <- merge(DEXSeq, whippet, by = "genomic_coord")


#Extract event and dPSI values
rMATS_SUPPA_common <- rMATS_SUPPA_common[, c(1, 8, 12)]
rMATS_MAJIQ_common <- rMATS_MAJIQ_common[, c(1, 8, 12)] ##
rMATS_leafcutter_common <- rMATS_leafcutter_common[, c(1, 8, 15)]
rMATS_DEXSeq_common <- rMATS_DEXSeq_common[, c(1, 8, 51)] ##
rMATS_whippet_common <- rMATS_whippet_common[, c(1, 8, 15)] ##

SUPPA2_MAJIQ_common <- SUPPA2_MAJIQ_common[, c(1, 5, 9)]
SUPPA2_leafcutter_common <- SUPPA2_leafcutter_common[, c(1, 5, 13)] ##
SUPPA2_DEXSeq_common <- SUPPA2_DEXSeq_common[, c(1, 5, 44)]
SUPPA2_whippet_common <- SUPPA2_whippet_common[, c(1, 5, 12)]

MAJIQ_leafcutter_common <- MAJIQ_leafcutter_common[, c(1, 5, 12)]
MAJIQ_DEXSeq_common <- MAJIQ_DEXSeq_common[, c(1, 5, 62)] ##
MAJIQ_whippet_common <- MAJIQ_whippet_common[, c(1, 5, 26)] ##

leafcutter_DEXSeq_common <- leafcutter_DEXSeq_common[, c(1, 7, 46)]
leafcutter_whippet_common <- leafcutter_whippet_common[, c(1, 7, 14)]

DEXSeq_whippet_common <- DEXSeq_whippet_common[, c(1, 48, 55)] ##

#Select only the largest abs value for MAJIQ
MAJIQ_DEXSeq_common$dPSI_MAJIQ <- sapply(
  strsplit(MAJIQ_DEXSeq_common$dPSI_MAJIQ, ";"),
  function(x) x <- as.numeric(x)[which.max(abs(as.numeric(x)))]
)
MAJIQ_leafcutter_common$dPSI_MAJIQ <- sapply(
  strsplit(MAJIQ_leafcutter_common$dPSI_MAJIQ, ";"),
  function(x) x <- as.numeric(x)[which.max(abs(as.numeric(x)))]
)
MAJIQ_whippet_common$dPSI_MAJIQ <- sapply(
  strsplit(MAJIQ_whippet_common$dPSI_MAJIQ, ";"),
  function(x) x <- as.numeric(x)[which.max(abs(as.numeric(x)))]
)
 rMATS_MAJIQ_common$dPSI_MAJIQ <- sapply(
  strsplit(rMATS_MAJIQ_common$dPSI_MAJIQ, ";"),
  function(x) x <- as.numeric(x)[which.max(abs(as.numeric(x)))]
)
SUPPA2_MAJIQ_common$dPSI_MAJIQ <- sapply(
  strsplit(SUPPA2_MAJIQ_common$dPSI_MAJIQ, ";"),
  function(x) x <- as.numeric(x)[which.max(abs(as.numeric(x)))]
)

rMATS_MAJIQ_common$dPSI_rMATS <- abs(rMATS_MAJIQ_common$dPSI_rMATS)
rMATS_MAJIQ_common$dPSI_MAJIQ <- abs(rMATS_MAJIQ_common$dPSI_MAJIQ)

rMATS_DEXSeq_common$dPSI_rMATS <- abs(rMATS_DEXSeq_common$dPSI_rMATS)
rMATS_DEXSeq_common$dPSI_DEXSeq <- abs(rMATS_DEXSeq_common$dPSI_DEXSeq)

rMATS_whippet_common$dPSI_rMATS <- abs(rMATS_whippet_common$dPSI_rMATS)
rMATS_whippet_common$dPSI_whippet <- abs(rMATS_whippet_common$dPSI_whippet)


SUPPA2_leafcutter_common$dPSI_SUPPA2 <- abs(SUPPA2_leafcutter_common$dPSI_SUPPA2)
SUPPA2_leafcutter_common$dPSI_leafcutter <- abs(SUPPA2_leafcutter_common$dPSI_leafcutter)


MAJIQ_DEXSeq_common$dPSI_MAJIQ <- abs(MAJIQ_DEXSeq_common$dPSI_MAJIQ)
MAJIQ_DEXSeq_common$dPSI_DEXSeq <- abs(MAJIQ_DEXSeq_common$dPSI_DEXSeq)

MAJIQ_whippet_common$dPSI_MAJIQ <- abs(MAJIQ_whippet_common$dPSI_MAJIQ)
MAJIQ_whippet_common$dPSI_whippet <- abs(MAJIQ_whippet_common$dPSI_whippet)


DEXSeq_whippet_common$dPSI_DEXSeq <- abs(DEXSeq_whippet_common$dPSI_DEXSeq)
DEXSeq_whippet_common$dPSI_whippet <- abs(DEXSeq_whippet_common$dPSI_whippet)




##### DEXSEQ and whippet

# Convert to numeric vectors if they aren't already
dpsi_dexseq <- as.numeric(DEXSeq_whippet_common$dPSI_DEXSeq)
dpsi_whippet <- as.numeric(DEXSeq_whippet_common$dPSI_whippet)

# Check normality of the differences
differences <- dpsi_dexseq - dpsi_whippet
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_dexseq, dpsi_whippet, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_dexseq, dpsi_whippet, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_dexseq)),
  DEXSeq = dpsi_dexseq,
  Whippet = dpsi_whippet
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("DEXSeq" = "#7AD151FF", "Whippet" = "#2A788EFF")

# Plot
png(filename = file.path(output_dir,"DEXSeq_whippet_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "DEXSeq vs Whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 0.1584

png(filename = file.path(output_dir,"DEXSeq_whippet_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "DEXSeq vs Whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()



##### MAJIQ and whippet

# Convert to numeric vectors if they aren't already
dpsi_MAJIQ <- as.numeric(MAJIQ_whippet_common$dPSI_MAJIQ)
dpsi_whippet <- as.numeric(MAJIQ_whippet_common$dPSI_whippet)

# Check normality of the differences
differences <- dpsi_MAJIQ - dpsi_whippet
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_MAJIQ, dpsi_whippet, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_MAJIQ, dpsi_whippet, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_MAJIQ)),
  MAJIQ = dpsi_MAJIQ,
  Whippet = dpsi_whippet
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("MAJIQ" = "#414487FF", "Whippet" = "#2A788EFF")

# Plot
png(filename = file.path(output_dir,"MAJIQ_whippet_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "MAJIQ vs Whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 0.9054

png(filename = file.path(output_dir,"MAJIQ_whippet_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "MAJIQ vs Whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()





##### MAJIQ and DEXSeq

# Convert to numeric vectors if they aren't already
dpsi_MAJIQ <- as.numeric(MAJIQ_DEXSeq_common$dPSI_MAJIQ)
dpsi_dexseq <- as.numeric(MAJIQ_DEXSeq_common$dPSI_DEXSeq)

# Check normality of the differences
differences <- dpsi_dexseq - dpsi_MAJIQ
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_dexseq, dpsi_MAJIQ, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_dexseq, dpsi_MAJIQ, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_MAJIQ)),
  DEXSeq = dpsi_dexseq,
  MAJIQ = dpsi_MAJIQ
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("MAJIQ" = "#414487FF", "DEXSeq" = "#7AD151FF")

# Plot
png(filename = file.path(output_dir,"MAJIQ_DEXSeq_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "MAJIQ vs DEXSeq",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 0.01563

png(filename = file.path(output_dir,"MAJIQ_DEXSeq_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "MAJIQ vs DEXSeq",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()





##### SUPPA2 and leafcutter

# Convert to numeric vectors if they aren't already
dpsi_SUPPA2 <- as.numeric(SUPPA2_leafcutter_common$dPSI_SUPPA2)
dpsi_leafcutter <- as.numeric(SUPPA2_leafcutter_common$dPSI_leafcutter)

# Check normality of the differences
differences <- dpsi_SUPPA2 - dpsi_leafcutter
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_SUPPA2, dpsi_leafcutter, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_SUPPA2, dpsi_leafcutter, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_SUPPA2)),
  SUPPA2 = dpsi_SUPPA2,
  leafcutter = dpsi_leafcutter
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("SUPPA2" = "#FDE725FF", "leafcutter" = "#22A884FF")

# Plot
png(filename = file.path(output_dir,"SUPPA2_leafcutter_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "SUPPA2 vs leafcutter",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 1.05e-06

png(filename = file.path(output_dir,"SUPPA2_leafcutter_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "SUPPA2 vs leafcutter",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()






##### rMATS and whippet

# Convert to numeric vectors if they aren't already
dpsi_rMATS <- as.numeric(rMATS_whippet_common$dPSI_rMATS)
dpsi_whippet <- as.numeric(rMATS_whippet_common$dPSI_whippet)

# Check normality of the differences
differences <- dpsi_rMATS - dpsi_whippet
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_rMATS, dpsi_whippet, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_rMATS, dpsi_whippet, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_rMATS)),
  rMATS = dpsi_rMATS,
  whippet = dpsi_whippet
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("rMATS" = "#440154FF", "whippet" = "#2A788EFF")

# Plot
png(filename = file.path(output_dir,"rMATS_whippet_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 0.277

png(filename = file.path(output_dir,"rMATS_whippet_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs whippet",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()





##### rMATS and DEXSeq

# Convert to numeric vectors if they aren't already
dpsi_rMATS <- as.numeric(rMATS_DEXSeq_common$dPSI_rMATS)
dpsi_dexseq <- as.numeric(rMATS_DEXSeq_common$dPSI_DEXSeq)

# Check normality of the differences
differences <- dpsi_rMATS - dpsi_dexseq
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_rMATS, dpsi_dexseq, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_rMATS, dpsi_whippet, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_rMATS)),
  rMATS = dpsi_rMATS,
  DEXSeq = dpsi_dexseq
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("rMATS" = "#440154FF", "DEXSeq" = "#7AD151FF")

# Plot
png(filename = file.path(output_dir,"rMATS_DEXSeq_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs DEXSeq",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 7.236e-05
  
  png(filename = file.path(output_dir,"rMATS_DEXSeq_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs DEXseq",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()





##### rMATS and MAJIQ

# Convert to numeric vectors if they aren't already
dpsi_rMATS <- as.numeric(rMATS_MAJIQ_common$dPSI_rMATS)
dpsi_MAJIQ <- as.numeric(rMATS_MAJIQ_common$dPSI_MAJIQ)

# Check normality of the differences
differences <- dpsi_rMATS - dpsi_MAJIQ
shapiro.test(differences)  # If p > 0.05, normality assumed

# Paired t-test
t.test(dpsi_rMATS, dpsi_MAJIQ, paired = TRUE)

# Or Wilcoxon signed-rank test
wilcox.test(dpsi_rMATS, dpsi_MAJIQ, paired = TRUE)


# Combine into long format for ggplot2
library(reshape2)
library(ggplot2)

# Prepare the data in long format
paired_df <- data.frame(
  Event = paste0("Event_", 1:length(dpsi_rMATS)),
  rMATS = dpsi_rMATS,
  MAJIQ = dpsi_MAJIQ
)
long_df <- melt(paired_df, id.vars = "Event", variable.name = "Tool", value.name = "dPSI")

# Custom color palette
colors <- c("rMATS" = "#440154FF", "MAJIQ" = "#414487FF")

# Plot
png(filename = file.path(output_dir,"rMATS_MAJIQ_overlapping_events.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs MAJIQ",
       y = expression(Delta*"PSI"),
       x = NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Example: result from your previous test
p_val <- 0.09424

png(filename = file.path(output_dir,"rMATS_MAJIQ_overlapping_events_stats.png"), width = 700, height = 800, res = 150)
ggplot(long_df, aes(x = Tool, y = dPSI, group = Event)) +
  geom_line(color = "gray60", size = 0.5, alpha = 0.6) +
  geom_point(aes(color = Tool), size = 5) +
  scale_color_manual(values = colors) +
  theme_classic(base_size = 14) +
  labs(title = "rMATS vs MAJIQ",
       y = expression(Delta*"PSI"),
       x = NULL) +
  annotate("text", x = 1.5, y = max(long_df$dPSI) + 0.05,
           label = paste0("Paired t-test: p = ", signif(p_val, 3)),
           size = 5, fontface = "italic") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  )
dev.off()
