#Load required libraries
library(ggplot2)

#create dataframe containing detection numbers per tool per dataset
detection_numbers <- data.frame(tool = c("rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2", "rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2", "rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2", "rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2", "rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2", "rMATS", "MAJIQ", "whippet", "leafcutter", "DEXSeq", "SUPPA2"),
                                dataset = c("HD2 (3,3)", "HD2 (3,3)", "HD2 (3,3)", "HD2 (3,3)", "HD2 (3,3)", "HD2 (3,3)", "FTD (4,4)", "FTD (4,4)", "FTD (4,4)", "FTD (4,4)", "FTD (4,4)", "FTD (4,4)", "HD (7,7)", "HD (7,7)", "HD (7,7)", "HD (7,7)", "HD (7,7)", "HD (7,7)", "DS (11,9)", "DS (11,9)", "DS (11,9)", "DS (11,9)", "DS (11,9)", "DS (11,9)", "DM (21,8)", "DM (21,8)", "DM (21,8)", "DM (21,8)", "DM (21,8)", "DM (21,8)", "ASD (12,12)", "ASD (12,12)", "ASD (12,12)", "ASD (12,12)", "ASD (12,12)", "ASD (12,12)"),
                                detection_number = c(3363, 8007, 5142, 127, 149, 1174, 2504, 5875, 4512, 181, 62, 489, 3034, 5652, 1743, 203, 354, 475, 5075, 6791, 3757, 2324, 1407, 253, 4090, 4133, 1522, 1957, 1793, 316, 857, 2555, 1999, 54, 2101, 218))

detection_numbers$dataset <- factor(detection_numbers$dataset, 
                                    levels = c("HD2 (3,3)", "FTD (4,4)", "HD (7,7)", "DS (11,9)", "DM (21,8)", "ASD (12,12)"))

#Create bar plot
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/detection_number_comparison_colours_reordered_with_asd.png", width = 1800, height = 900, res = 150)
ggplot(detection_numbers, aes(x = dataset, y = detection_number, fill = tool)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = detection_number),
            position = position_dodge(width = 0.75),
            vjust = -0.3, size = 3) +
  labs(title = "AS Event Detection No. Per Tool",
       x = "Dataset + Sample No.",
       y = "AS Event Detection No.") +
  scale_fill_manual(values = c(
    "rMATS" = "#440154FF",
    "MAJIQ" = "#414487FF",
    "whippet" = "#2A788EFF",
    "leafcutter" = "#22A884FF",
    "DEXSeq" = "#7AD151FF",
    "SUPPA2" = "#FDE725FF"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8500)) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid = element_blank()
  )

dev.off()




# Custom tool colors
tool_colors <- c(
  "rMATS" = "#440154FF",
  "MAJIQ" = "#414487FF",
  "whippet" = "#2A788EFF",
  "leafcutter" = "#22A884FF",
  "DEXSeq" = "#7AD151FF",
  "SUPPA2" = "#FDE725FF"
)


detection_numbers$dataset <- factor(detection_numbers$dataset,
                                    levels = c("HD2 (3,3)", "FTD (4,4)", "HD (7,7)", "DS (11,9)", "DM (21,8)", "ASD (12,12)"))

#Create line graph
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/detection_number_comparison_linegraph_no_grid_reordered_with_asd.png", width = 1800, height = 900, res = 150)
ggplot(detection_numbers, aes(x = dataset, y = detection_number, group = tool, color = tool)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = tool_colors) +
  labs(
    x = "Dataset + Sample No.",
    y = "DS Event Detection No.",
    color = "Tool",
    title = "DS Event Detection by Tool and Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black")
  )
dev.off()
