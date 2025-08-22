# Create a symmetric matrix of Jaccard indices
#HD
jaccard_matrix <- matrix(c(
  1.00, 0.026023392, 0.122657361, 0.014097744, 0.022021116, 0.033758927,
  0.026023392, 1.00, 0.021507169, 0.011940299, 0.012210012, 0.008181818,
  0.122657361, 0.021507169, 1.00, 0.010702572, 0.016071731, 0.02381282,
  0.014097744, 0.011940299, 0.010702572, 1.00, 0.018281536, 0.002059732,
  0.022021116, 0.012210012, 0.016071731, 0.018281536, 1.00, 0.008658009,
  0.033758927, 0.008181818, 0.02381282, 0.002059732, 0.008658009, 1.00
), nrow = 6, byrow = TRUE)

#DS
jaccard_matrix <- matrix(c(
  1.00, 0.011773642, 0.164474975, 0.078414225, 0.070166749, 0.032499415,
  0.011773642, 1.00, 0.013233602, 0.007427678, 0.021538462, 0.003754693,
  0.164474975, 0.013233602, 1.00, 0.074501945, 0.065782631, 0.031085044,
  0.078414225, 0.007427678, 0.074501945, 1.00, 0.03783032, 0.011140672,
  0.070166749, 0.021538462, 0.065782631, 0.03783032, 1.00, 0.019344651,
  0.032499415, 0.003754693, 0.031085044, 0.011140672, 0.019344651, 1.00
), nrow = 6, byrow = TRUE)

#DM
jaccard_matrix <- matrix(c(
  1.00, 0.022511023, 0.140974053, 0.082527748, 0.06576087, 0.027274391,
  0.022511023, 1.00, 0.025824303, 0.031306715, 0.029785156, 0.0137893,
  0.140974053, 0.025824303, 1.00, 0.0935536, 0.077454545, 0.018551873,
  0.082527748, 0.031306715, 0.0935536, 1.00, 0.050420168, 0.015766423,
  0.06576087, 0.029785156, 0.077454545, 0.050420168, 1.00, 0.023148148,
  0.027274391, 0.0137893, 0.018551873, 0.015766423, 0.023148148, 1.00
), nrow = 6, byrow = TRUE)

#ASD
jaccard_matrix <- matrix(c(
  1.00, 0.005834306, 0.080987029, 0.00770077, 0.036376355, 0.010240113,
  0.005834306, 1.00, 0.010200364, 0, 0.021136063, 0.002260398,
  0.080987029, 0.010200364, 1.00, 0.00307574, 0.072316905, 0.012224939,
  0.00770077, 0, 0.00307574, 1.00, 0.002325581, 0.000487329,
  0.036376355, 0.021136063, 0.072316905, 0.002325581, 1.00, 0.007371007,
  0.010240113, 0.002260398, 0.012224939, 0.000487329, 0.007371007, 1.00
), nrow = 6, byrow = TRUE)

#average across all datasets
jaccard_matrix <- matrix(c(
  1.00, 0.017223507, 0.127489823, 0.032519077, 0.034155683, 0.04131882,
  0.017223507, 1.00, 0.018674241, 0.010607774, 0.017472031, 0.011140216,
  0.127489823, 0.018674241, 1.00, 0.1956942, 0.040007657, 0.037459778,
  0.032519077, 0.010607774, 0.1956942, 1.00, 0.021357562, 0.006281453,
  0.034155683, 0.017472031, 0.040007657, 0.021357562, 1.00, 0.011370448,
  0.04131882, 0.011140216, 0.037459778, 0.006281453, 0.011370448, 1.00
), nrow = 6, byrow = TRUE)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("rMATS", "SUPPA2", "MAJIQ", "leafcutter", "DEXSeq", "whippet")

library(reshape2)
library(ggplot2)
library(dplyr)

# Melt the matrix
melted_jaccard <- melt(jaccard_matrix)

# Filter out the upper triangle
melted_jaccard <- melted_jaccard %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))

# Heatmap with custom gradient and theme
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/average_across_all_jaccard_heatmap_avg_all_datasets.png", width = 1100, height = 800, res = 150)
ggplot(melted_jaccard, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3)), size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "#08306b", name = "Jaccard Index") +
  theme_minimal(base_size = 14) +
  labs(title = "Jaccard Index Overlap Between AS Tools", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_fixed()  # Keep tiles square
dev.off()

#Nice Design
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/average_across_all_jaccard_heatmap_colours.png", width = 1100, height = 800, res = 150)
ggplot(melted_jaccard, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3), color = value > 0.5), size = 4) +
  scale_fill_viridis_c(option = "C", direction = -1, name = "Jaccard Index") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Pairwise Jaccard Index Between AS Tools",
    subtitle = "Overlap of Detected DS Genes",
    x = "Tool B",
    y = "Tool A"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed()
dev.off()











#repeat for DS events
# Create a symmetric matrix of Jaccard indices
jaccard_matrix <- matrix(c(
  1.00, 0, 0.003819036, 0, 0.000675447, 0,
  0, 1.00, 0, 0, 0, 0,
  0.003819036, 0, 1.00, 0, 0.0015057, 0,
  0, 0, 0, 1.00, 0, 0,
  0.000675447, 0, 0.0015057, 0, 1.00, 0,
  0, 0, 0, 0, 0, 1.00
), nrow = 6, byrow = TRUE)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("rMATS", "SUPPA2", "MAJIQ", "leafcutter", "DEXSeq", "whippet")

library(reshape2)
library(ggplot2)
library(dplyr)

# Melt the matrix
melted_jaccard <- melt(jaccard_matrix)

# Filter out the upper triangle
melted_jaccard <- melted_jaccard %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))

# Heatmap with custom gradient and theme
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/ASD_event_jaccard_heatmap.png", width = 1100, height = 800, res = 150)
ggplot(melted_jaccard, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3)), size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "#003300", name = "Jaccard Index") +
  theme_minimal(base_size = 14) +
  labs(title = "Jaccard Index Overlap Between AS Events", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_fixed()  # Keep tiles square
dev.off()

#Nicer design
png("~/Documents/Big_Data_Biology/Research_Project/All_datasets_graphs/event_jaccard_heatmap_avg_all_datasets.png", width = 1100, height = 800, res = 150)
ggplot(melted_jaccard, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(value, 3), color = value > 0.5), size = 4) +
  scale_fill_viridis_c(option = "C", direction = -1, name = "Jaccard Index") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Pairwise Jaccard Index Between AS Tools",
    subtitle = "Overlap of Detected DS Events",
    x = "Tool B",
    y = "Tool A"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed()
dev.off()
