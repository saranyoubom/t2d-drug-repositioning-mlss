rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(factoextra)
library(ggfortify)
library(ggrepel)

data <- read_excel("Fig_3g_PCA_Relative_mRNA_expression.xlsx", sheet = 1)
data <- data %>% select(-any_of("Ryr1"))
attach(data)

pca_data <- data %>% select(-Group)
x_order <- c("CTRL", "DAPT", "DKK-1")
data$Group <- factor(data$Group, levels = x_order)

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)
summary(pca_result)

# ── CSV Export ──────────────────────────────────────────────────────────────

# 1. PCA scores (individual coordinates)
scores <- as.data.frame(pca_result$x)
scores$Group <- data$Group
write.csv(scores, "Fig_3g_PCA_scores.csv", row.names = TRUE)

# 2. Loadings (variable coordinates)
loadings <- as.data.frame(pca_result$rotation)
write.csv(loadings, "Fig_3g_PCA_loadings.csv", row.names = TRUE)

# 3. Variance explained
var_explained <- summary(pca_result)$importance
write.csv(as.data.frame(var_explained), "Fig_3g_PCA_variance_explained.csv", row.names = TRUE)

# ── Biplot with labels pushed to arrow tips ─────────────────────────────────

# Extract variable coordinates for manual label placement
var_coords <- as.data.frame(pca_result$rotation[, 1:2])
colnames(var_coords) <- c("Dim1", "Dim2")

# Scale factor (same as factoextra default: sqrt of eigenvalue * n)
n <- nrow(pca_data)
eig <- pca_result$sdev^2
scale_factor <- sqrt(eig[1] * n) * 0.7   # 0.7 = factoextra default scaling

var_coords$Dim1_scaled <- var_coords$Dim1 * scale_factor
var_coords$Dim2_scaled <- var_coords$Dim2 * sqrt(eig[2] * n) * 0.7
var_coords$Gene <- rownames(var_coords)

# Push labels beyond arrow tip by a nudge factor
nudge <- 1.25   # increase to push labels further out
var_coords$label_x <- var_coords$Dim1_scaled * nudge
var_coords$label_y <- var_coords$Dim2_scaled * nudge

biplot <- fviz_pca_biplot(pca_result,
                          col.ind = data$Group,
                          palette = "npg",
                          addEllipses = TRUE,
                          ellipse.type = "confidence",
                          ellipse.level = 0.95,
                          label = "none",          # suppress default labels
                          col.var = "grey20",
                          repel = FALSE,
                          arrowsize = 0.2,
                          legend.title = "Treatment") +
  # Add gene labels manually at arrow tips with repel
  geom_text_repel(
    data = var_coords,
    aes(x = label_x, y = label_y, label = Gene),
    size = 2.5,
    color = "grey20",
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.size = 0.3,
    segment.color = "grey60",
    segment.alpha = 0.6,
    min.segment.length = 0.2
  ) +
  theme_minimal(base_size = 8) +
  lims(y = c(-5, 5), x = c(-5, 5)) +
  theme(legend.position = "top",
        title = element_blank(),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))

biplot

ggsave("Fig_3g_PCA_Relative_mRNA_expression.png", plot = biplot,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3g_PCA_Relative_mRNA_expression.svg", plot = biplot,
       scale = 1, width = 8, height = 5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

