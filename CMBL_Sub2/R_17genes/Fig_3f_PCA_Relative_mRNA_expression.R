# Clear workspace
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

# Load required libraries
library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(factoextra)
library(ggfortify)

# Load the dataset
data <- read_excel("Fig_3f_PCA_Relative_mRNA_expression.xlsx", sheet = 1)
attach(data)

# Remove non-numeric columns for PCA
pca_data <- data %>% select(-Group)

x_order <- c("CTRL", "DAPT", "DKK-1")
data$Group <- factor(data$Group, levels = x_order)

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Summarize PCA results
summary(pca_result)

# Scree Plot
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# Loading Plot
loadings <- as.data.frame(pca_result$rotation)
loadings$Variables <- rownames(loadings)

ggplot(loadings, aes(x = PC1, y = PC2, label = Variables)) +
  geom_point(size = 4, color = "blue") +
  geom_text(vjust = -0.5, size = 4) +
  labs(title = "Loading Plot", x = "PC1", y = "PC2") +
  theme_minimal()

# Individual Samples Plot
fviz_pca_ind(pca_result,
             #geom.ind = "point", 
             col.ind = data$Group, 
             palette = "npg", 
             addEllipses = TRUE, 
             ellipse.level = 0.95, 
             legend.title = "Treatment")

biplot <- fviz_pca_biplot(pca_result, 
                col.ind = data$Group, palette = "npg", 
                addEllipses = TRUE, label = "var",
                col.var = "grey20", repel = TRUE, 
                arrowsize = 0.4,
                legend.title = "Treatment") +
  theme_minimal(base_size = 8)+
  lims (y = c(-5, 4), x = c(-6, 8))+
  theme(legend.position = "top", title = element_blank(),
        legend.key.width  = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = margin(5),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, vjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))


# Gene (variable) labels in italic; the biplot has no other text layers
for (i in seq_along(biplot$layers)) {
  if (inherits(biplot$layers[[i]]$geom, c("GeomText", "GeomTextRepel", "GeomLabel")))
    biplot$layers[[i]]$aes_params$fontface <- "italic"
}

# Display the biplot
biplot

# Save the PCA biplot
ggsave("Fig_3f_PCA_Relative_mRNA_expression.png", plot = biplot,
       scale = 1, width = 7.6, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_3f_PCA_Relative_mRNA_expression.svg", plot = biplot,
       scale = 1, width = 7.6, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
