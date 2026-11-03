# Clear workspace
rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

# Load libraries
library(ggplot2)
library(dplyr)
library(ggsci)
library(readxl)

dat2 <- read_excel("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.xlsx", sheet = 1)
attach(dat2)

dat3 <- filter(dat2, Marker == "Intracellular markers")
dat3 <- filter(dat3, Treatment != "CTRL")

# Calculate mean and SEM
dat_summary <- dat3 %>%
  group_by(Gene, Treatment) %>%
  summarise(
    mean_LogFC = mean(LogFC),
    sem_LogFC = sd(LogFC) / sqrt(n())
  )

# Get D3 palette colors
npg_colors <- pal_npg()(10)  # Pull first 10 colors from D3 palette

# Pick only the 2nd and 3rd colors
custom_colors <- c(npg_colors[2], npg_colors[3])

dat_summary$Gene <- factor(dat_summary$Gene, levels = c("Ptbp1", "Itpr1", 
                                                "Ryr2", "Ryr3"))

# Create bar plot
plot_bar <- ggplot(dat_summary, aes(x = Gene, y = mean_LogFC, fill = Treatment, color = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_LogFC - sem_LogFC, ymax = mean_LogFC + sem_LogFC),
                width = 0.2, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = custom_colors) +   # <--- updated here!
  scale_color_manual(values = custom_colors) +
  labs(x = "Gene", y = "Log2(Fold change)\nof mRNA expression") +
  lims(y = c(-1.5, 0.5))+
  theme_minimal() +
  theme(legend.position = "top",
        #strip.text = element_text(face = "italic"),
        legend.key.width  = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.spacing = margin(5),
        axis.title.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5, face = "italic"),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8))
plot_bar

# Save the plot
ggsave("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.png", plot = plot_bar,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_3e_Barplot_Log_Fold_Change_mRNA_expression.svg", plot = plot_bar,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

