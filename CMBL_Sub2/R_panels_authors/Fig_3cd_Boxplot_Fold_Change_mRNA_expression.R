rm(list = ls())
# setwd() to SBA/Figures removed on 2 Nov 2026: that folder holds a newer, different Fig_3cd xlsx and the script overwrote the archived npj SBA panels. Run from the folder that holds the April data copy.

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library("readxl")

dat1 <- read_excel("Fig_3cd_Boxplot_Fold_change_mRNA_expression.xlsx", sheet = 1)
attach(dat1)

my_comparison1 <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))

dat.Glp1r <- filter(dat1, Gene == "Glp1r")
dat.Ins1 <- filter(dat1, Gene == "Ins1")

#------------------------------Glp1r-------------------------------
dat1.Glp1r <- dat.Glp1r
shapiro.test(dat1.Glp1r$Fold_change) #W = 0.88173, p-value = 0.09224 #Normal distribution
Glp1r.labs <- c("Glp1r (Insulin synthesis and secretion)")
names(Glp1r.labs) <- c("Glp1r")

plot.Glp1r <- ggboxplot(dat1.Glp1r, x = "Treatment", y = "Fold_change", 
                        palette = "npg",
                       color = "Treatment", fill = "Treatment", alpha = 0.5, 
                       add = "jitter",
                       ylim = c(-0.1, 2)
) +
  labs(y = "Fold change in\nmRNA expression (AU)")+
  facet_wrap(vars(Gene), labeller = labeller(Gene = Glp1r.labs))+
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(aes(label = after_stat(paste0('"Anova, "*italic(p)*" = ', p.format, '"'))), parse = TRUE,
                     label.y= 2 , label.x.npc = "center", size=2.5, 
                     method= "anova", hjust = 0.5)+
  stat_compare_means(comparisons = my_comparison1, method = "t.test", 
                     tip.length = 0.1, step.increase = 0.5,
                     aes(label = after_stat(p.signif)),paired = FALSE, size = 2.5)+ 
  # Annotate mean ± SEM at the base
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Glp1r

plot.Glp1r$coordinates$clip <- "off"  # the omnibus label sits at the top of the y range; do not clip it (limits unchanged)
ggsave("Fig_3c_Boxplot_Fold_change_Glp1r_expression.png", plot = plot.Glp1r,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_3c_Boxplot_Fold_change_Glp1r_expression.svg", plot = plot.Glp1r,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)

#------------------------------Ins1-------------------------------
dat1.Ins1 <- dat.Ins1
shapiro.test(dat1.Ins1$Fold_change) #W = 0.90338, p-value = 0.1754 #Normal distribution
Ins1.labs <- c("Ins1 (Insulin synthesis)")
names(Ins1.labs) <- c("Ins1")

plot.Ins1 <- ggboxplot(dat1.Ins1, x = "Treatment", y = "Fold_change", 
                        palette = "npg",
                        color = "Treatment", fill = "Treatment", alpha = 0.5, 
                        add = "jitter",
                        ylim = c(-0.1, 2)
) +
  labs(y = "Fold change in\nmRNA expression (AU)")+
  facet_wrap(vars(Gene), labeller = labeller(Gene = Ins1.labs))+
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, colour = "black", face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(aes(label = after_stat(paste0('"Anova, "*italic(p)*" = ', p.format, '"'))), parse = TRUE,
                     label.y= 2 , label.x.npc = "center", size=2.5, 
                     method= "anova", hjust = 0.5)+
  stat_compare_means(comparisons = my_comparison1, method = "t.test", 
                     tip.length = 0.1, step.increase = 0.4,
                     aes(label = after_stat(p.signif)),paired = FALSE, size = 2.5)+ 
  # Annotate mean ± SEM at the base
  stat_summary(fun.data = function(x) {
    mean_value <- mean(x)
    sem_value <- sd(x) / sqrt(length(x))
    data.frame(y = 0, label = paste0(round(mean_value, 2), " ± ", round(sem_value, 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot.Ins1

plot.Ins1$coordinates$clip <- "off"  # the omnibus label sits at the top of the y range; do not clip it (limits unchanged)
ggsave("Fig_3d_Boxplot_Fold_change_Ins1_expression.png", plot = plot.Ins1,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_3d_Boxplot_Fold_change_Ins1_expression.svg", plot = plot.Ins1,
       scale = 1, width = 8, height =5.3, units = "cm",
       dpi = 600, limitsize = TRUE)
