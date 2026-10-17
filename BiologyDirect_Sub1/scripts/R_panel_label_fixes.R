# Corrected Fig. 1c panel (the only panel left in this script; the Fig. 3a section was removed, see the end of the file). Same data, same code and size as the
# revised original in SBA/Figures/github (Fig_1c_Boxplot_proliferation_revised.r); only the text labels change:
#   Fig. 1c: strip label "INS-1 proliferation" -> "Total DNA content"; axis title typo "concent" -> "content"
# The corrected PNG/SVG files must be dropped into the authors' composite (PowerPoint) by hand.
# Usage: Rscript R_panel_label_fixes.R <SBA/Figures dir> <out dir>
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr); library(ggsci); library(dplyr); library(readxl) })
args <- commandArgs(trailingOnly = TRUE); SRC <- args[1]; OUT <- args[2]

# ---- Fig. 1c (statistics and annotation exactly as in SBA/Figures/github/Fig_1c_Boxplot_proliferation_revised.r: Bonferroni-adjusted pairwise t-tests, "ns#" flag) ----
dat1 <- read_excel(file.path(SRC, "Fig_1c_Boxplot_INS1_proliferation.xlsx"), sheet = 1)
dat1$Assay <- "Total DNA content"
pw    <- pairwise.t.test(dat1$DNA, dat1$Treatment, p.adjust.method = "bonferroni")
pw_df <- as.data.frame(as.table(pw$p.value)); colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>% mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>% select(-temp)
pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***", ifelse(pw_df$p_adjusted < 0.01, "**", ifelse(pw_df$p_adjusted < 0.05, "*", "ns")))
pw_raw    <- pairwise.t.test(dat1$DNA, dat1$Treatment, p.adjust.method = "none")
pw_raw_df <- as.data.frame(as.table(pw_raw$p.value)); colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>% mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>% select(-temp)
pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")], by = c("Group1", "Group2"), all.x = TRUE)
pw_plot <- pw_df %>% rename(group1 = Group1, group2 = Group2) %>%
  mutate(p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05, "ns#", Significance), y.position = c(0.44, 0.52, 0.60))
plot1 <- ggboxplot(dat1, x = "Treatment", y = "DNA", palette = "npg", color = "Treatment", fill = "Treatment", alpha = 0.5, add = "jitter", ylim = c(0, 0.65)) +
  labs(y = "
Total DNA content (µg)") + facet_wrap(vars(Assay)) + theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(size = 8, colour = "black", face = "bold"), axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5), axis.text.y = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 0.64, label.x.npc = "center", size = 2.5, method = "anova", hjust = 0.5) +
  stat_pvalue_manual(pw_plot, label = "p.adj.signif", tip.length = 0.03, size = 2.5, step.increase = 0.001) +
  stat_summary(fun.data = function(x) data.frame(y = 0, label = paste0(round(mean(x), 2), " ± ", round(sd(x) / sqrt(length(x)), 2))), geom = "text", vjust = 0, size = 2.5)
ggsave(file.path(OUT, "Fig_1c_Boxplot_DNA_content.png"), plot1, width = 8.08, height = 4.57, units = "cm", dpi = 600)
ggsave(file.path(OUT, "Fig_1c_Boxplot_DNA_content.svg"), plot1, width = 8.08, height = 4.57, units = "cm")

# ---- Fig. 3a: REMOVED (15 Oct 2026) ------------------------------------------------------------------------------------------------
# The former Fig. 3a (LOESS curve with 95% bands of C-peptide in ng/µg DNA) treated the four ELISA readings of ONE culture well per group and
# glucose concentration as independent points. Fig. 3a and 3b are now drawn by R_Fig3ab_cpeptide_descriptive.R from the raw ELISA workbook
# (Cpeptide_raw_elisa_tidy.csv, made by elisa_raw_audit.py): raw signal, no smoothing, no bands, no tests.
cat("done\n")
