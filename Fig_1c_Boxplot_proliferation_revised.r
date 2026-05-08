# INS1 Cell Proliferation Analysis - Figure 1c

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)

# ── LOAD DATA ────────────────────────────────────────────────────────────────

dat1 <- read_excel("Fig_1c_Boxplot_INS1_proliferation.xlsx", sheet = 1)
attach(dat1)

# ── STATISTICS ───────────────────────────────────────────────────────────────

sw         <- shapiro.test(dat1$DNA)          # W = 0.8916, p = 0.1236 (normal)
aov_result <- aov(DNA ~ Treatment, data = dat1)
aov_sum    <- summary(aov_result)
aov_F      <- aov_sum[[1]][["F value"]][1]
aov_p      <- aov_sum[[1]][["Pr(>F)"]][1]

# Bonferroni-corrected pairwise comparisons
pw    <- pairwise.t.test(dat1$DNA, dat1$Treatment, p.adjust.method = "bonferroni")
pw_df <- as.data.frame(as.table(pw$p.value))
colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
  mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
  select(-temp)
pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                             ifelse(pw_df$p_adjusted < 0.01,  "**",
                                    ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))

# Raw (uncorrected) pairwise comparisons
pw_raw    <- pairwise.t.test(dat1$DNA, dat1$Treatment, p.adjust.method = "none")
pw_raw_df <- as.data.frame(as.table(pw_raw$p.value))
colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>%
  mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
  select(-temp)

# Merge and flag trends
pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")],
               by = c("Group1", "Group2"), all.x = TRUE)
pw_df$Trend <- ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted >= 0.05, "\u2020",
                      ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted < 0.05,
                             pw_df$Significance, "ns"))

desc_stats <- dat1 %>%
  group_by(Treatment) %>%
  summarise(Mean = mean(DNA), SEM = sd(DNA) / sqrt(n()),
            N = n(), .groups = "drop")

stat_summary_df <- data.frame(
  Test           = c("Shapiro-Wilk (all groups)", "One-way ANOVA"),
  Statistic      = c(round(sw$statistic, 4), round(aov_F, 4)),
  p_value        = c(round(sw$p.value, 4),   round(aov_p, 4)),
  Interpretation = c(
    ifelse(sw$p.value > 0.05, "Normal distribution", "Non-normal distribution"),
    ifelse(aov_p      < 0.05, "Significant",         "Not significant")
  )
)

# ── EXPORT CSV ───────────────────────────────────────────────────────────────

sec1 <- data.frame(
  Section      = c("Descriptive Statistics", rep("", nrow(desc_stats))),
  Treatment    = c("Treatment",  as.character(desc_stats$Treatment)),
  Mean         = c("Mean (µg)", as.character(round(desc_stats$Mean, 4))),
  SEM          = c("SEM",       as.character(round(desc_stats$SEM,  4))),
  N            = c("N",         as.character(desc_stats$N)),
  Significance = NA,
  Trend        = NA
)
sec2 <- data.frame(
  Section      = c("Normality & ANOVA", rep("", nrow(stat_summary_df))),
  Treatment    = c("Test",           as.character(stat_summary_df$Test)),
  Mean         = c("Statistic",      as.character(stat_summary_df$Statistic)),
  SEM          = c("p_value",        as.character(stat_summary_df$p_value)),
  N            = c("Interpretation", as.character(stat_summary_df$Interpretation)),
  Significance = NA,
  Trend        = NA
)
sec3 <- data.frame(
  Section      = c("Pairwise Comparisons (Bonferroni + Raw)", rep("", nrow(pw_df))),
  Treatment    = c("Group1",       as.character(pw_df$Group1)),
  Mean         = c("Group2",       as.character(pw_df$Group2)),
  SEM          = c("p_raw",        as.character(round(pw_df$p_raw,      6))),
  N            = c("p_adjusted",   as.character(round(pw_df$p_adjusted, 6))),
  Significance = c("Significance", as.character(pw_df$Significance)),
  Trend        = c("Trend_flag",   as.character(pw_df$Trend))
)
spacer   <- data.frame(Section = "---", Treatment = NA, Mean = NA,
                       SEM = NA, N = NA, Significance = NA, Trend = NA)
write.csv(rbind(sec1, spacer, sec2, spacer, sec3), "Fig_1c_stats.csv",
          row.names = FALSE)

# ── PLOT ANNOTATIONS ─────────────────────────────────────────────────────────

# Embed "ns#" for trend pairs — single bracket layer, matches all other figures
pw_plot <- pw_df %>%
  rename(group1 = Group1, group2 = Group2) %>%
  mutate(
    p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                          "ns#",
                          Significance),
    y.position   = c(0.44, 0.52, 0.60)
  )

# ── PLOT ─────────────────────────────────────────────────────────────────────

plot1 <- ggboxplot(dat1, x = "Treatment", y = "DNA",
                   palette = "npg", color = "Treatment", fill = "Treatment",
                   alpha = 0.5, add = "jitter", ylim = c(0, 0.65)) +
  labs(y = "\nTotal DNA concent (µg)") +
  facet_wrap(vars(Assay)) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text   = element_text(size = 8, colour = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 8, colour = "black",
                                    angle = 0, hjust = 0.5),
        axis.text.y  = element_text(size = 8, colour = "black",
                                    angle = 90, hjust = 0.5),
        axis.title.y = element_text(size = 8)) +
  stat_compare_means(label.y = 0.64, label.x.npc = "center", size = 2.5,
                     method = "anova", hjust = 0.5) +
  stat_pvalue_manual(pw_plot, label = "p.adj.signif",
                     tip.length = 0.03, size = 2.5, step.increase = 0.001) +
  stat_summary(fun.data = function(x) {
    data.frame(y = 0,
               label = paste0(round(mean(x), 2), " \u00b1 ",
                              round(sd(x) / sqrt(length(x)), 2)))
  }, geom = "text", vjust = 0, size = 2.5)

plot1

ggsave("Fig_1c_Boxplot_INS1_proliferation.png", plot = plot1,
       scale = 1, width = 8.08, height = 4.57, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_1c_Boxplot_INS1_proliferation.svg", plot = plot1,
       scale = 1, width = 8.08, height = 4.57, units = "cm",
       dpi = 600, limitsize = TRUE)

