# Functional Gene Expression - Figure 2c-f

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(outliers)

dat1 <- read_excel("Fig_2cdef_Fold_change_mRNA_expression.xlsx", sheet = 1)
attach(dat1)

my_comparison1 <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))

dat.Glut2   <- filter(dat1, Gene == "Glut2")
dat.Kcnj11  <- filter(dat1, Gene == "Kcnj11")
dat.Cacna1c <- filter(dat1, Gene == "Cacna1c")
dat.Cacna1d <- filter(dat1, Gene == "Cacna1d")

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: Grubbs' test per treatment group, remove outlier if significant
# ═══════════════════════════════════════════════════════════════════════════════

grubbs_clean <- function(df, value_col = "Fold_change", group_col = "Treatment",
                         gene_name = "") {
  groups   <- unique(df[[group_col]])
  removed  <- list()
  
  for (grp in groups) {
    idx  <- which(df[[group_col]] == grp)
    vals <- df[[value_col]][idx]
    
    if (length(unique(vals)) < 3 || length(vals) < 7) {
      cat(gene_name, "|", grp, "— skipped (n <7 or insufficient variance)\n")
      next
    }
    
    gt <- grubbs.test(vals)
    cat(gene_name, "|", grp, "— Grubbs G =", round(gt$statistic["G"], 4),
        ", p =", round(gt$p.value, 4), "\n")
    
    if (gt$p.value < 0.05) {
      outlier_val <- ifelse(grepl("highest", gt$alternative), max(vals), min(vals))
      cat("  → Outlier removed:", outlier_val, "\n")
      removed[[grp]] <- outlier_val
      # Remove only first instance of the outlier value in that group
      first_hit <- idx[which(vals == outlier_val)[1]]
      df <- df[-first_hit, ]
    }
  }
  return(df)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: compute Bonferroni + raw stats, export CSV for one gene
# ═══════════════════════════════════════════════════════════════════════════════

bonferroni_stats <- function(vals, groups, is_normal, csv_name) {
  sw <- shapiro.test(vals)
  
  if (is_normal) {
    aov_res   <- summary(aov(vals ~ groups))
    omni_name <- "One-way ANOVA"
    omni_stat <- round(aov_res[[1]][["F value"]][1], 4)
    omni_p    <- round(aov_res[[1]][["Pr(>F)"]][1], 4)
    pw        <- pairwise.t.test(vals, groups, p.adjust.method = "bonferroni")
    pw_raw    <- pairwise.t.test(vals, groups, p.adjust.method = "none")
  } else {
    kw        <- kruskal.test(vals ~ groups)
    omni_name <- "Kruskal-Wallis"
    omni_stat <- round(kw$statistic, 4)
    omni_p    <- round(kw$p.value, 4)
    pw        <- pairwise.wilcox.test(vals, groups, p.adjust.method = "bonferroni")
    pw_raw    <- pairwise.wilcox.test(vals, groups, p.adjust.method = "none")
  }
  
  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
  pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                               ifelse(pw_df$p_adjusted < 0.01,  "**",
                                      ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))
  
  pw_raw_df <- as.data.frame(as.table(pw_raw$p.value))
  colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
  pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  
  pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")],
                 by = c("Group1", "Group2"), all.x = TRUE)
  pw_df$Trend <- ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted >= 0.05, "†",
                        ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted < 0.05,
                               pw_df$Significance, "ns"))
  
  desc <- data.frame(Treatment = levels(factor(groups))) %>%
    rowwise() %>%
    mutate(Mean = mean(vals[groups == Treatment]),
           SEM  = sd(vals[groups == Treatment]) / sqrt(sum(groups == Treatment)),
           N    = sum(groups == Treatment)) %>%
    ungroup()
  
  stat_df <- data.frame(
    Test      = c("Shapiro-Wilk", omni_name),
    Statistic = c(round(sw$statistic, 4), omni_stat),
    p_value   = c(round(sw$p.value,   4), omni_p),
    Interpretation = c(
      ifelse(sw$p.value > 0.05, "Normal", "Non-normal"),
      ifelse(omni_p     < 0.05, "Significant", "Not significant")
    )
  )
  
  sec1 <- data.frame(
    Section      = c("Descriptive Statistics", rep("", nrow(desc))),
    Treatment    = c("Treatment", as.character(desc$Treatment)),
    Mean         = c("Mean (AU)", as.character(round(desc$Mean, 4))),
    SEM          = c("SEM",       as.character(round(desc$SEM,  4))),
    N            = c("N",         as.character(desc$N)),
    Significance = NA, Trend = NA
  )
  sec2 <- data.frame(
    Section      = c("Normality & Omnibus Test", rep("", nrow(stat_df))),
    Treatment    = c("Test",           as.character(stat_df$Test)),
    Mean         = c("Statistic",      as.character(stat_df$Statistic)),
    SEM          = c("p_value",        as.character(stat_df$p_value)),
    N            = c("Interpretation", as.character(stat_df$Interpretation)),
    Significance = NA, Trend = NA
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
  spacer <- data.frame(Section = "---", Treatment = NA, Mean = NA, SEM = NA,
                       N = NA, Significance = NA, Trend = NA)
  write.csv(rbind(sec1, spacer, sec2, spacer, sec3), csv_name, row.names = FALSE)
  
  return(pw_df)
}

# ── Shared plot builder ───────────────────────────────────────────────────────

make_boxplot_gene <- function(dat, pw_df, gene_lab_vec, omni_method,
                              y_positions, png_name, svg_name) {
  
  # Embed "ns#" for trend pairs — single bracket layer, no alignment issue
  pw_plot <- pw_df %>%
    rename(group1 = Group1, group2 = Group2) %>%
    mutate(
      p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                            "ns#",
                            Significance),
      y.position   = y_positions
    )
  
  p <- ggboxplot(dat, x = "Treatment", y = "Fold_change",
                 palette = "npg",
                 color = "Treatment", fill = "Treatment", alpha = 0.5,
                 add = "jitter",
                 ylim = c(-0.1, 2.3)) +
    labs(y = "Fold change in\nmRNA expression (AU)") +
    facet_wrap(vars(Gene), labeller = labeller(Gene = gene_lab_vec)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text   = element_text(size = 8, colour = "black", face = "bold.italic"),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 8, colour = "black", angle = 0, hjust = 0.5),
      axis.text.y  = element_text(size = 8, colour = "black", angle = 90, hjust = 0.5),
      axis.title.y = element_text(size = 8)
    ) +
    stat_compare_means(label.y = 2.3, label.x.npc = "center", size = 2.5,
                       method = omni_method, hjust = 0.5) +
    stat_pvalue_manual(pw_plot, label = "p.adj.signif",
                       tip.length = 0.03, size = 2.5, step.increase = 0.08) +
    stat_summary(fun.data = function(x) {
      data.frame(y = 0,
                 label = paste0(round(mean(x), 2), " \u00b1 ",
                                round(sd(x) / sqrt(length(x)), 2)))
    }, geom = "text", vjust = 1, size = 2.5)
  
  ggsave(png_name, plot = p, scale = 1, width = 8, height = 5.3,
         units = "cm", dpi = 600, limitsize = TRUE)
  ggsave(svg_name, plot = p, scale = 1, width = 8, height = 5.3,
         units = "cm", dpi = 600, limitsize = TRUE)
  
  return(p)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Glut2 (Glucose uptake and sensing) — normal → ANOVA + t-test
# ═══════════════════════════════════════════════════════════════════════════════

dat1.Glut2 <- dat.Glut2
dat1.Glut2$Fold_change <- as.numeric(dat1.Glut2$Fold_change)

cat("\n── Glut2: Grubbs outlier screen ──\n")
dat1.Glut2 <- grubbs_clean(dat1.Glut2, gene_name = "Glut2")
shapiro.test(dat1.Glut2$Fold_change)

pw_glut2 <- bonferroni_stats(dat1.Glut2$Fold_change, dat1.Glut2$Treatment,
                             is_normal = TRUE, csv_name = "Fig_2c_stats.csv")

Glut2.labs <- c("Glut2" = "Glut2 (Glucose uptake and sensing)")

plot.Glut2 <- make_boxplot_gene(
  dat          = dat1.Glut2,
  pw_df        = pw_glut2,
  gene_lab_vec = Glut2.labs,
  omni_method  = "anova",
  y_positions  = c(1.55, 1.70, 1.85),
  png_name     = "Fig_2c_Boxplot_Fold_change_Glut2_expression.png",
  svg_name     = "Fig_2c_Boxplot_Fold_change_Glut2_expression.svg"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Kcnj11 (ATP-sensitive potassium channel) — normal → ANOVA + t-test
# ═══════════════════════════════════════════════════════════════════════════════

dat1.Kcnj11 <- dat.Kcnj11
dat1.Kcnj11$Fold_change <- as.numeric(dat1.Kcnj11$Fold_change)

cat("\n── Kcnj11: Grubbs outlier screen ──\n")
dat1.Kcnj11 <- grubbs_clean(dat1.Kcnj11, gene_name = "Kcnj11")
shapiro.test(dat1.Kcnj11$Fold_change)

pw_kcnj11 <- bonferroni_stats(dat1.Kcnj11$Fold_change, dat1.Kcnj11$Treatment,
                              is_normal = TRUE, csv_name = "Fig_2d_stats.csv")

Kcnj11.labs <- c("Kcnj11" = "Kcnj11 (ATP-sensitive potassium channel)")

plot.Kcnj11 <- make_boxplot_gene(
  dat          = dat1.Kcnj11,
  pw_df        = pw_kcnj11,
  gene_lab_vec = Kcnj11.labs,
  omni_method  = "anova",
  y_positions  = c(1.55, 1.70, 1.85),
  png_name     = "Fig_2d_Boxplot_Fold_change_Kcnj11_expression.png",
  svg_name     = "Fig_2d_Boxplot_Fold_change_Kcnj11_expression.svg"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Cacna1c (Calcium voltage-gated channel) — normal → ANOVA + t-test
# ═══════════════════════════════════════════════════════════════════════════════

dat1.Cacna1c <- dat.Cacna1c
dat1.Cacna1c$Fold_change <- as.numeric(dat1.Cacna1c$Fold_change)

cat("\n── Cacna1c: Grubbs outlier screen ──\n")
dat1.Cacna1c <- grubbs_clean(dat1.Cacna1c, gene_name = "Cacna1c")
shapiro.test(dat1.Cacna1c$Fold_change)

pw_cacna1c <- bonferroni_stats(dat1.Cacna1c$Fold_change, dat1.Cacna1c$Treatment,
                               is_normal = TRUE, csv_name = "Fig_2e_stats.csv")

Cacna1c.labs <- c("Cacna1c" = "Cacna1c (Calcium voltage-gated channel)")

plot.Cacna1c <- make_boxplot_gene(
  dat          = dat1.Cacna1c,
  pw_df        = pw_cacna1c,
  gene_lab_vec = Cacna1c.labs,
  omni_method  = "anova",
  y_positions  = c(1.55, 1.70, 1.85),
  png_name     = "Fig_2e_Boxplot_Fold_change_Cacna1c_expression.png",
  svg_name     = "Fig_2e_Boxplot_Fold_change_Cacna1c_expression.svg"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Cacna1d (Calcium voltage-gated channel) — normal → ANOVA + t-test
# ═══════════════════════════════════════════════════════════════════════════════

dat1.Cacna1d <- dat.Cacna1d
dat1.Cacna1d$Fold_change <- as.numeric(dat1.Cacna1d$Fold_change)

cat("\n── Cacna1d: Grubbs outlier screen ──\n")
dat1.Cacna1d <- grubbs_clean(dat1.Cacna1d, gene_name = "Cacna1d")
shapiro.test(dat1.Cacna1d$Fold_change)

pw_cacna1d <- bonferroni_stats(dat1.Cacna1d$Fold_change, dat1.Cacna1d$Treatment,
                               is_normal = TRUE, csv_name = "Fig_2f_stats.csv")

Cacna1d.labs <- c("Cacna1d" = "Cacna1d (Calcium voltage-gated channel)")

plot.Cacna1d <- make_boxplot_gene(
  dat          = dat1.Cacna1d,
  pw_df        = pw_cacna1d,
  gene_lab_vec = Cacna1d.labs,
  omni_method  = "anova",
  y_positions  = c(1.55, 1.70, 1.85),
  png_name     = "Fig_2f_Boxplot_Fold_change_Cacna1d_expression.png",
  svg_name     = "Fig_2f_Boxplot_Fold_change_Cacna1d_expression.svg"
)

cat("\nAll Fig 2c-f files saved.\n")

