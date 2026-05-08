# ════════════════════════════════════════════════════════════════════════════════
#
#  COMBINED ANALYSIS SCRIPTS
#  Oontawee et al. — Notch-Wnt Pathway Inhibition in β-Cells (npj SBA, 2026)
#  Generated: 2026-04-23 18:58:41.470861
#  Total scripts: 21
#
#  SCRIPT ORDER:
#     1. Fig_1c_Boxplot_proliferation_revised.r — Fig 1c: Total DNA content — proliferation boxplot (ANOVA/KW + Bonferroni)
#     2. Fig_1d_Boxplot_Notch_Wnt_genes.r — Fig 1d: Notch/Wnt gene expression boxplots — Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Tcf7, Lef1, Tcf7l2
#     3. Fig_2b_glucose_uptake_revised.r — Fig 2b: 2-NBDG glucose uptake fluorescence boxplot
#     4. Fig_2cdef_functional_genes_revised.r — Fig 2c–2f: Glut2, Kcnj11, Cacna1c, Cacna1d fold change boxplots
#     5. Fig_3a_C-peptide_curve.r — Fig 3a: Glucose-stimulated C-peptide secretion LOESS curve
#     6. Fig_3b_C-peptide_5.5mm_revised.r — Fig 3b: C-peptide secretion at 5.5 mM glucose — boxplot
#     7. Fig_3ef_Boxplot_script.r — Fig 3c–3f: Glp1r, Ins1, Ptbp1, Itpr1, Ryr1–3 gene expression boxplots
#     8. Fig_3g_PCA_script.r — Fig 3g: PCA of 11 secretory-pathway genes (biplot)
#     9. Fig_4b_correlation_matrix.r — Fig 4b: Pearson correlation matrix heatmap (Hes1-normalized)
#    10. Fig_4c_normalization_methods.r — Fig 4c: Multi-method normalization comparison — 4 strategies
#    11. Fig_4d_Venn_diagram.r — Fig 4d: Four-way Venn diagram — 28-pair consensus network
#    12. Fig_4e_Hes1-normalized_correlation_matrix.R — Fig 4e: Hes1-normalized annotated correlation heatmap (consensus pairs marked)
#    13. Fig_S1_validation_analysis.r — Fig S1: Validation analysis — observed vs null distribution, false positive rates, sign-flip check
#    14. Fig_5abcd_network_analysis.r — Fig 5a–5d: Gene interaction network — 15 nodes, STRING classification, novel pair ranking
#    15. Fig_6b_consensus_hub_network.r — Fig 6b: Consensus hub gene network (circular layout, igraph/ggraph)
#    16. Fig_6c_MLSS_drug_combinations.r — Fig 6c: MLSS v4.0 scoring — all drug combinations ranked by synergy score
#    17. Fig_6d_Fig_S2_MLSS_validation.r — Fig 6d + Fig S2: MLSS robustness — weight sensitivity, ablation analysis (8 scenarios)
#    18. Fig_6e_drug_combination_network.r — Fig 6e: Drug combination network — top 15 combinations (circular layout)
#    19. Fig_6f_regulatory_mechanistic_sankey.r — Fig 6f: Six-layer regulatory-mechanistic Sankey diagram
#    20. 0_Analysis Script Combination.R
#    21. 0_Directory File Checker.R
#
# ════════════════════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 01/21  —  Fig_1c_Boxplot_proliferation_revised.r                      ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 1c: Total DNA content — proliferation boxplot (ANOVA/KW + Bonferroni)
# Lines: 146  |  File: Fig_1c_Boxplot_proliferation_revised.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

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

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 01: Fig_1c_Boxplot_proliferation_revised.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 02/21  —  Fig_1d_Boxplot_Notch_Wnt_genes.r                            ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 1d: Notch/Wnt gene expression boxplots — Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Tcf7, Lef1, Tcf7l2
# Lines: 374  |  File: Fig_1d_Boxplot_Notch_Wnt_genes.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Figure 1d — Boxplots: Notch targets & Wnt pathway gene expression
# Theme: matched exactly to Fig_2cdef | Layout: 5 columns × 2 rows | 16.2 × 12.8 cm

rm(list = ls())

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(rstatix)
library(patchwork)
library(rlang)

npg_colors <- pal_npg("nrc")(10)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
dat <- read_excel("Fig_1d_Fold_change_mRNA_expression.xlsx", sheet = 1)

if (!"Gene" %in% colnames(dat)) {
  dat <- dat %>%
    pivot_longer(cols = -c(Treatment, Sample),
                 names_to  = "Gene",
                 values_to = "Fold_change")
}

gene_order   <- c("Hes1", "Hey1", "Wnt2", "Wnt2b", "Wnt5a",
                  "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2")
treat_levels <- c("CTRL", "DAPT", "DKK-1")

dat <- dat %>%
  filter(Gene %in% gene_order) %>%
  mutate(Gene        = factor(Gene,      levels = gene_order),
         Treatment   = factor(Treatment, levels = treat_levels),
         Fold_change = as.numeric(Fold_change))

# ── HELPER: STATISTICS + CSV EXPORT ──────────────────────────────────────────
run_stats <- function(df, gene_name, parametric = TRUE) {
  
  x   <- df$Fold_change
  grp <- df$Treatment
  
  sw <- tryCatch(shapiro.test(x),
                 error   = function(e) list(statistic = NA, p.value = NA),
                 warning = function(w) list(statistic = NA, p.value = NA))
  
  if (parametric) {
    fit           <- lm(Fold_change ~ Treatment, data = df)
    fit_anv       <- anova(fit)
    omni_stat     <- round(as.numeric(fit_anv[["F value"]][1]), 4)
    omni_p        <- round(as.numeric(fit_anv[["Pr(>F)"]][1]), 4)
    omni_name     <- "One-way ANOVA (F)"
    pw            <- pairwise.t.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.t.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise t-test (Bonferroni)"
  } else {
    omni          <- kruskal.test(Fold_change ~ Treatment, data = df)
    omni_stat     <- round(as.numeric(omni$statistic), 4)
    omni_p        <- round(as.numeric(omni$p.value), 4)
    omni_name     <- "Kruskal-Wallis (H)"
    pw            <- pairwise.wilcox.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.wilcox.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise Wilcoxon (Bonferroni)"
  }
  
  # Bonferroni-corrected p-values
  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
  pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                               ifelse(pw_df$p_adjusted < 0.01,  "**",
                                      ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))
  
  # Raw (uncorrected) p-values
  pw_raw_df <- as.data.frame(as.table(pw_raw$p.value))
  colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
  pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  
  # Merge raw p into pw_df and flag trends
  pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")],
                 by = c("Group1", "Group2"), all.x = TRUE)
  pw_df$Trend <- ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted >= 0.05, "#",
                        ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted < 0.05,
                               pw_df$Significance, "ns"))
  
  desc <- df %>%
    group_by(Treatment) %>%
    summarise(Mean = mean(Fold_change), SEM = sd(Fold_change) / sqrt(n()),
              N = n(), .groups = "drop")
  
  sec1 <- data.frame(
    Section = c(paste0("Descriptive Statistics - ", gene_name), rep("", nrow(desc))),
    Col1    = c("Treatment",  as.character(desc$Treatment)),
    Col2    = c("Mean (AU)",  as.character(round(desc$Mean, 4))),
    Col3    = c("SEM",        as.character(round(desc$SEM,  4))),
    Col4    = c("N",          as.character(desc$N)),
    Col5    = NA,
    Col6    = NA
  )
  
  stat_df <- data.frame(
    Test = c(paste0("Shapiro-Wilk (", gene_name, ")"), omni_name),
    Statistic = c(round(as.numeric(sw$statistic), 4), omni_stat),
    p_value   = c(round(as.numeric(sw$p.value),   4), omni_p),
    Interpretation = c(
      ifelse(is.na(sw$p.value), "Undetermined",
             ifelse(sw$p.value > 0.05, "Normal", "Non-normal")),
      ifelse(is.na(omni_p), "Undetermined",
             ifelse(omni_p < 0.05, "Significant", "Not significant"))
    )
  )
  
  sec2 <- data.frame(
    Section = c(paste0("Normality & Omnibus - ", gene_name), rep("", nrow(stat_df))),
    Col1    = c("Test",           as.character(stat_df$Test)),
    Col2    = c("Statistic",      as.character(stat_df$Statistic)),
    Col3    = c("p_value",        as.character(stat_df$p_value)),
    Col4    = c("Interpretation", as.character(stat_df$Interpretation)),
    Col5    = NA,
    Col6    = NA
  )
  
  sec3 <- data.frame(
    Section = c(paste0("Pairwise - ", gene_name, " (", post_hoc_name, ")"),
                rep("", nrow(pw_df))),
    Col1    = c("Group1",       as.character(pw_df$Group1)),
    Col2    = c("Group2",       as.character(pw_df$Group2)),
    Col3    = c("p_raw",        as.character(round(pw_df$p_raw,      6))),
    Col4    = c("p_adjusted",   as.character(round(pw_df$p_adjusted, 6))),
    Col5    = c("Significance", as.character(pw_df$Significance)),
    Col6    = c("Trend_flag",   as.character(pw_df$Trend))
  )
  
  spacer <- data.frame(Section = "---", Col1 = NA, Col2 = NA, Col3 = NA,
                       Col4 = NA, Col5 = NA, Col6 = NA)
  out    <- rbind(sec1, spacer, sec2, spacer, sec3, spacer)
  
  list(csv = out, pw = pw_df)
}

# ── NORMALITY SCREEN (non-CTRL only) ─────────────────────────────────────────
safe_shapiro_p <- function(x) {
  x <- x[!is.na(x)]
  if (length(unique(x)) < 3 || length(x) < 3) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value,
           error   = function(e) NA_real_,
           warning = function(w) NA_real_)
}

parametric_map <- sapply(gene_order, function(g) {
  pvals <- dat %>%
    filter(Gene == g, Treatment != "CTRL") %>%
    group_by(Treatment) %>%
    summarise(p = safe_shapiro_p(Fold_change), .groups = "drop") %>%
    pull(p)
  all(pvals > 0.05, na.rm = TRUE)
})

cat("Normality screen (non-CTRL groups):\n")
print(data.frame(
  Gene       = gene_order,
  Parametric = unname(parametric_map),
  Test       = ifelse(unname(parametric_map),
                      "One-way ANOVA + Bonferroni",
                      "Kruskal-Wallis + Wilcoxon Bonferroni")
))

# ── RUN STATISTICS PER GENE ──────────────────────────────────────────────────
gene_dfs <- lapply(gene_order, function(g) dat %>% filter(Gene == g))
names(gene_dfs) <- gene_order

res <- lapply(gene_order, function(g) {
  cat("Processing:", g, "| parametric =", unname(parametric_map[g]), "\n")
  run_stats(gene_dfs[[g]], g, unname(parametric_map[g]))
})
names(res) <- gene_order

combined_csv <- do.call(rbind, lapply(res, `[[`, "csv"))
write.csv(combined_csv, "Fig_1d_stats.csv", row.names = FALSE)

# ── HELPER: build pw_plot with ns# embedded ───────────────────────────────────
make_pw_plot <- function(pw_df, y_positions) {
  pw_df %>%
    rename(group1 = Group1, group2 = Group2) %>%
    mutate(
      p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                            "ns#",
                            Significance),
      y.position   = y_positions
    )
}

# ── HELPER: boxplot builder ───────────────────────────────────────────────────
make_boxplot <- function(df, pw_plot, gene_lab, gene_key,
                         omnibus_method, y_lab_pos, y_lim_top,
                         show_y_title = TRUE,
                         show_y_text  = TRUE,
                         summary_y    = 0.25) {
  labs_vec <- setNames(gene_lab, gene_key)
  
  p <- ggboxplot(df,
                 x       = "Treatment",
                 y       = "Fold_change",
                 palette = "npg",
                 color   = "Treatment",
                 fill    = "Treatment",
                 alpha   = 0.5,
                 add     = "jitter",
                 ylim    = c(-0.1, y_lim_top)) +
    labs(y = if (show_y_title) "Fold change in\nmRNA expression (AU)" else NULL) +
    facet_wrap(vars(Gene), labeller = labeller(Gene = labs_vec)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 8, colour = "black", face = "bold.italic"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = 8, colour = "black",
                                     angle = 45, hjust = 1, vjust = 1),
      axis.text.y     = if (show_y_text)
        element_text(size = 8, colour = "black", angle = 90, hjust = 0.5)
      else
        element_blank(),
      axis.title.y    = if (show_y_title) element_text(size = 8) else element_blank(),
      plot.margin     = margin(0, 0, 0, 0)
    ) +
    stat_compare_means(
      label.y     = y_lab_pos,
      label.x.npc = "center",
      size        = 2.5,
      method      = omnibus_method,
      hjust       = 0.5
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          y     = summary_y,
          label = paste0(round(mean(x), 2), "\n\u00b1 ",
                         round(sd(x) / sqrt(length(x)), 2))
        )
      },
      geom  = "text",
      vjust = 1,
      size  = 2.5,
      angle = 0
    )
  
  if (nrow(pw_plot) > 0) {
    p <- p + stat_pvalue_manual(pw_plot,
                                label      = "p.adj.signif",
                                tip.length = 0.03,
                                size       = 2.5,
                                step.increase = 0.001)
  }
  
  return(p)
}

save_plot <- function(plot, filename) {
  ggsave(paste0(filename, ".png"), plot = plot, scale = 1,
         width = 8, height = 5.3, units = "cm", dpi = 600, limitsize = TRUE)
  ggsave(paste0(filename, ".svg"), plot = plot, scale = 1,
         width = 8, height = 5.3, units = "cm", dpi = 600, limitsize = TRUE)
}

# ── GENE LABELS ──────────────────────────────────────────────────────────────
gene_labels <- c(
  "Hes1"   = "Hes1 (Notch target)",
  "Hey1"   = "Hey1 (Notch target)",
  "Wnt2"   = "Wnt2 (Wnt ligand)",
  "Wnt2b"  = "Wnt2b (Wnt ligand)",
  "Wnt5a"  = "Wnt5a (NC-Wnt)",
  "Wnt5b"  = "Wnt5b (NC-Wnt)",
  "Wnt9a"  = "Wnt9a (Wnt ligand)",
  "Lef1"   = "Lef1 (Wnt TF)",
  "Tcf7"   = "Tcf7 (Wnt TF)",
  "Tcf7l2" = "Tcf7l2 (Wnt TF)"
)

y_title_genes <- c("Hes1", "Wnt5b")
y_text_genes  <- c("Hes1", "Wnt5b")

# ── PASS 1: per-gene y_lim_top candidates ────────────────────────────────────
gene_y_info <- lapply(gene_order, function(g) {
  df      <- gene_dfs[[g]]
  pw_df   <- res[[g]]$pw
  n_pairs <- nrow(pw_df)
  y_max   <- max(df$Fold_change, na.rm = TRUE)
  
  if (n_pairs > 0) {
    y_pos     <- seq(y_max * 1.12, y_max * 1.50, length.out = n_pairs)
    y_lim_top <- max(y_pos) * 1.15
  } else {
    y_pos     <- numeric(0)
    y_lim_top <- y_max * 1.65
  }
  list(y_lim_top = y_lim_top)
})
names(gene_y_info) <- gene_order

# ── PASS 2: shared ceiling & bracket positions ────────────────────────────────
global_y_lim_top <- max(sapply(gene_y_info, `[[`, "y_lim_top"))
max_n_pairs      <- max(sapply(gene_order, function(g) nrow(res[[g]]$pw)))

global_y_pos_all <- seq(global_y_lim_top * 0.62,
                        global_y_lim_top * 0.87,
                        length.out = max(max_n_pairs, 1))
global_y_lab_pos <- global_y_lim_top * 0.97

cat("Global y_lim_top:", round(global_y_lim_top, 3),
    "| bracket range:", round(global_y_pos_all[1], 3), "-",
    round(global_y_pos_all[length(global_y_pos_all)], 3),
    "| label pos:", round(global_y_lab_pos, 3), "\n")

# ── BUILD PLOTS ──────────────────────────────────────────────────────────────
plot_list <- list()

for (g in gene_order) {
  df       <- gene_dfs[[g]]
  pw_df    <- res[[g]]$pw
  is_param <- unname(parametric_map[g])
  omni_m   <- ifelse(is_param, "anova", "kruskal.test")
  n_pairs  <- nrow(pw_df)
  y_pos_g  <- if (n_pairs > 0) global_y_pos_all[seq_len(n_pairs)] else numeric(0)
  
  pw_plot <- make_pw_plot(pw_df, y_pos_g)
  
  plot_list[[g]] <- make_boxplot(
    df             = df,
    pw_plot        = pw_plot,
    gene_lab       = gene_labels[g],
    gene_key       = g,
    omnibus_method = omni_m,
    y_lab_pos      = global_y_lab_pos,
    y_lim_top      = global_y_lim_top,
    show_y_title   = g %in% y_title_genes,
    show_y_text    = g %in% y_text_genes,
    summary_y      = global_y_lim_top * 0.07
  )
  
  cat("Plot built:", g,
      "| n_brackets =", n_pairs,
      "| bracket y =", paste(round(y_pos_g, 2), collapse = ", "), "\n")
}

# ── SAVE INDIVIDUAL PANELS ───────────────────────────────────────────────────
panel_files <- c(
  "Hes1"   = "Fig_1d_1_Boxplot_Hes1",
  "Hey1"   = "Fig_1d_2_Boxplot_Hey1",
  "Wnt2"   = "Fig_1d_3_Boxplot_Wnt2",
  "Wnt2b"  = "Fig_1d_4_Boxplot_Wnt2b",
  "Wnt5a"  = "Fig_1d_5_Boxplot_Wnt5a",
  "Wnt5b"  = "Fig_1d_6_Boxplot_Wnt5b",
  "Wnt9a"  = "Fig_1d_7_Boxplot_Wnt9a",
  "Lef1"   = "Fig_1d_8_Boxplot_Lef1",
  "Tcf7"   = "Fig_1d_9_Boxplot_Tcf7",
  "Tcf7l2" = "Fig_1d_10_Boxplot_Tcf7l2"
)

for (g in gene_order) save_plot(plot_list[[g]], panel_files[g])

# ── ASSEMBLE 5×2 COMBINED FIGURE ─────────────────────────────────────────────
combined <- wrap_plots(plot_list, ncol = 5, nrow = 2)

ggsave("Fig_1d_Boxplot_combined_expression.png", plot = combined, scale = 1,
       width = 16.6, height = 12.8, units = "cm", dpi = 600, limitsize = TRUE)
ggsave("Fig_1d_Boxplot_combined_expression.svg", plot = combined, scale = 1,
       width = 16.6, height = 12.8, units = "cm", dpi = 600, limitsize = TRUE)

cat("\nAll Fig 1d files saved.\nStats: Fig_1d_stats.csv\n")

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 02: Fig_1d_Boxplot_Notch_Wnt_genes.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 03/21  —  Fig_2b_glucose_uptake_revised.r                             ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 2b: 2-NBDG glucose uptake fluorescence boxplot
# Lines: 165  |  File: Fig_2b_glucose_uptake_revised.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(outliers)

# ── LOAD DATA ────────────────────────────────────────────────────────────────

dat1 <- read_excel("Fig_2b_Boxplot_Glucose_uptake.xlsx", sheet = 1)

# ── OUTLIER DETECTION (Grubbs' test on CTRL) ─────────────────────────────────

ctrl_vals   <- dat1$Relative_intensity[dat1$Treatment == "CTRL"]
grubbs_ctrl <- grubbs.test(ctrl_vals)
print(grubbs_ctrl)  # G = 1.4987, p = 0.0017 → highest value is outlier

ctrl_outlier_val <- max(ctrl_vals)
dat1_clean <- dat1 %>%
  filter(!(Treatment == "CTRL" & Relative_intensity == ctrl_outlier_val))

# ── NORMALIZE TO MEAN OF CLEANED CTRL ────────────────────────────────────────

ctrl_mean  <- mean(dat1_clean$Relative_intensity[dat1_clean$Treatment == "CTRL"])
dat1_clean <- dat1_clean %>%
  mutate(Fold_change = Relative_intensity / ctrl_mean)

# ── STATISTICS ───────────────────────────────────────────────────────────────

sw        <- shapiro.test(dat1_clean$Fold_change)
kw_result <- kruskal.test(Fold_change ~ Treatment, data = dat1_clean)
kw_stat   <- kw_result$statistic
kw_p      <- kw_result$p.value

# Bonferroni-corrected pairwise comparisons
pw    <- pairwise.wilcox.test(dat1_clean$Fold_change, dat1_clean$Treatment,
                              p.adjust.method = "bonferroni")
pw_df <- as.data.frame(as.table(pw$p.value))
colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
  mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
  select(-temp)
pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                             ifelse(pw_df$p_adjusted < 0.01,  "**",
                                    ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))

# Raw (uncorrected) pairwise comparisons
pw_raw    <- pairwise.wilcox.test(dat1_clean$Fold_change, dat1_clean$Treatment,
                                  p.adjust.method = "none")
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

# Descriptive statistics
desc_stats <- dat1_clean %>%
  group_by(Treatment) %>%
  summarise(Mean = mean(Fold_change),
            SEM  = sd(Fold_change) / sqrt(n()),
            N    = n(), .groups = "drop")

stat_summary_df <- data.frame(
  Test           = c("Shapiro-Wilk (all groups)", "Kruskal-Wallis"),
  Statistic      = c(round(sw$statistic, 4), round(kw_stat, 4)),
  p_value        = c(round(sw$p.value, 4),   round(kw_p,   4)),
  Interpretation = c(
    ifelse(sw$p.value > 0.05, "Normal distribution", "Non-normal distribution"),
    ifelse(kw_p       < 0.05, "Significant",         "Not significant")
  )
)

# ── EXPORT CSV ───────────────────────────────────────────────────────────────

sec1 <- data.frame(
  Section      = c("Descriptive Statistics", rep("", nrow(desc_stats))),
  Treatment    = c("Treatment",  as.character(desc_stats$Treatment)),
  Mean         = c("Mean (AU)", as.character(round(desc_stats$Mean, 4))),
  SEM          = c("SEM",       as.character(round(desc_stats$SEM,  4))),
  N            = c("N",         as.character(desc_stats$N)),
  Significance = NA,
  Trend        = NA
)
sec2 <- data.frame(
  Section      = c("Normality & Omnibus Test", rep("", nrow(stat_summary_df))),
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
write.csv(rbind(sec1, spacer, sec2, spacer, sec3), "Fig_2b_stats.csv",
          row.names = FALSE)

# ── PLOT ANNOTATIONS ─────────────────────────────────────────────────────────

# Embed "ns#" directly for trend pairs — single bracket layer, matches Fig 2c-f
pw_plot <- pw_df %>%
  rename(group1 = Group1, group2 = Group2) %>%
  mutate(
    p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                          "ns#",
                          Significance),
    y.position   = c(1.55, 1.70, 1.85)
  )

# ── PLOT ─────────────────────────────────────────────────────────────────────

plot1 <- ggboxplot(dat1_clean, x = "Treatment", y = "Fold_change",
                   palette = "npg", color = "Treatment", fill = "Treatment",
                   alpha = 0.5, add = "jitter", ylim = c(-0.1, 2)) +
  labs(y = "Fold change in\n2-NBDG intensity (AU)") +
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
  stat_compare_means(label.y = 2, label.x.npc = "center", size = 2.5,
                     method = "kruskal.test", hjust = 0.5) +
  stat_pvalue_manual(pw_plot, label = "p.adj.signif",
                     tip.length = 0.03, size = 2.5, step.increase = 0.001) +
  stat_summary(fun.data = function(x) {
    data.frame(y = 0,
               label = paste0(round(mean(x), 2), " \u00b1 ",
                              round(sd(x) / sqrt(length(x)), 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot1

# ── SAVE ─────────────────────────────────────────────────────────────────────

ggsave("Fig_2b_Boxplot_Glucose_uptake.png", plot = plot1,
       scale = 1, width = 8, height = 5.3, units = "cm", dpi = 600,
       limitsize = TRUE)
ggsave("Fig_2b_Boxplot_Glucose_uptake.svg", plot = plot1,
       scale = 1, width = 8, height = 5.3, units = "cm", dpi = 600,
       limitsize = TRUE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 03: Fig_2b_glucose_uptake_revised.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 04/21  —  Fig_2cdef_functional_genes_revised.r                        ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 2c–2f: Glut2, Kcnj11, Cacna1c, Cacna1d fold change boxplots
# Lines: 304  |  File: Fig_2cdef_functional_genes_revised.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

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

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 04: Fig_2cdef_functional_genes_revised.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 05/21  —  Fig_3a_C-peptide_curve.r                                    ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 3a: Glucose-stimulated C-peptide secretion LOESS curve
# Lines: 82  |  File: Fig_3a_C-peptide_curve.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

rm(list = ls())

library(ggplot2)
library(readr)
library(ggsci)
library(readxl)
library(dplyr)

# ── LOAD DATA ────────────────────────────────────────────────────────────────

dat2 <- read_excel("Fig_3a_Scatterplot_C_peptide.xlsx", sheet = 1)

dat2$Glucose <- as.numeric(dat2$Glucose)
dat2 <- filter(dat2, Glucose != 44)

# ── EXPORT LOESS SMOOTH + CI PER TREATMENT ───────────────────────────────────

# Generate predicted values across the glucose range for each Treatment
glucose_seq <- seq(min(dat2$Glucose), max(dat2$Glucose), length.out = 200)

loess_export <- lapply(unique(dat2$Treatment), function(trt) {
  sub   <- dat2 %>% filter(Treatment == trt)
  model <- loess(C_peptide ~ Glucose, data = sub, span = 0.75)  # default ggplot span
  pred  <- predict(model, newdata = data.frame(Glucose = glucose_seq), se = TRUE)
  
  data.frame(
    Treatment  = trt,
    Assay      = unique(sub$Assay)[1],
    Glucose    = glucose_seq,
    Fitted     = pred$fit,
    SE         = pred$se.fit,
    CI_lower   = pred$fit - 1.96 * pred$se.fit,
    CI_upper   = pred$fit + 1.96 * pred$se.fit
  )
})

loess_df <- bind_rows(loess_export)

# Clip CI to plot limits (y >= 0) to match the plot
loess_df <- loess_df %>%
  mutate(CI_lower = pmax(CI_lower, 0),
         CI_upper = pmin(CI_upper, 320))

write.csv(loess_df, "Fig_3a_LOESS_smooth_CI.csv", row.names = FALSE)

# ── PLOT ─────────────────────────────────────────────────────────────────────

plot_loess <- ggplot(dat2, aes(x = Glucose, y = C_peptide, color = Treatment)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.5) +
  facet_grid(cols = vars(Assay)) +
  scale_color_npg() +
  labs(x = "Glucose concentration (mM)",
       y = "\nRelative C-peptide release (AU)") +
  lims(y = c(0, 320)) +
  theme_minimal() +
  theme(legend.position      = "top",
        legend.key.width     = unit(0.5, "cm"),
        legend.key.height    = unit(0.5, "cm"),
        strip.text           = element_text(size = 8, colour = "black", face = "bold"),
        legend.margin        = margin(0, 0, 0, 0),
        legend.box.spacing   = margin(5),
        legend.text          = element_text(size = 8),
        legend.title         = element_text(size = 8),
        axis.title.x         = element_text(size = 8),
        axis.text.x          = element_text(size = 8, colour = "black",
                                            angle = 0, vjust = 0.5),
        axis.text.y          = element_text(size = 8, colour = "black",
                                            angle = 90, hjust = 0.5),
        axis.title.y         = element_text(size = 8))

plot_loess

# ── SAVE ─────────────────────────────────────────────────────────────────────

ggsave("Fig_3a_Scatterplot_C_peptide.png", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

ggsave("Fig_3a_Scatterplot_C_peptide.svg", plot = plot_loess,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 05: Fig_3a_C-peptide_curve.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 06/21  —  Fig_3b_C-peptide_5.5mm_revised.r                            ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 3b: C-peptide secretion at 5.5 mM glucose — boxplot
# Lines: 164  |  File: Fig_3b_C-peptide_5.5mm_revised.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# C-peptide at Physiological Glucose - Figure 3b (Revised: Bonferroni-corrected p-values)

rm(list = ls())

library(ggplot2)
library(data.table)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)

# ── LOAD DATA ────────────────────────────────────────────────────────────────

dat <- read_excel("Fig_3b_Boxplot_C_peptide.xlsx", sheet = 1)
attach(dat)

dat1 <- dat
dat1 <- filter(dat1, Glucose == "5.5 mM glucose")

# ── STATISTICS ───────────────────────────────────────────────────────────────

sw        <- shapiro.test(dat1$C_peptide)   # W = 0.6920, p = 0.0007 (non-normal)
kw_result <- kruskal.test(C_peptide ~ Treatment, data = dat1)
kw_stat   <- kw_result$statistic
kw_p      <- kw_result$p.value

# Bonferroni-corrected pairwise comparisons
pw    <- pairwise.wilcox.test(dat1$C_peptide, dat1$Treatment,
                              p.adjust.method = "bonferroni")
pw_df <- as.data.frame(as.table(pw$p.value))
colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
  mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
  select(-temp)
pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                             ifelse(pw_df$p_adjusted < 0.01,  "**",
                                    ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))

# Raw (uncorrected) pairwise comparisons
pw_raw    <- pairwise.wilcox.test(dat1$C_peptide, dat1$Treatment,
                                  p.adjust.method = "none")
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
  summarise(Mean = mean(C_peptide), SEM = sd(C_peptide) / sqrt(n()),
            N = n(), .groups = "drop")

stat_summary_df <- data.frame(
  Test           = c("Shapiro-Wilk (all groups)", "Kruskal-Wallis"),
  Statistic      = c(round(sw$statistic, 4), round(kw_stat, 4)),
  p_value        = c(round(sw$p.value, 4),   round(kw_p,   4)),
  Interpretation = c(
    ifelse(sw$p.value > 0.05, "Normal distribution", "Non-normal distribution"),
    ifelse(kw_p       < 0.05, "Significant",         "Not significant")
  )
)

# ── EXPORT CSV ───────────────────────────────────────────────────────────────

sec1 <- data.frame(
  Section      = c("Descriptive Statistics", rep("", nrow(desc_stats))),
  Treatment    = c("Treatment",        as.character(desc_stats$Treatment)),
  Mean         = c("Mean (ng/µg DNA)", as.character(round(desc_stats$Mean, 4))),
  SEM          = c("SEM",              as.character(round(desc_stats$SEM,  4))),
  N            = c("N",                as.character(desc_stats$N)),
  Significance = NA,
  Trend        = NA
)
sec2 <- data.frame(
  Section      = c("Normality & Omnibus Test", rep("", nrow(stat_summary_df))),
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
write.csv(rbind(sec1, spacer, sec2, spacer, sec3), "Fig_3b_stats.csv",
          row.names = FALSE)

# ── PLOT ANNOTATIONS ─────────────────────────────────────────────────────────

glucose.labs <- c("5.5 mM glucose" = "5.5 mM glucose (Physiological level)")

dat1$Glucose <- factor(dat1$Glucose, levels = c("Basal medium", "2.8 mM glucose",
                                                "5.5 mM glucose", "22 mM glucose",
                                                "44 mM glucose"))

# Embed "ns#" for trend pairs — single bracket layer, matches Fig 2b/2c-f/1c/1d
pw_plot <- pw_df %>%
  rename(group1 = Group1, group2 = Group2) %>%
  mutate(
    p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                          "ns#",
                          Significance),
    y.position   = c(350, 380, 410)
  )

# ── PLOT ─────────────────────────────────────────────────────────────────────

plot1 <- ggboxplot(dat1, x = "Treatment", y = "C_peptide",
                   palette = "npg",
                   color = "Treatment", fill = "Treatment", alpha = 0.5,
                   add = "jitter",
                   ylim = c(0, 450)) +
  labs(y = "\nC-peptide release (ng/µg DNA)") +
  facet_grid(cols = vars(Glucose), labeller = labeller(Glucose = glucose.labs)) +
  theme_minimal() +
  theme(legend.position    = "top",
        legend.key.width   = unit(0.5, "cm"),
        legend.key.height  = unit(0.5, "cm"),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        legend.text        = element_text(size = 8),
        legend.title       = element_text(size = 8),
        strip.text         = element_text(size = 8, face = "bold"),
        axis.title.x       = element_text(size = 8, colour = "black"),
        axis.text.x        = element_text(size = 8, colour = "black",
                                          angle = 0, hjust = 0.5),
        axis.text.y        = element_text(size = 8, colour = "black",
                                          angle = 90, hjust = 0.5),
        axis.title.y       = element_text(size = 8)) +
  stat_compare_means(label.y = 448, label.x.npc = "center", size = 2.5,
                     method = "kruskal.test", hjust = 0.5) +
  stat_pvalue_manual(pw_plot, label = "p.adj.signif",
                     tip.length = 0.03, size = 2.5, step.increase = 0.001) +
  stat_summary(fun.data = function(x) {
    data.frame(y = 0,
               label = paste0(round(mean(x), 2), " \u00b1 ",
                              round(sd(x) / sqrt(length(x)), 2)))
  }, geom = "text", vjust = 1, size = 2.5)

plot1

ggsave("Fig_3b_Boxplot_C_peptide.png", plot = plot1,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)
ggsave("Fig_3b_Boxplot_C_peptide.svg", plot = plot1,
       scale = 1, width = 8, height = 6.5, units = "cm",
       dpi = 600, limitsize = TRUE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 06: Fig_3b_C-peptide_5.5mm_revised.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 07/21  —  Fig_3ef_Boxplot_script.r                                    ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 3c–3f: Glp1r, Ins1, Ptbp1, Itpr1, Ryr1–3 gene expression boxplots
# Lines: 359  |  File: Fig_3ef_Boxplot_script.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Calcium Signalling & RNA-Binding Gene Expression — Figure 3e–f
# Fig 3e: Ryr2, Ryr3, Itpr1 (Ryr1 excluded — near-undetectable CTRL expression)
# Fig 3f: Ptbp1
# Dimensions: Fig 3e = 16.2 × 5.3 cm | Fig 3f = 8 × 5.3 cm

rm(list = ls())

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggsci)
library(readxl)
library(rstatix)
library(patchwork)

`%||%` <- function(a, b) if (!is.null(a)) a else b

npg_colors <- pal_npg("nrc")(10)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

dat_3e <- read_excel("Fig_3ef_Boxplot_Fold_change_mRNA_expression.xlsx", sheet = 1)
dat_3f <- read_excel("Fig_3ef_Boxplot_Fold_change_mRNA_expression.xlsx", sheet = 2)

gene_order_3e <- c("Ryr2", "Ryr3", "Itpr1")
treat_levels  <- c("CTRL", "DAPT", "DKK-1")

dat_3e <- dat_3e %>%
  mutate(Gene        = factor(Gene,      levels = gene_order_3e),
         Treatment   = factor(Treatment, levels = treat_levels),
         Fold_change = as.numeric(Fold_change))

dat_3f <- dat_3f %>%
  mutate(Treatment   = factor(Treatment, levels = treat_levels),
         Fold_change = as.numeric(Fold_change))

# ── HELPER: NORMALITY SCREEN ──────────────────────────────────────────────────

safe_shapiro_p <- function(x) {
  x <- x[!is.na(x)]
  if (length(unique(x)) < 3 || length(x) < 3) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value,
           error   = function(e) NA_real_,
           warning = function(w) NA_real_)
}

check_parametric <- function(dat, genes) {
  sapply(genes, function(g) {
    pvals <- dat %>%
      filter(Gene == g, Treatment != "CTRL") %>%
      group_by(Treatment) %>%
      summarise(p = safe_shapiro_p(Fold_change), .groups = "drop") %>%
      pull(p)
    all(pvals > 0.05, na.rm = TRUE)
  })
}

# ── HELPER: STATISTICS + CSV EXPORT ──────────────────────────────────────────

run_stats <- function(df, gene_name, parametric = TRUE) {
  
  x   <- df$Fold_change
  grp <- df$Treatment
  
  sw <- tryCatch(shapiro.test(x),
                 error   = function(e) list(statistic = NA, p.value = NA),
                 warning = function(w) list(statistic = NA, p.value = NA))
  
  if (parametric) {
    fit           <- lm(Fold_change ~ Treatment, data = df)
    fit_anv       <- anova(fit)
    omni_stat     <- round(as.numeric(fit_anv[["F value"]][1]), 4)
    omni_p        <- round(as.numeric(fit_anv[["Pr(>F)"]][1]),  4)
    omni_name     <- "One-way ANOVA (F)"
    pw            <- pairwise.t.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.t.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise t-test (Bonferroni)"
  } else {
    omni          <- kruskal.test(Fold_change ~ Treatment, data = df)
    omni_stat     <- round(as.numeric(omni$statistic), 4)
    omni_p        <- round(as.numeric(omni$p.value),   4)
    omni_name     <- "Kruskal-Wallis (H)"
    pw            <- pairwise.wilcox.test(x, grp, p.adjust.method = "bonferroni")
    pw_raw        <- pairwise.wilcox.test(x, grp, p.adjust.method = "none")
    post_hoc_name <- "Pairwise Wilcoxon (Bonferroni)"
  }
  
  # Bonferroni-corrected
  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("Group1", "Group2", "p_adjusted")
  pw_df <- pw_df[!is.na(pw_df$p_adjusted), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  pw_df$Significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                               ifelse(pw_df$p_adjusted < 0.01,  "**",
                                      ifelse(pw_df$p_adjusted < 0.05,  "*", "ns")))
  
  # Raw (uncorrected)
  pw_raw_df <- as.data.frame(as.table(pw_raw$p.value))
  colnames(pw_raw_df) <- c("Group1", "Group2", "p_raw")
  pw_raw_df <- pw_raw_df[!is.na(pw_raw_df$p_raw), ] %>%
    mutate(temp = Group1, Group1 = Group2, Group2 = temp) %>%
    select(-temp)
  
  # Merge and flag trends
  pw_df <- merge(pw_df, pw_raw_df[, c("Group1", "Group2", "p_raw")],
                 by = c("Group1", "Group2"), all.x = TRUE)
  pw_df$Trend <- ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted >= 0.05, "#",
                        ifelse(pw_df$p_raw < 0.05 & pw_df$p_adjusted < 0.05,
                               pw_df$Significance, "ns"))
  
  desc <- df %>%
    group_by(Treatment) %>%
    summarise(Mean = mean(Fold_change), SEM = sd(Fold_change) / sqrt(n()),
              N = n(), .groups = "drop")
  
  stat_df <- data.frame(
    Test = c(paste0("Shapiro-Wilk (", gene_name, ")"), omni_name),
    Statistic      = c(round(as.numeric(sw$statistic), 4), omni_stat),
    p_value        = c(round(as.numeric(sw$p.value),   4), omni_p),
    Interpretation = c(
      ifelse(is.na(sw$p.value), "Undetermined",
             ifelse(sw$p.value > 0.05, "Normal", "Non-normal")),
      ifelse(is.na(omni_p), "Undetermined",
             ifelse(omni_p < 0.05, "Significant", "Not significant"))
    )
  )
  
  sec1 <- data.frame(
    Section = c(paste0("Descriptive Statistics - ", gene_name), rep("", nrow(desc))),
    Col1    = c("Treatment",  as.character(desc$Treatment)),
    Col2    = c("Mean (AU)",  as.character(round(desc$Mean, 4))),
    Col3    = c("SEM",        as.character(round(desc$SEM,  4))),
    Col4    = c("N",          as.character(desc$N)),
    Col5    = NA, Col6 = NA
  )
  sec2 <- data.frame(
    Section = c(paste0("Normality & Omnibus - ", gene_name), rep("", nrow(stat_df))),
    Col1    = c("Test",           as.character(stat_df$Test)),
    Col2    = c("Statistic",      as.character(stat_df$Statistic)),
    Col3    = c("p_value",        as.character(stat_df$p_value)),
    Col4    = c("Interpretation", as.character(stat_df$Interpretation)),
    Col5    = NA, Col6 = NA
  )
  sec3 <- data.frame(
    Section = c(paste0("Pairwise - ", gene_name, " (", post_hoc_name, ")"),
                rep("", nrow(pw_df))),
    Col1    = c("Group1",       as.character(pw_df$Group1)),
    Col2    = c("Group2",       as.character(pw_df$Group2)),
    Col3    = c("p_raw",        as.character(round(pw_df$p_raw,      6))),
    Col4    = c("p_adjusted",   as.character(round(pw_df$p_adjusted, 6))),
    Col5    = c("Significance", as.character(pw_df$Significance)),
    Col6    = c("Trend_flag",   as.character(pw_df$Trend))
  )
  spacer <- data.frame(Section = "---", Col1 = NA, Col2 = NA,
                       Col3 = NA, Col4 = NA, Col5 = NA, Col6 = NA)
  out    <- rbind(sec1, spacer, sec2, spacer, sec3, spacer)
  
  list(csv = out, pw = pw_df)
}

# ── HELPER: build pw_plot with ns# embedded ──────────────────────────────────

make_pw_plot <- function(pw_df, y_positions) {
  pw_df %>%
    rename(group1 = Group1, group2 = Group2) %>%
    mutate(
      p.adj.signif = ifelse(p_raw < 0.05 & p_adjusted >= 0.05,
                            "ns#",
                            Significance),
      y.position   = y_positions
    )
}

# ── HELPER: BOXPLOT BUILDER ───────────────────────────────────────────────────

make_boxplot <- function(df, pw_plot, gene_lab, gene_key,
                         omnibus_method, y_lab_pos, y_lim_top,
                         show_y_title = TRUE,
                         summary_y    = 0.10) {
  
  labs_vec <- setNames(gene_lab, gene_key)
  
  p <- ggboxplot(df,
                 x = "Treatment", y = "Fold_change",
                 palette = "npg",
                 color = "Treatment", fill = "Treatment", alpha = 0.5,
                 add = "jitter",
                 ylim = c(-0.1, y_lim_top)) +
    labs(y = if (show_y_title) "Fold change in\nmRNA expression (AU)" else NULL) +
    facet_wrap(vars(Gene), labeller = labeller(Gene = labs_vec)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 8, colour = "black", face = "bold.italic"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = 8, colour = "black",
                                     angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y     = element_text(size = 8, colour = "black",
                                     angle = 90, hjust = 0.5),
      axis.title.y    = if (show_y_title) element_text(size = 8) else element_blank()
    ) +
    stat_compare_means(
      label.y = y_lab_pos, label.x.npc = "center",
      size = 2.5, method = omnibus_method, hjust = 0.5
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(y     = summary_y,
                   label = paste0(round(mean(x), 2), " \u00b1 ",
                                  round(sd(x) / sqrt(length(x)), 2)))
      },
      geom = "text", vjust = 1, size = 2.5, angle = 0
    )
  
  if (nrow(pw_plot) > 0) {
    p <- p + stat_pvalue_manual(pw_plot, label = "p.adj.signif",
                                tip.length = 0.03, size = 2.5)
  }
  return(p)
}

save_plot <- function(plot, filename, width_cm, height_cm) {
  ggsave(paste0(filename, ".png"), plot = plot, scale = 1,
         width = width_cm, height = height_cm, units = "cm",
         dpi = 600, limitsize = TRUE)
  ggsave(paste0(filename, ".svg"), plot = plot, scale = 1,
         width = width_cm, height = height_cm, units = "cm",
         dpi = 600, limitsize = TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3e — Ryr2, Ryr3, Itpr1  (16.2 × 5.3 cm)  FREE Y-AXIS
# Note: Ryr1 excluded — CTRL 2^-dCt ~1.7-2.0e-5 (near-undetectable);
#       fold-change quantification unreliable. Raw data retained in Table S2.
# ══════════════════════════════════════════════════════════════════════════════

parametric_map_3e <- check_parametric(dat_3e, gene_order_3e)

cat("Normality screen — Fig 3e:\n")
print(data.frame(
  Gene       = gene_order_3e,
  Parametric = unname(parametric_map_3e),
  Test       = ifelse(unname(parametric_map_3e),
                      "One-way ANOVA + Bonferroni",
                      "Kruskal-Wallis + Wilcoxon Bonferroni")
))

gene_dfs_3e <- lapply(gene_order_3e, function(g) filter(dat_3e, Gene == g))
names(gene_dfs_3e) <- gene_order_3e

res_3e <- lapply(gene_order_3e, function(g) {
  run_stats(gene_dfs_3e[[g]], g, unname(parametric_map_3e[g]))
})
names(res_3e) <- gene_order_3e

combined_csv_3e <- do.call(rbind, lapply(res_3e, `[[`, "csv"))
write.csv(combined_csv_3e, "Fig_3e_stats.csv", row.names = FALSE)

gene_labels_3e <- c(
  "Ryr2"  = "Ryr2 (ER Ca\u00b2\u207a release)",
  "Ryr3"  = "Ryr3 (ER Ca\u00b2\u207a release)",
  "Itpr1" = "Itpr1 (ER Ca\u00b2\u207a release)"
)

plot_list_3e <- list()

for (g in gene_order_3e) {
  df       <- gene_dfs_3e[[g]]
  pw_df    <- res_3e[[g]]$pw
  n_pairs  <- nrow(pw_df)
  omni_m   <- ifelse(unname(parametric_map_3e[g]), "anova", "kruskal.test")
  y_max    <- max(df$Fold_change, na.rm = TRUE)
  
  y_pos_g   <- if (n_pairs > 0)
    seq(y_max * 1.12, y_max * 1.50, length.out = n_pairs) else numeric(0)
  y_lim_top <- if (n_pairs > 0) max(y_pos_g) * 1.15 else y_max * 1.65
  y_lab_pos <- y_lim_top * 0.97
  summary_y <- y_max * 0.05
  
  pw_plot <- make_pw_plot(pw_df, y_pos_g)
  
  plot_list_3e[[g]] <- make_boxplot(
    df             = df,
    pw_plot        = pw_plot,
    gene_lab       = gene_labels_3e[g],
    gene_key       = g,
    omnibus_method = omni_m,
    y_lab_pos      = y_lab_pos,
    y_lim_top      = y_lim_top,
    show_y_title   = (g == "Ryr2"),
    summary_y      = summary_y
  )
  
  cat("Plot built:", g, "| y_max =", round(y_max, 3),
      "| y_lim_top =", round(y_lim_top, 3),
      "| n_brackets =", n_pairs, "\n")
}

for (g in gene_order_3e) {
  save_plot(plot_list_3e[[g]],
            paste0("Fig_3e_Boxplot_", g, "_expression"),
            width_cm = 4.05, height_cm = 5.3)
}

combined_3e <- wrap_plots(plot_list_3e, ncol = 3, nrow = 1) &
  theme(plot.margin = margin(2.5, 2.5, 2.5, 2.5))

save_plot(combined_3e,
          "Fig_3e_Boxplot_Ryr2_Ryr3_Itpr1_expression",
          width_cm = 16.2, height_cm = 5.3)

cat("\nFig 3e saved (Ryr2, Ryr3, Itpr1 — Ryr1 excluded).\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3f — Ptbp1  (8 × 5.3 cm)
# ══════════════════════════════════════════════════════════════════════════════

sw_ptbp1_dapt <- safe_shapiro_p(filter(dat_3f, Treatment == "DAPT")$Fold_change)
sw_ptbp1_dkk1 <- safe_shapiro_p(filter(dat_3f, Treatment == "DKK-1")$Fold_change)
is_param_3f   <- all(c(sw_ptbp1_dapt, sw_ptbp1_dkk1) > 0.05, na.rm = TRUE)

cat("\nNormality screen — Fig 3f (Ptbp1):",
    ifelse(is_param_3f, "Parametric (ANOVA)", "Non-parametric (KW)"), "\n")

res_3f <- run_stats(dat_3f, "Ptbp1", is_param_3f)
write.csv(res_3f$csv, "Fig_3f_Ptbp1_stats.csv", row.names = FALSE)

pw_df_3f   <- res_3f$pw
y_max_3f   <- max(dat_3f$Fold_change, na.rm = TRUE)
n_pairs_3f <- nrow(pw_df_3f)
y_pos_3f   <- if (n_pairs_3f > 0)
  seq(y_max_3f * 1.12, y_max_3f * 1.50, length.out = n_pairs_3f) else numeric(0)
y_lim_3f   <- if (n_pairs_3f > 0) max(y_pos_3f) * 1.15 else y_max_3f * 1.65
y_lab_3f   <- y_lim_3f * 0.97

pw_plot_3f <- make_pw_plot(pw_df_3f, y_pos_3f)

Ptbp1.labs <- c("Ptbp1" = "Ptbp1 (mRNA stability & translation)")

plot.3f <- make_boxplot(
  df             = dat_3f,
  pw_plot        = pw_plot_3f,
  gene_lab       = Ptbp1.labs["Ptbp1"],
  gene_key       = "Ptbp1",
  omnibus_method = ifelse(is_param_3f, "anova", "kruskal.test"),
  y_lab_pos      = y_lab_3f,
  y_lim_top      = y_lim_3f,
  show_y_title   = TRUE,
  summary_y      = y_max_3f * 0.05
)

save_plot(plot.3f,
          "Fig_3f_Boxplot_Ptbp1_expression",
          width_cm = 8, height_cm = 5.3)

cat("Fig 3f saved.\n")
cat("\nAll files saved.\nStats: Fig_3e_stats.csv | Fig_3f_Ptbp1_stats.csv\n")

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 07: Fig_3ef_Boxplot_script.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 08/21  —  Fig_3g_PCA_script.r                                         ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 3g: PCA of 11 secretory-pathway genes (biplot)
# Lines: 110  |  File: Fig_3g_PCA_script.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

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

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 08: Fig_3g_PCA_script.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 09/21  —  Fig_4b_correlation_matrix.r                                 ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 4b: Pearson correlation matrix heatmap (Hes1-normalized)
# Lines: 59  |  File: Fig_4b_correlation_matrix.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Gene Correlation Matrix - Figure 4b
 
rm(list = ls())

library(ggplot2)
library(readxl)
library(ggcorrplot)
library(dplyr)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load data
dat <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

# Calculate correlations
corr <- round(cor(dat, use = "pairwise.complete.obs"), 2)
p.mat <- round(cor_pmat(dat), 2)

# Extract significant correlations only
sig_corr <- corr
sig_corr[p.mat > 0.05] <- NA

# Export
write.csv(corr, "Fig_4b_full_correlation_matrix.csv", row.names = TRUE)
write.csv(p.mat, "Fig_4b_pvalue_matrix.csv", row.names = TRUE)
write.csv(sig_corr, "Fig_4b_significant_correlations.csv", row.names = TRUE)

# Plot
panel_b <- ggcorrplot(sig_corr,
                      hc.order = FALSE,
                      outline.col = "white",
                      type = "full",
                      show.legend = TRUE) +
  scale_fill_gradientn(colors = c(npg_colors[4], "white", npg_colors[1]),
                       values = scales::rescale(c(-1, 0, 1)),
                       limits = c(-1, 1), 
                       name = "Pearson correlation\ncoefficient (r)",
                       na.value = "grey85") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = margin(5),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8, family = "Arial"),
        legend.title = element_text(size = 8, colour = "black", angle = 0,
                                   vjust = 0.95, hjust = 1, family = "Arial"),
        axis.text.x = element_text(size = 8, colour = "black", face = "italic",
                                   angle = 90, hjust = 1, vjust = 0.5, family = "Arial"),
        axis.text.y = element_text(size = 8, colour = "black", face = "italic",
                                   family = "Arial"),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

ggsave("Fig_4b_Correlation_matrix.png", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

ggsave("Fig_4b_Correlation_matrix.svg", plot = panel_b,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 09: Fig_4b_correlation_matrix.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 10/21  —  Fig_4c_normalization_methods.r                              ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 4c: Multi-method normalization comparison — 4 strategies
# Lines: 118  |  File: Fig_4c_normalization_methods.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# =============================================================================
# Normalization Methods Comparison - Figure 4c
# Revised: Hmisc replaced with base R cor.test(); Ryr1 excluded (QC failure)
# =============================================================================

rm(list = ls())

library(ggplot2)
library(readxl)
library(reshape2)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load data ──────────────────────────────────────────────────────────────────
gene_expr_data <- read_excel("Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx", sheet = 1)

# ── Exclude Ryr1 (non-specific melt curve; QC failure) ────────────────────────
gene_expr_data <- gene_expr_data[, colnames(gene_expr_data) != "Ryr1"]

# ── Correlation analysis function (base R, no Hmisc) ──────────────────────────
perform_correlation_analysis <- function(data, method_name) {
  data_clean <- as.matrix(na.omit(data))
  n_genes <- ncol(data_clean)
  sig_count <- 0
  strong_count <- 0
  very_strong_count <- 0
  
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      test  <- cor.test(data_clean[, i], data_clean[, j], method = "pearson")
      p_val <- test$p.value
      r_val <- abs(test$estimate)
      if (!is.na(p_val) && p_val < 0.05) {
        sig_count <- sig_count + 1
        if (r_val > 0.7) strong_count <- strong_count + 1
        if (r_val > 0.8) very_strong_count <- very_strong_count + 1
      }
    }
  }
  
  return(data.frame(
    Method    = method_name,
    Significant  = sig_count,
    Strong       = strong_count,
    VeryStrong   = very_strong_count
  ))
}

# ── Hes1-normalized ───────────────────────────────────────────────────────────
hes1_norm <- gene_expr_data
for (col in colnames(hes1_norm)) {
  if (col != "Hes1") hes1_norm[[col]] <- hes1_norm[[col]] / hes1_norm$Hes1
}
hes1_norm <- hes1_norm[, colnames(hes1_norm) != "Hes1"]
results_list <- list()
results_list[[1]] <- perform_correlation_analysis(hes1_norm, "Hes1-\nnormalized")
write.csv(hes1_norm, "Fig_4c_Hes1_normalized_data.csv", row.names = FALSE)

# ── Kcnj11-normalized ─────────────────────────────────────────────────────────
kcnj11_norm <- gene_expr_data
for (col in colnames(kcnj11_norm)) {
  if (col != "Kcnj11") kcnj11_norm[[col]] <- kcnj11_norm[[col]] / kcnj11_norm$Kcnj11
}
kcnj11_norm <- kcnj11_norm[, colnames(kcnj11_norm) != "Kcnj11"]
results_list[[2]] <- perform_correlation_analysis(kcnj11_norm, "Kcnj11-\nnormalized")
write.csv(kcnj11_norm, "Fig_4c_Kcnj11_normalized_data.csv", row.names = FALSE)

# ── Ins1-normalized ───────────────────────────────────────────────────────────
ins1_norm <- gene_expr_data
for (col in colnames(ins1_norm)) {
  if (col != "Ins1") ins1_norm[[col]] <- ins1_norm[[col]] / ins1_norm$Ins1
}
ins1_norm <- ins1_norm[, colnames(ins1_norm) != "Ins1"]
results_list[[3]] <- perform_correlation_analysis(ins1_norm, "Ins1-\nnormalized")
write.csv(ins1_norm, "Fig_4c_Ins1_normalized_data.csv", row.names = FALSE)

# ── Raw expression ────────────────────────────────────────────────────────────
results_list[[4]] <- perform_correlation_analysis(gene_expr_data, "Raw\nexpression")
write.csv(gene_expr_data, "Fig_4c_Raw_expression_data.csv", row.names = FALSE)

# ── Combine and export results ────────────────────────────────────────────────
comparison_data <- do.call(rbind, results_list)
write.csv(comparison_data, "Fig_4c_normalization_comparison_data.csv", row.names = FALSE)

# ── Reshape for plotting ──────────────────────────────────────────────────────
comparison_long <- melt(comparison_data, id.vars = "Method")
comparison_long$Method <- factor(comparison_long$Method, levels = comparison_data$Method)

# ── Plot ──────────────────────────────────────────────────────────────────────
panel_c <- ggplot(comparison_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = value),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 2.5, family = "Arial") +
  scale_fill_manual(
    values = c(npg_colors[4], npg_colors[2], npg_colors[3]),
    labels = c("Significant\n(p<0.05)", "Strong\n(|r|>0.7)", "Very Strong\n(|r|>0.8)")
  ) +
  labs(x = NULL, y = "Count", fill = NULL) +
  theme_minimal() +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(size = 8, color = "black", family = "Arial"),
    axis.text.y      = element_text(size = 8, color = "black", family = "Arial"),
    axis.title.y     = element_text(size = 8, face = "bold", family = "Arial"),
    legend.text      = element_text(size = 8, family = "Arial"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# ── Save outputs ──────────────────────────────────────────────────────────────
ggsave("Fig_4c_Normalization_comparison.png", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

ggsave("Fig_4c_Normalization_comparison.svg", plot = panel_c,
       width = 8.1, height = 8.7, units = "cm", dpi = 600)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 10: Fig_4c_normalization_methods.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 11/21  —  Fig_4d_Venn_diagram.r                                       ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 4d: Four-way Venn diagram — 28-pair consensus network
# Lines: 98  |  File: Fig_4d_Venn_diagram.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# =============================================================================
# Four-way Venn Diagram - Figure 4d
# Revised: Hmisc replaced with base R cor.test(); Ryr1 excluded (QC failure)
# =============================================================================

rm(list = ls())

library(VennDiagram)
library(grid)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load normalized data (already Ryr1-excluded from Fig 4c script) ───────────
hes1_norm   <- read.csv("Fig_4c_Hes1_normalized_data.csv")
kcnj11_norm <- read.csv("Fig_4c_Kcnj11_normalized_data.csv")
ins1_norm   <- read.csv("Fig_4c_Ins1_normalized_data.csv")
raw_data    <- read.csv("Fig_4c_Raw_expression_data.csv")

# ── Safety check: drop Ryr1 if still present ──────────────────────────────────
drop_ryr1 <- function(df) df[, colnames(df) != "Ryr1"]
hes1_norm   <- drop_ryr1(hes1_norm)
kcnj11_norm <- drop_ryr1(kcnj11_norm)
ins1_norm   <- drop_ryr1(ins1_norm)
raw_data    <- drop_ryr1(raw_data)

# ── Function to extract significant pairs (base R, no Hmisc) ─────────────────
get_pairs <- function(data) {
  data_clean <- as.matrix(na.omit(data))
  n_genes    <- ncol(data_clean)
  genes      <- colnames(data_clean)
  pairs      <- character()
  
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      test  <- cor.test(data_clean[, i], data_clean[, j], method = "pearson")
      p_val <- test$p.value
      if (!is.na(p_val) && p_val < 0.05) {
        pairs <- c(pairs, paste(genes[i], genes[j], sep = "_"))
      }
    }
  }
  return(pairs)
}

# ── Extract pairs from each method ────────────────────────────────────────────
pairs_list <- list(
  "Hes1-\nnormalized"   = get_pairs(hes1_norm),
  "Kcnj11-\nnormalized" = get_pairs(kcnj11_norm),
  "Ins1-\nnormalized"   = get_pairs(ins1_norm),
  "Raw\nexpression"     = get_pairs(raw_data)
)

# ── Identify consensus pairs (present in all 4 methods) ───────────────────────
consensus_pairs <- Reduce(intersect, pairs_list)

# ── Export pair lists and consensus ───────────────────────────────────────────
write.csv(data.frame(Pairs = pairs_list[[1]]),
          "Fig_4d_Hes1_normalized_pairs.csv",   row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[2]]),
          "Fig_4d_Kcnj11_normalized_pairs.csv", row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[3]]),
          "Fig_4d_Ins1_normalized_pairs.csv",   row.names = FALSE)
write.csv(data.frame(Pairs = pairs_list[[4]]),
          "Fig_4d_Raw_expression_pairs.csv",    row.names = FALSE)
write.csv(data.frame(Consensus_Pairs = consensus_pairs),
          "Fig_4d_consensus_pairs.csv",         row.names = FALSE)

# ── Build Venn diagram object ─────────────────────────────────────────────────
venn.plot <- venn.diagram(
  x              = pairs_list,
  category.names = names(pairs_list),
  filename       = NULL,
  fill           = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  col            = c(npg_colors[4], npg_colors[2], npg_colors[3], npg_colors[7]),
  lwd            = 1.5,
  alpha          = 0.5,
  cex            = 0.8,
  fontface       = "bold",
  fontfamily     = "Arial",
  cat.cex        = 0.7,
  cat.fontface   = "bold",
  cat.fontfamily = "Arial",
  cat.default.pos = "outer",
  cat.dist       = c(0.08, 0.08, 0.08, 0.08)
)

# ── Save PNG ──────────────────────────────────────────────────────────────────
png("Fig_4d_Venn_diagram.png",
    width = 9.16, height = 9.16, units = "cm", res = 600)
grid.draw(venn.plot)
dev.off()

# ── Save SVG ──────────────────────────────────────────────────────────────────
svg("Fig_4d_Venn_diagram.svg",
    width = 9.16 / 2.54, height = 9.16 / 2.54)
grid.draw(venn.plot)
dev.off()

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 11: Fig_4d_Venn_diagram.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 12/21  —  Fig_4e_Hes1-normalized_correlation_matrix.R                 ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 4e: Hes1-normalized annotated correlation heatmap (consensus pairs marked)
# Lines: 213  |  File: Fig_4e_Hes1-normalized_correlation_matrix.R  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

################################################################################
# FIGURE 4e - HES1-NORMALIZED CORRELATION MATRIX
# SIGNIFICANT PAIRS ONLY + ANNOTATE 28 CONSENSUS PAIRS WITH ASTERISKS
# Revised: Hmisc removed; Ryr1 excluded (QC failure); cor/cor_pmat used only
################################################################################

rm(list = ls())

library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

################################################################################
# LOAD DATA
################################################################################

cat("\nLoading Hes1-normalized data...\n")
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv")

# Safety: drop Ryr1 if present
hes1_norm <- hes1_norm[, colnames(hes1_norm) != "Ryr1"]

cat("Genes:", ncol(hes1_norm), "\n")
cat("Samples:", nrow(hes1_norm), "\n")

cat("\nLoading consensus pairs...\n")
consensus_df    <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
consensus_pairs <- consensus_df$Consensus_Pairs
cat("Consensus pairs loaded:", length(consensus_pairs), "\n")
print(consensus_pairs)

################################################################################
# CALCULATE CORRELATIONS
################################################################################

cat("\nCalculating correlation matrix...\n")

corr_e <- cor(hes1_norm, use = "pairwise.complete.obs")
p_mat_e <- cor_pmat(hes1_norm)

cat("Correlation matrix dimensions:", dim(corr_e), "\n")

################################################################################
# FILTER TO SIGNIFICANT CORRELATIONS ONLY
################################################################################

cat("\nFiltering significant correlations...\n")

sig_corr_e <- corr_e
sig_corr_e[p_mat_e >= 0.05] <- NA

n_sig <- sum(p_mat_e < 0.05, na.rm = TRUE) / 2
cat("Significant pairs (p < 0.05):", n_sig, "\n")

################################################################################
# CREATE CONSENSUS PAIR ANNOTATION MATRIX
################################################################################

cat("\nCreating consensus pair annotations...\n")

genes <- rownames(corr_e)
consensus_labels <- matrix("", nrow = nrow(corr_e), ncol = ncol(corr_e))
rownames(consensus_labels) <- genes
colnames(consensus_labels) <- genes

for (pair in consensus_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  gene1 <- pair_genes[1]
  gene2 <- pair_genes[2]
  idx1  <- which(genes == gene1)
  idx2  <- which(genes == gene2)
  if (length(idx1) > 0 && length(idx2) > 0) {
    consensus_labels[idx1, idx2] <- "*"
    consensus_labels[idx2, idx1] <- "*"
  }
}

n_consensus_annotated <- sum(consensus_labels == "*") / 2
cat("Consensus pairs annotated:", n_consensus_annotated, "\n")

################################################################################
# EXPORT DATA
################################################################################

cat("\nExporting data...\n")

write.csv(as.data.frame(sig_corr_e),
          "Fig_4e_Hes1_significant_correlations_only.csv", row.names = TRUE)

consensus_matrix_df <- as.data.frame(consensus_labels)
write.csv(consensus_matrix_df,
          "Fig_4e_consensus_annotations.csv", row.names = TRUE)

cat("Exported: Fig_4e_Hes1_significant_correlations_only.csv\n")
cat("Exported: Fig_4e_consensus_annotations.csv\n")

################################################################################
# VERIFY KEY PAIRS
################################################################################

cat("\n=== VERIFICATION: Sample correlation values ===\n")
test_pairs <- c("Cacna1d_Wnt5a", "Cacna1d_Ptbp1", "Ptbp1_Wnt5a")
for (pair in test_pairs) {
  pair_genes <- strsplit(pair, "_")[[1]]
  if (pair_genes[1] %in% genes && pair_genes[2] %in% genes) {
    val   <- corr_e[pair_genes[1], pair_genes[2]]
    p_val <- p_mat_e[pair_genes[1], pair_genes[2]]
    cat(sprintf("%s: r = %.10f (p = %.4e)\n", pair, val, p_val))
  }
}
cat("===============================================\n\n")

################################################################################
# PREPARE PLOT DATA
################################################################################

cat("\nPreparing plot data...\n")

corr_e_display          <- round(corr_e, 2)
sig_corr_e_display      <- corr_e_display
sig_corr_e_display[p_mat_e >= 0.05] <- NA

corr_long     <- melt(sig_corr_e_display)
consensus_long <- melt(consensus_labels)

plot_data <- data.frame(
  Gene1       = corr_long$Var1,
  Gene2       = corr_long$Var2,
  Correlation = corr_long$value,
  Consensus   = consensus_long$value
)
plot_data <- plot_data[!is.na(plot_data$Correlation), ]

cat("Total data points to plot:", nrow(plot_data), "\n")

################################################################################
# CREATE PLOT
################################################################################

cat("\nGenerating Figure 4e...\n")

panel_e <- ggplot(plot_data, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.2) +
  geom_text(aes(label = Consensus),
            size = 3, color = "black", fontface = "bold", family = "Arial") +
  scale_fill_gradient2(
    low      = npg_colors[4],
    mid      = "white",
    high     = npg_colors[1],
    midpoint = 0,
    limits   = c(-1, 1),
    na.value = "grey85",
    name     = "Pearson correlation\ncoefficient (r)"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    legend.position      = "top",
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = margin(5),
    legend.justification = c("right", "top"),
    legend.text          = element_text(size = 8, family = "Arial"),
    legend.title         = element_text(size = 8, colour = "black", angle = 0,
                                        vjust = 0.95, hjust = 1, family = "Arial"),
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    axis.text.x          = element_text(size = 8, colour = "black", face = "italic",
                                        angle = 90, hjust = 1, vjust = 0.5,
                                        family = "Arial"),
    axis.text.y          = element_text(size = 8, colour = "black", face = "italic",
                                        family = "Arial"),
    panel.grid           = element_blank()
  )

################################################################################
# SAVE FIGURE
################################################################################

cat("\nSaving figure...\n")

ggsave("Fig_4e_Significant_pairs_consensus_annotated.png", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)
ggsave("Fig_4e_Significant_pairs_consensus_annotated.svg", plot = panel_e,
       width = 9.16, height = 9.16, units = "cm", dpi = 600)

################################################################################
# SUMMARY
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("FIGURE 4e COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")
cat("Data source    : Hes1-normalized gene expression\n")
cat("Genes analysed :", ncol(hes1_norm), "(Ryr1 excluded - QC failure)\n")
cat("Total pairs    :", nrow(corr_e) * (nrow(corr_e) - 1) / 2, "\n")
cat("Significant    :", n_sig, "(p < 0.05)\n")
cat("Consensus      :", n_consensus_annotated, "(all 4 methods)\n")
cat("% consensus/sig:",
    round(100 * n_consensus_annotated / n_sig, 1), "%\n\n")
cat("Exports:\n")
cat("  Fig_4e_Hes1_significant_correlations_only.csv (full precision)\n")
cat("  Fig_4e_consensus_annotations.csv\n")
cat("  Fig_4e_Significant_pairs_consensus_annotated.png (600 DPI)\n")
cat("  Fig_4e_Significant_pairs_consensus_annotated.svg\n\n")
cat("NOTE: CSV exported at full precision; figure displays rounded (2 dp).\n")
cat(rep("=", 60), "\n\n", sep = "")

################################################################################
# END OF SCRIPT
################################################################################

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 12: Fig_4e_Hes1-normalized_correlation_matrix.R
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 13/21  —  Fig_S1_validation_analysis.r                                ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig S1: Validation analysis — observed vs null distribution, false positive rates, sign-flip check
# Lines: 491  |  File: Fig_S1_validation_analysis.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# =============================================================================
# Fig. S1 — Multi-Normalization Consensus Validation (Simplified, 3 panels)
#
# DATA DEPENDENCIES:
#   Fig_4c_Raw_expression_data.csv
#   Fig_4c_Ins1_normalized_data.csv
#   Fig_4c_Kcnj11_normalized_data.csv
#   Fig_4c_Hes1_normalized_data.csv
#   Fig_4d_consensus_pairs.csv
#
# Run order: Fig_4c → Fig_4d → this script
# =============================================================================

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_Fig_S1 <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 9),
    axis.line          = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks         = element_line(linewidth = 0.3, colour = "grey30"),
    axis.text          = element_text(colour = "grey20", family = "Arial"),
    axis.title         = element_text(colour = "grey20", family = "Arial"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.key.size    = unit(0.35, "cm"),
    legend.text        = element_text(size = 8, family = "Arial"),
    legend.title       = element_text(size = 8, face = "bold", family = "Arial"),
    plot.title         = element_text(face = "bold", size = 9,
                                      colour = "grey15", family = "Arial"),
    plot.subtitle      = element_text(size = 7.5, colour = "grey40",
                                      family = "Arial")
  )


# =============================================================================
# 0. LOAD DATA
# =============================================================================

raw_data  <- read.csv("Fig_4c_Raw_expression_data.csv",    check.names = FALSE)
ins1_norm <- read.csv("Fig_4c_Ins1_normalized_data.csv",   check.names = FALSE)
kcnj_norm <- read.csv("Fig_4c_Kcnj11_normalized_data.csv", check.names = FALSE)
hes1_norm <- read.csv("Fig_4c_Hes1_normalized_data.csv",   check.names = FALSE)

drop_ryr1 <- function(df) df[, colnames(df) != "Ryr1", drop = FALSE]
raw_data  <- drop_ryr1(raw_data)
ins1_norm <- drop_ryr1(ins1_norm)
kcnj_norm <- drop_ryr1(kcnj_norm)
hes1_norm <- drop_ryr1(hes1_norm)

consensus_df  <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
consensus_ids <- consensus_df[[1]]

n_obs <- nrow(raw_data)   # 12


# =============================================================================
# 1. HELPERS
# =============================================================================

compute_pairs <- function(df, method_label) {
  df_mat  <- as.matrix(na.omit(df))
  genes   <- colnames(df_mat)
  n_genes <- length(genes)
  out     <- vector("list", choose(n_genes, 2))
  k <- 0L
  for (i in seq_len(n_genes - 1)) {
    for (j in seq(i + 1, n_genes)) {
      ct  <- cor.test(df_mat[, i], df_mat[, j], method = "pearson")
      k   <- k + 1L
      out[[k]] <- data.frame(
        Gene1   = genes[i], Gene2 = genes[j],
        pair_id = paste(genes[i], genes[j], sep = "_"),
        r       = as.numeric(ct$estimate),
        p_val   = ct$p.value,
        CI_low  = ct$conf.int[1],
        CI_high = ct$conf.int[2],
        Method  = method_label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

fisher_pool <- function(r_vec, n_per_method = 12, alpha = 0.05) {
  r_vec  <- pmin(pmax(r_vec, -0.9999), 0.9999)
  zs     <- atanh(r_vec)
  ws     <- rep(n_per_method - 3, length(r_vec))
  z_bar  <- weighted.mean(zs, ws)
  se     <- 1 / sqrt(sum(ws))
  z_crit <- qnorm(1 - alpha / 2)
  list(r_pool = tanh(z_bar),
       ci_lo  = tanh(z_bar - z_crit * se),
       ci_hi  = tanh(z_bar + z_crit * se))
}


# =============================================================================
# 2. COMPUTE CORRELATIONS
# =============================================================================

method_levels <- c("Raw\nexpression", "Ins1-\nnormalized",
                   "Kcnj11-\nnormalized", "Hes1-\nnormalized")

all_pairs <- bind_rows(
  compute_pairs(raw_data,  "Raw\nexpression"),
  compute_pairs(ins1_norm, "Ins1-\nnormalized"),
  compute_pairs(kcnj_norm, "Kcnj11-\nnormalized"),
  compute_pairs(hes1_norm, "Hes1-\nnormalized")
) %>%
  mutate(
    Method    = factor(Method, levels = method_levels),
    sig       = p_val < 0.05,
    consensus = pair_id %in% consensus_ids
  )


# =============================================================================
# 3. PERMUTATION NULL (1,000 permutations)
# =============================================================================

set.seed(2025)
run_permutation <- function(df_expr, n_perm = 1000, alpha = 0.05) {
  df_mat  <- as.matrix(na.omit(df_expr))
  n       <- nrow(df_mat)
  n_genes <- ncol(df_mat)
  null_counts <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    perm_mat <- apply(df_mat, 2, function(x) x[sample(n)])
    cnt <- 0L
    for (i in seq_len(n_genes - 1))
      for (j in seq(i + 1, n_genes)) {
        ct <- cor.test(perm_mat[, i], perm_mat[, j], method = "pearson")
        if (!is.na(ct$p.value) && ct$p.value < alpha) cnt <- cnt + 1L
      }
    null_counts[p] <- cnt
  }
  list(mean = mean(null_counts), sd = sd(null_counts),
       q95  = quantile(null_counts, 0.95))
}

cat("Running permutations...\n")
nulls <- list(
  "Raw\nexpression"     = run_permutation(raw_data),
  "Ins1-\nnormalized"   = run_permutation(ins1_norm),
  "Kcnj11-\nnormalized" = run_permutation(kcnj_norm),
  "Hes1-\nnormalized"   = run_permutation(hes1_norm)
)


# =============================================================================
# 4. PANEL a — Observed vs Permutation Null
# =============================================================================

obs_df <- all_pairs %>%
  filter(sig) %>%
  count(Method, name = "Observed") %>%
  mutate(
    Null_Mean  = sapply(as.character(Method), function(m) nulls[[m]]$mean),
    Null_SD    = sapply(as.character(Method), function(m) nulls[[m]]$sd),
    Emp_p      = pnorm(Observed, mean = Null_Mean, sd = Null_SD,
                       lower.tail = FALSE),
    sig_label  = case_when(
      Emp_p < 0.001 ~ "***",
      Emp_p < 0.01  ~ "**",
      Emp_p < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    Fill_group = "Single",
    label_y    = pmax(Observed, Null_Mean + Null_SD) + max(Observed) * 0.04
  )

# Add consensus row
consensus_row <- tibble(
  Method     = factor("4-Method\nConsensus",
                      levels = c(method_levels, "4-Method\nConsensus")),
  Observed   = length(consensus_ids),
  Null_Mean  = 0.05^4 * choose(ncol(raw_data), 2),
  Null_SD    = NA_real_,
  Emp_p      = NA_real_,
  sig_label  = "",
  Fill_group = "Consensus",
  label_y    = length(consensus_ids) + max(obs_df$Observed) * 0.04
)

summary_df <- bind_rows(obs_df, consensus_row) %>%
  mutate(Method = factor(Method,
                         levels = c(method_levels, "4-Method\nConsensus")))

fill_colors <- c("Single" = npg_colors[4], "Consensus" = npg_colors[1])

p_s1a <- ggplot(summary_df,
                aes(x = Method, y = Observed, fill = Fill_group)) +
  
  geom_col(width = 0.62, alpha = 0.7, color = NA) +
  
  geom_errorbar(
    data = filter(summary_df, !is.na(Null_SD)),
    aes(x = Method,
        ymin = Null_Mean - Null_SD,
        ymax = Null_Mean + Null_SD),
    width = 0.22, linewidth = 0.7, color = "grey25",
    inherit.aes = FALSE
  ) +
  
  geom_point(
    data = filter(summary_df, !is.na(Null_SD)),
    aes(x = Method, y = Null_Mean),
    shape = 18, size = 2.8, color = "grey25",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    aes(y = Observed * 0.92, label = Observed),
    size = 3.0, family = "Arial", fontface = "bold",
    color = "white", vjust = 1
  ) +
  
  geom_text(
    aes(y = label_y, label = sig_label),
    size = 3.5, family = "Arial", color = "grey20"
  ) +
  
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title    = "(a)  Observed vs Permutation Null",
    subtitle = "\u25c6 = null mean \u00b1 SD (1,000 permutations, \u03b1 = 0.05)",
    x = NULL, y = "Significant Gene Pairs"
  ) +
  theme_Fig_S1 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))


# =============================================================================
# 5. PANEL b — Fisher-z pooled CI: single vs consensus (sign-stable pairs)
# =============================================================================

hey1_pairs <- c("Hey1_Cacna1c", "Tcf7l2_Hey1", "Hey1_Cacna1d",
                "Wnt5b_Hey1",   "Hey1_Ptbp1",  "Wnt5a_Hey1")

# Sign-stable consensus pairs (21)
ci_df <- all_pairs %>%
  filter(consensus, !pair_id %in% hey1_pairs) %>%
  group_by(pair_id) %>%
  summarise(
    pool    = list(fisher_pool(r, n_per_method = n_obs)),
    r_pool  = pool[[1]]$r_pool,
    lo_pool = pool[[1]]$ci_lo,
    hi_pool = pool[[1]]$ci_hi,
    r_raw   = r[Method == "Raw\nexpression"],
    lo_raw  = CI_low[Method  == "Raw\nexpression"],
    hi_raw  = CI_high[Method == "Raw\nexpression"],
    .groups = "drop"
  ) %>%
  mutate(
    width_single     = hi_raw  - lo_raw,
    width_pool       = hi_pool - lo_pool,
    CI_reduction_pct = (width_single - width_pool) / width_single * 100,
    pair_label       = gsub("_", "\u2013", pair_id),
    pair_label       = reorder(pair_label, r_pool)
  ) %>%
  select(-pool)

mean_reduction <- mean(ci_df$CI_reduction_pct, na.rm = TRUE)

ci_long <- ci_df %>%
  select(pair_label, width_single, width_pool) %>%
  pivot_longer(c(width_single, width_pool),
               names_to = "Approach", values_to = "CI_width") %>%
  mutate(Approach = recode(Approach,
                           width_single = "Single method (Raw)",
                           width_pool   = "4-method pooled"))

p_s1b <- ggplot(ci_long,
                aes(x = pair_label, y = CI_width,
                    colour = Approach, shape = Approach)) +
  geom_line(aes(group = pair_label), colour = "grey75", linewidth = 0.35) +
  geom_point(size = 2.0, alpha = 0.88) +
  scale_colour_manual(
    values = c("Single method (Raw)" = npg_colors[1],
               "4-method pooled"     = npg_colors[4])
  ) +
  scale_shape_manual(
    values = c("Single method (Raw)" = 16,
               "4-method pooled"     = 18)
  ) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.08))) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Mean CI reduction: %.1f%%", mean_reduction),
           hjust = 1.05, vjust = 1.6,
           size = 2.8, family = "Arial",
           colour = npg_colors[4], fontface = "bold") +
  coord_flip() +
  labs(
    title    = "(b)  95% CI Width: Single vs 4-Method Pooled",
    subtitle = "21 sign-stable consensus pairs | Fisher z meta-analytic pooling",
    x = NULL, y = "95% CI width",
    colour = NULL, shape = NULL
  ) +
  theme_Fig_S1 +
  theme(
    axis.text.y     = element_text(size = 6.5, face = "italic"),
    legend.position = "bottom"
  )


# =============================================================================
# 6. PANEL c — Hey1 sign-flip: self-correcting property of consensus filter
# =============================================================================

flip_df <- all_pairs %>%
  filter(consensus, pair_id %in% hey1_pairs) %>%
  group_by(pair_id) %>%
  summarise(
    pool    = list(fisher_pool(r, n_per_method = n_obs)),
    r_pool  = pool[[1]]$r_pool,
    lo_pool = pool[[1]]$ci_lo,
    hi_pool = pool[[1]]$ci_hi,
    .groups = "drop"
  ) %>%
  select(-pool) %>%
  left_join(
    all_pairs %>%
      filter(pair_id %in% hey1_pairs) %>%
      select(pair_id, Method, r) %>%
      pivot_wider(names_from = Method, values_from = r,
                  names_glue = "r_{Method}"),
    by = "pair_id"
  ) %>%
  mutate(
    pair_label = gsub("_", "\u2013", pair_id),
    pair_label = reorder(pair_label, r_pool)
  )

# Reshape for dot plot: per-method r values + pooled
flip_long <- all_pairs %>%
  filter(pair_id %in% hey1_pairs) %>%
  mutate(pair_label = gsub("_", "\u2013", pair_id)) %>%
  bind_rows(
    flip_df %>%
      transmute(pair_id, pair_label,
                Method = factor("4-Method\nPooled"),
                r = r_pool,
                CI_low  = lo_pool,
                CI_high = hi_pool)
  ) %>%
  mutate(
    Method     = factor(Method,
                        levels = c(method_levels, "4-Method\nPooled")),
    pair_label = reorder(pair_label, ifelse(Method == "Raw\nexpression", r, NA),
                         FUN = function(x) mean(x, na.rm = TRUE))
  )

method_dot_colors <- c(
  "Raw\nexpression"     = npg_colors[7],
  "Ins1-\nnormalized"   = npg_colors[3],
  "Kcnj11-\nnormalized" = npg_colors[2],
  "Hes1-\nnormalized"   = npg_colors[4],
  "4-Method\nPooled"    = "grey30"
)

p_s1c <- ggplot(flip_long,
                aes(x = pair_label, y = r,
                    colour = Method, shape = Method)) +
  
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  
  # ── FIXED: single aes() with x, ymin, ymax ──────────────────────────────
  geom_errorbar(
    data = filter(flip_long, Method == "4-Method\nPooled"),
    aes(x = pair_label, ymin = CI_low, ymax = CI_high),
    width = 0.25, linewidth = 0.6, colour = "grey30",
    inherit.aes = FALSE
  ) +
  
  geom_point(size = 2.2, alpha = 0.90) +
  
  scale_colour_manual(values = method_dot_colors) +
  scale_shape_manual(
    values = c(16, 16, 16, 16, 18)
  ) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.5)) +
  coord_flip() +
  labs(
    title    = "(c)  Hey1 Pairs: Sign Reversal Across Normalization Methods",
    subtitle = "Pooled CI crosses zero \u2192 correctly excluded by consensus filter",
    x = NULL, y = "Pearson r",
    colour = NULL, shape = NULL
  ) +
  theme_Fig_S1 +
  theme(
    axis.text.y     = element_text(size = 7.5, face = "italic"),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm")
  )


# =============================================================================
# 7. CSV EXPORTS
# =============================================================================

write.csv(
  summary_df %>%
    mutate(Method = gsub("\n", " ", as.character(Method))) %>%
    select(Method, Observed, Null_Mean, Null_SD, Emp_p, sig_label),
  "Fig_S1a_Observed_vs_Null.csv", row.names = FALSE, na = ""
)

write.csv(
  ci_df %>%
    mutate(pair_label = as.character(pair_label)) %>%
    select(pair_label, r_pool, lo_pool, hi_pool,
           r_raw, lo_raw, hi_raw,
           width_single, width_pool, CI_reduction_pct) %>%
    arrange(desc(r_pool)),
  "Fig_S1b_CI_width_pooled_vs_single.csv", row.names = FALSE, na = ""
)

write.csv(
  flip_df %>%
    mutate(pair_label = as.character(pair_label)) %>%
    select(pair_label, r_pool, lo_pool, hi_pool) %>%
    left_join(
      all_pairs %>%
        filter(pair_id %in% hey1_pairs) %>%
        mutate(Method = gsub("\n", " ", as.character(Method))) %>%
        select(pair_id, Method, r) %>%
        pivot_wider(names_from = Method, values_from = r),
      by = c("pair_label" = "pair_id")
    ),
  "Fig_S1c_Hey1_signflip.csv", row.names = FALSE, na = ""
)

cat("CSVs exported.\n")


# =============================================================================
# 8. ASSEMBLE — 3 rows, 15 cm × 20 cm
# =============================================================================

p_Fig_S1 <- p_s1a / p_s1b / p_s1c +
  plot_layout(heights = c(1, 1.2, 1)) +
  plot_annotation(
    #title    = "Fig. S1 \u2502 Multi-Normalization Consensus: Statistical Validation",
    subtitle = paste0(
      #"n = ", n_obs,
      #" (4 biological replicates \u00d7 3 conditions) \u2502 ",
      #length(consensus_ids),
      #" consensus pairs significant across all 4 normalization methods"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 10,
                                   family = "Arial", colour = "grey10"),
      plot.subtitle = element_text(size = 8, family = "Arial",
                                   colour = "grey40")
    )
  )


# =============================================================================
# 9. SAVE
# =============================================================================

ggsave("Fig_S1_Validation.svg",
       plot   = p_Fig_S1,
       width  = 15, height = 20, units = "cm",
       dpi    = 300)

ggsave("Fig_S1_Validation.png",
       plot   = p_Fig_S1,
       width  = 15, height = 20, units = "cm",
       dpi    = 300)

cat(sprintf(
  "\nFig. S1 saved (3 rows, 15 \u00d7 20 cm).\n  n=%d | consensus=%d | mean CI reduction=%.1f%%\n",
  n_obs, length(consensus_ids), mean_reduction
))

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 13: Fig_S1_validation_analysis.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 14/21  —  Fig_5abcd_network_analysis.r                                ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 5a–5d: Gene interaction network — 15 nodes, STRING classification, novel pair ranking
# Lines: 345  |  File: Fig_5abcd_network_analysis.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Gene Interaction Network Analysis - Figure 5

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(ggraph)
library(patchwork)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Helper function: create canonical pair names
create_pair_name <- function(gene1, gene2) {
  pairs <- mapply(function(g1, g2) {
    sorted <- sort(c(g1, g2))
    paste(sorted, collapse = "-")
  }, gene1, gene2, SIMPLIFY = TRUE)
  return(pairs)
}

# Load STRING interactions
string_data <- read.csv("string_interactions.csv", stringsAsFactors = FALSE) %>%
  mutate(Pair = create_pair_name(node1, node2)) %>%
  select(Pair, node1, node2, combined_score) %>%
  rename(STRING_Score = combined_score) %>%
  distinct(Pair, .keep_all = TRUE)

# Load consensus pairs from Fig 4d
validated_pairs_raw <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)
validated_pairs <- validated_pairs_raw %>%
  mutate(
    Gene1 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 1),
    Gene2 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 2),
    Pair = create_pair_name(Gene1, Gene2)
  ) %>%
  select(Pair, Gene1, Gene2) %>%
  distinct(Pair, .keep_all = TRUE)

# Load Hes1-normalized correlations from Fig 4e
correlation_data <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                             stringsAsFactors = FALSE, row.names = 1) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene1") %>%
  pivot_longer(-Gene1, names_to = "Gene2", values_to = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  filter(Gene1 < Gene2) %>%
  mutate(
    Pair = create_pair_name(Gene1, Gene2),
    r = Correlation
  ) %>%
  select(Pair, Gene1, Gene2, r)

# Classify pairs
string_set <- string_data$Pair
validated_set <- validated_pairs$Pair

all_pairs <- correlation_data %>%
  filter(abs(r) >= 0.7) %>%
  mutate(
    In_STRING = Pair %in% string_set,
    Four_Method_Validated = Pair %in% validated_set,
    Category = case_when(
      In_STRING & Four_Method_Validated ~ "Known Validated",
      !In_STRING & Four_Method_Validated ~ "Novel Discovery",
      In_STRING & !Four_Method_Validated ~ "STRING Only",
      TRUE ~ "Not Validated"
    )
  )

novel_discoveries <- all_pairs %>%
  filter(Category == "Novel Discovery") %>%
  arrange(desc(abs(r)))

# Summary statistics
summary_stats <- all_pairs %>%
  group_by(Category) %>%
  summarise(
    Count = n(),
    Mean_r = mean(abs(r)),
    Median_r = median(abs(r)),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

# Custom theme
theme_Fig_5 <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = "Arial"),
      axis.title = element_text(size = 10, face = "bold", color = "black", family = "Arial"),
      axis.text = element_text(size = 10, color = "black", family = "Arial"),
      legend.text = element_text(size = 10, family = "Arial"),
      legend.title = element_text(size = 10, face = "bold", family = "Arial"),
      panel.grid.major = element_line(color = "gray70", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# Panel A: Gene Interaction Network
string_validated <- all_pairs %>%
  filter(Category == "Known Validated") %>%
  arrange(desc(abs(r))) %>%
  head(30)

novel_top <- novel_discoveries %>% head(20)

remaining_needed <- 50 - nrow(string_validated) - nrow(novel_top)
if (remaining_needed > 0) {
  supplement <- all_pairs %>%
    filter(Category %in% c("STRING Only", "Not Validated")) %>%
    arrange(desc(abs(r))) %>%
    head(remaining_needed)
} else {
  supplement <- data.frame()
}

# Create edge list
edges_list <- list()
if (nrow(string_validated) > 0) {
  edges_list[[1]] <- data.frame(
    from = string_validated$Gene1,
    to = string_validated$Gene2,
    weight = abs(string_validated$r),
    EdgeType = "STRING Validated",
    stringsAsFactors = FALSE
  )
}
if (nrow(novel_top) > 0) {
  edges_list[[2]] <- data.frame(
    from = novel_top$Gene1,
    to = novel_top$Gene2,
    weight = abs(novel_top$r),
    EdgeType = "Novel Discovery",
    stringsAsFactors = FALSE
  )
}
if (nrow(supplement) > 0) {
  edges_list[[3]] <- data.frame(
    from = supplement$Gene1,
    to = supplement$Gene2,
    weight = abs(supplement$r),
    EdgeType = supplement$Category,
    stringsAsFactors = FALSE
  )
}

edges <- do.call(rbind, edges_list)

# Create graph
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, weights = E(g)$weight, directed = FALSE)

novel_genes <- unique(c(novel_top$Gene1, novel_top$Gene2))
V(g)$node_type <- ifelse(V(g)$name %in% novel_genes, "Novel", "Known")

edge_colors <- case_when(
  E(g)$EdgeType == "Novel Discovery" ~ npg_colors[1],
  E(g)$EdgeType == "STRING Validated" ~ npg_colors[4],
  E(g)$EdgeType == "STRING Only" ~ npg_colors[2],
  TRUE ~ "gray70"
)

E(g)$color <- edge_colors
E(g)$width <- E(g)$weight * 1.5

set.seed(42)
layout_net <- layout_with_fr(g, weights = E(g)$weight)

panel_a <- ggraph(g, layout = layout_net) +
  geom_edge_link(aes(edge_color = EdgeType, edge_width = weight),
                 alpha = 0.6, show.legend = TRUE) +
  geom_node_point(aes(size = degree, color = node_type, fill = node_type),
                  alpha = 0.8, shape = 21) +
  geom_node_text(aes(label = name, color = node_type),
                 size = 2.5, repel = TRUE, max.overlaps = 20,
                 fontface = "italic", family = "Arial") +
  scale_edge_color_manual(
    name = "Interaction",
    values = c("STRING Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  scale_edge_width_continuous(name = "Abs. Correlation (|r|)", range = c(0.3, 2)) +
  scale_color_manual(
    name = "Gene Status",
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_fill_manual(
    name = "Gene Status",
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 10)) +
  labs(title = "Gene Interaction Network") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, face = "bold", family = "Arial")
  )

ggsave("Fig_5a_Network.png",
       plot = panel_a, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("Fig_5a_Network.svg",
       plot = panel_a, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

# Panel B: STRING Coverage Distribution
panel_b_data <- all_pairs %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(
    Percentage = 100 * Count / sum(Count),
    Category = factor(Category,
                      levels = c("Known Validated", "Novel Discovery",
                                 "STRING Only", "Not Validated"))
  )

panel_b <- ggplot(panel_b_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", color = "white", family = "Arial") +
  scale_fill_manual(
    name = NULL,
    values = c("Known Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  labs(title = "STRING Database Coverage") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "bottom",
    legend.text = element_text(size = 8, family = "Arial")
  )

ggsave("Fig_5b_Coverage.png",
       plot = panel_b, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("Fig_5b_Coverage.svg",
       plot = panel_b, scale = 1.5, width = 8, height = 7.5,
       dpi = 600, units = "cm", bg = "white")

write.csv(panel_b_data,
          "Fig_5b_Coverage_Summary.csv",
          row.names = FALSE)

# Panel C: Correlation Distribution
panel_c <- ggplot(all_pairs, aes(x = abs(r), fill = Category)) +
  geom_histogram(bins = 15, alpha = 0.7, position = "identity",
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black",
             linewidth = 0.4, alpha = 0.5) +
  facet_wrap(~factor(Category,
                     levels = c("Known Validated", "Novel Discovery",
                                "STRING Only", "Not Validated")),
             ncol = 2, scales = "free_y") +
  scale_fill_manual(
    values = c("Known Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  labs(title = "Correlation Strength Distribution",
       x = "Absolute Correlation (|r|)",
       y = "Frequency") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    axis.title = element_text(size = 10, family = "Arial", color = "black")
  )

ggsave("Fig_5c_Distribution.png",
       plot = panel_c, scale = 1.5, width = 8, height = 8.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("Fig_5c_Distribution.svg",
       plot = panel_c, scale = 1.5, width = 8, height = 8.5,
       dpi = 600, units = "cm", bg = "white")

# Panel D: Top Novel Discoveries
top_novel <- novel_discoveries %>%
  head(15) %>%
  mutate(
    Pair_Label = paste(Gene1, Gene2, sep = "-"),
    Rank = row_number()
  ) %>%
  arrange(desc(abs(r)))

panel_d <- ggplot(top_novel, aes(x = reorder(Pair_Label, abs(r)), y = abs(r))) +
  geom_col(fill = npg_colors[1], alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", r)),
            hjust = -0.2, size = 2.5, family = "Arial") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Top 15 Novel Gene Pair Discoveries",
       x = NULL,
       y = "Absolute Correlation (|r|)") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    axis.text.y = element_text(size = 8, face = "italic", family = "Arial"),
    axis.text.x = element_text(size = 10, family = "Arial"),
    axis.title = element_text(size = 10, family = "Arial"),
    panel.grid.major.y = element_blank()
  )

ggsave("Fig_5d_TopNovel.png",
       plot = panel_d, scale = 1.5, width = 8, height = 8.5,
       dpi = 600, units = "cm", bg = "white")

ggsave("Fig_5d_TopNovel.svg",
       plot = panel_d, scale = 1.5, width = 8, height = 8.5,
       dpi = 600, units = "cm", bg = "white")

write.csv(top_novel,
          "Fig_5d_Top_Novel_Discoveries.csv",
          row.names = FALSE)

# Export summary
write.csv(summary_stats,
          "Fig_5_Summary_Statistics.csv",
          row.names = FALSE)

write.csv(all_pairs,
          "Fig_5_All_Classified_Pairs.csv",
          row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 14: Fig_5abcd_network_analysis.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 15/21  —  Fig_6b_consensus_hub_network.r                              ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 6b: Consensus hub gene network (circular layout, igraph/ggraph)
# Lines: 164  |  File: Fig_6b_consensus_hub_network.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Consensus Hub Gene Network Analysis

rm(list = ls())
set.seed(42)

library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load data
consensus_pairs <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE) %>%
  separate(Consensus_Pairs, into = c("Gene1", "Gene2"), sep = "_", remove = FALSE)

corr_mat_raw <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)
gene_names <- corr_mat_raw[[1]]
corr_mat <- corr_mat_raw[, -1]
rownames(corr_mat) <- gene_names

# Convert correlation matrix to edge list
corr_long <- as.data.frame(as.matrix(corr_mat)) %>%
  rownames_to_column(var = "Gene1") %>%
  pivot_longer(cols = -Gene1, names_to = "Gene2", values_to = "r") %>%
  filter(Gene1 != Gene2, !is.na(r)) %>%
  mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  distinct(pair_key, .keep_all = TRUE)

# Merge with consensus pairs
consensus_pairs <- consensus_pairs %>%
  mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  left_join(corr_long %>% select(pair_key, r), by = "pair_key") %>%
  mutate(edge_type = ifelse(r > 0, "Positive", "Negative"))

# Gene categories
calcium_signaling <- c("Cacna1c", "Cacna1d", "Itpr1", "Ryr1", "Ryr2", "Ryr3")
incretin <- c("Glp1r")
wnt_signaling <- c("Hey1", "Tcf7", "Tcf7l2", "Wnt5a", "Wnt5b", "Wnt9a")
metabolic <- c("Glut2", "Ptbp1")

all_genes <- unique(c(consensus_pairs$Gene1, consensus_pairs$Gene2))

gene_metadata <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
  mutate(Category = case_when(
    Gene %in% calcium_signaling ~ "Calcium Signaling",
    Gene %in% incretin ~ "Incretin Signaling",
    Gene %in% wnt_signaling ~ "Wnt Signaling",
    Gene %in% metabolic ~ "Metabolic Regulation",
    TRUE ~ "Other"
  ))

# Build network
g_consensus <- graph_from_data_frame(
  d = consensus_pairs %>% select(Gene1, Gene2, r, edge_type),
  vertices = gene_metadata,
  directed = FALSE
)

E(g_consensus)$correlation <- consensus_pairs$r
E(g_consensus)$abs_correlation <- abs(consensus_pairs$r)
E(g_consensus)$edge_type <- consensus_pairs$edge_type

V(g_consensus)$degree <- degree(g_consensus)
V(g_consensus)$betweenness <- betweenness(g_consensus, normalized = TRUE)
V(g_consensus)$closeness <- closeness(g_consensus, normalized = TRUE)
V(g_consensus)$eigen_centrality <- eigen_centrality(g_consensus)$vector

hub_threshold <- 6
V(g_consensus)$is_hub <- ifelse(V(g_consensus)$degree >= hub_threshold, "Hub", "Non-hub")

# Edge width scaling
min_corr <- min(E(g_consensus)$abs_correlation)
max_corr <- max(E(g_consensus)$abs_correlation)
E(g_consensus)$width <- 0.5 + (E(g_consensus)$abs_correlation - min_corr) / 
                         (max_corr - min_corr) * (3.0 - 0.5)

# NPG colors
npg_category_colors <- c(
  "Calcium Signaling" = npg_colors[4],
  "Incretin Signaling" = npg_colors[1],
  "Wnt Signaling" = npg_colors[2],
  "Metabolic Regulation" = npg_colors[3]
)

# Circular network plot
p_network <- ggraph(g_consensus, layout = "circle") +
  geom_edge_link(
    aes(width = width, linetype = edge_type),
    color = "grey70",
    alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(
    name = "Correlation",
    values = c("Positive" = "solid", "Negative" = "dotted"),
    guide = guide_legend(order = 3, override.aes = list(color = "black", width = 1))
  ) +
  scale_edge_width_continuous(
    name = "Correlation Strength",
    range = c(0.5, 3.0),
    breaks = c(0.5, 1.5, 2.5, 3.0),
    labels = c("Weak", "Moderate", "Strong", "Very Strong"),
    guide = guide_legend(order = 4, override.aes = list(color = "black"))
  ) +
  geom_node_point(
    aes(size = degree, color = Category, fill = Category),
    shape = 21,
    stroke = 1.2,
    alpha = 0.85
  ) +
  scale_color_manual(name = "Functional Category", values = npg_category_colors, 
                     guide = guide_legend(order = 1)) +
  scale_fill_manual(name = "Functional Category", values = npg_category_colors, 
                    guide = "none") +
  geom_node_text(
    aes(label = name, fontface = ifelse(degree >= hub_threshold, "bold.italic", "italic")),
    repel = TRUE, size = 5, max.overlaps = 25, family = "Arial"
  ) +
  scale_size_continuous(
    name = "Degree Centrality",
    range = c(5, 14),
    breaks = c(2, 4, 6),
    labels = c("Low (2)", "Medium (4)", "High (6+)"),
    guide = guide_legend(order = 2)
  ) +
  theme_void(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, family = "Arial"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.key.size = unit(0.4, "cm")
  ) +
  labs(title = "Consensus Hub Gene Network")

ggsave("Fig_6b_consensus_hub_network.png", p_network, scale = 2.3,
       width = 7.7, height = 6, units = "cm", dpi = 600, bg = "white")
ggsave("Fig_6b_consensus_hub_network.svg", p_network, scale = 2.3,
       width = 7.7, height = 6, units = "cm", dpi = 600, bg = "white")

# Export statistics
gene_stats <- data.frame(
  Gene = V(g_consensus)$name,
  Degree = V(g_consensus)$degree,
  Betweenness = round(V(g_consensus)$betweenness, 4),
  Closeness = round(V(g_consensus)$closeness, 4),
  Eigenvector = round(V(g_consensus)$eigen_centrality, 4),
  Category = V(g_consensus)$Category,
  Hub_Status = V(g_consensus)$is_hub
) %>% arrange(desc(Degree))

write.csv(gene_stats, "Fig_6b_gene_centrality_summary.csv", row.names = FALSE)

edge_stats <- igraph::as_data_frame(g_consensus, what = "edges") %>%
  arrange(desc(abs(correlation))) %>%
  select(from, to, correlation, abs_correlation, edge_type, width) %>%
  rename(Gene1 = from, Gene2 = to, `Correlation (r)` = correlation,
         `|Correlation|` = abs_correlation, `Edge Type` = edge_type,
         `Edge Width (pt)` = width)

write.csv(edge_stats, "Fig_6b_edge_correlation_statistics.csv", row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 15: Fig_6b_consensus_hub_network.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 16/21  —  Fig_6c_MLSS_drug_combinations.r                             ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 6c: MLSS v4.0 scoring — all drug combinations ranked by synergy score
# Lines: 340  |  File: Fig_6c_MLSS_drug_combinations.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# MLSS v4.0 — Multi-Layer Synergy Score for Drug Combinations
#
# Formula:
#   MLSS = (0.45·C + 0.30·B + 0.10·V + 0.15·P) × 9
#          + Σ(r_pos × potency_mult)
#          - Σ(|r_neg| × potency_mult)
#
# potency_mult = min(max(pIC50_i, pIC50_j) / 9, 1)
#
# Inputs:
#   Fig_6b_Gene_Centrality_Summary.csv
#   Fig_4e_Hes1_significant_correlations_only.csv
#   drug.target.interaction.csv
#
# v4.0 — January 2026

rm(list = ls())
set.seed(42)

# ---- packages ----------------------------------------------------------------

required_packages <- c("dplyr", "tidyr", "ggplot2", "ggsci", "tibble",
                       "igraph", "ggraph", "scales")

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}))

npg_colors <- pal_npg("nrc")(10)

# ---- load data ---------------------------------------------------------------

hub_data     <- read.csv("Fig_6b_Gene_Centrality_Summary.csv",         stringsAsFactors = FALSE)
drug_data    <- read.csv("drug.target.interaction.csv",                 stringsAsFactors = FALSE)
corr_mat_raw <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                         stringsAsFactors = FALSE, check.names = FALSE)

gene_names         <- corr_mat_raw[[1]]
corr_mat           <- corr_mat_raw[, -1]
rownames(corr_mat) <- toupper(gene_names)
colnames(corr_mat) <- toupper(gene_names)

# ---- potency lookup (best pIC50 per drug-gene pair) -------------------------

potency_lookup <- drug_data %>%
  dplyr::filter(ACT_TYPE %in% c("IC50", "Ki", "EC50"),
                ORGANISM  == "Homo sapiens",
                RELATION  == "=",
                !is.na(ACT_VALUE)) %>%
  dplyr::mutate(Gene_Standardized = toupper(GENE)) %>%
  dplyr::group_by(DRUG_NAME, Gene_Standardized) %>%
  dplyr::summarise(pIC50 = max(ACT_VALUE, na.rm = TRUE), .groups = "drop")

# ---- hub genes (degree >= 6) -------------------------------------------------

hub_genes_df <- hub_data %>% dplyr::filter(Degree >= 6)
hub_genes    <- unique(toupper(hub_genes_df$Gene))
hub_data_cat <- hub_genes_df %>% dplyr::mutate(Gene_Upper = toupper(Gene))

cat(sprintf("Hub genes (n=%d): %s\n", length(hub_genes), paste(hub_genes, collapse = ", ")))

# ---- pairwise correlation table ----------------------------------------------

correlations_slim <- as.data.frame(as.matrix(corr_mat)) %>%
  tibble::rownames_to_column("Gene1") %>%
  tidyr::pivot_longer(-Gene1, names_to = "Gene2", values_to = "r") %>%
  dplyr::filter(Gene1 != Gene2, !is.na(r)) %>%
  dplyr::mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  dplyr::distinct(pair_key, .keep_all = TRUE) %>%
  dplyr::select(pair_key, r)

# ---- drug-target table (hub genes only) -------------------------------------

drug_targets <- drug_data %>%
  dplyr::filter(toupper(GENE) %in% hub_genes, ORGANISM == "Homo sapiens") %>%
  dplyr::select(DRUG_NAME, GENE, TARGET_NAME, ACTION_TYPE) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Gene_Standardized = toupper(GENE))

unique_drugs <- sort(unique(drug_targets$DRUG_NAME))
cat(sprintf("Drugs targeting hub genes: %d\n", length(unique_drugs)))

# ---- helper functions --------------------------------------------------------

get_targets <- function(drug) {
  drug_targets %>% dplyr::filter(DRUG_NAME == drug) %>%
    dplyr::pull(Gene_Standardized) %>% unique()
}

get_pIC50 <- function(drug, gene) {
  v <- potency_lookup %>%
    dplyr::filter(DRUG_NAME == drug, Gene_Standardized == toupper(gene)) %>%
    dplyr::pull(pIC50)
  if (length(v) == 0) NA else v
}

potency_mult <- function(drug_a, gene_i, drug_b, gene_j) {
  p_i <- get_pIC50(drug_a, gene_i)
  p_j <- get_pIC50(drug_b, gene_j)
  if (is.na(p_i) || is.na(p_j)) return(0)
  min(max(p_i, p_j) / 9, 1.0)
}

mean_pIC50_drug <- function(drug, targets) {
  vals <- sapply(unique(toupper(targets)), function(g) get_pIC50(drug, g))
  m    <- mean(vals[!is.na(vals)], na.rm = TRUE)
  if (is.nan(m)) 0 else m
}

corr_stats <- function(genes_A, genes_B, corr_df) {
  genes_A <- unique(toupper(genes_A))
  genes_B <- unique(toupper(genes_B))
  empty   <- data.frame(Mean_Pos_r = 0, Max_Pos_r = 0, Mean_Neg_r = 0,
                        Sum_Neg_Abs_r = 0, Mean_Abs_r = 0)
  if (!length(genes_A) || !length(genes_B)) return(empty)
  
  pairs <- expand.grid(Gene1 = genes_A, Gene2 = genes_B, stringsAsFactors = FALSE) %>%
    dplyr::mutate(pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
    dplyr::left_join(corr_df, by = "pair_key") %>%
    dplyr::filter(!is.na(r))
  
  if (!nrow(pairs)) return(empty)
  
  pos <- pairs$r[pairs$r > 0]
  neg <- pairs$r[pairs$r < 0]
  data.frame(
    Mean_Pos_r    = if (length(pos)) mean(pos) else 0,
    Max_Pos_r     = if (length(pos)) max(pos)  else 0,
    Mean_Neg_r    = if (length(neg)) mean(neg) else 0,
    Sum_Neg_Abs_r = if (length(neg)) sum(abs(neg)) else 0,
    Mean_Abs_r    = mean(abs(pairs$r))
  )
}

synergy_score <- function(drug_a, ta, drug_b, tb, corr_df) {
  ta <- unique(toupper(ta)); tb <- unique(toupper(tb)); total <- 0
  for (i in seq_along(ta)) for (j in seq_along(tb)) {
    pk <- paste(pmin(ta[i], tb[j]), pmax(ta[i], tb[j]), sep = "_")
    r  <- corr_df %>% dplyr::filter(pair_key == pk) %>% dplyr::pull(r)
    if (!length(r) || is.na(r) || r <= 0) next
    total <- total + r * potency_mult(drug_a, ta[i], drug_b, tb[j])
  }
  total
}

antagonism_score <- function(drug_a, ta, drug_b, tb, corr_df) {
  ta <- unique(toupper(ta)); tb <- unique(toupper(tb)); total <- 0
  for (i in seq_along(ta)) for (j in seq_along(tb)) {
    pk <- paste(pmin(ta[i], tb[j]), pmax(ta[i], tb[j]), sep = "_")
    r  <- corr_df %>% dplyr::filter(pair_key == pk) %>% dplyr::pull(r)
    if (!length(r) || is.na(r) || r >= 0) next
    total <- total + abs(r) * potency_mult(drug_a, ta[i], drug_b, tb[j])
  }
  total
}

# ---- all drug combinations ---------------------------------------------------

all_combos <- expand.grid(Drug_A = unique_drugs, Drug_B = unique_drugs,
                          stringsAsFactors = FALSE) %>%
  dplyr::filter(Drug_A < Drug_B) %>%
  dplyr::mutate(Combination = paste(Drug_A, "+", Drug_B),
                Combo_ID    = dplyr::row_number())

total_combos <- nrow(all_combos)
cat(sprintf("Combinations to score: %d\n", total_combos))

# ---- enrich: targets, pathways -----------------------------------------------

all_combos <- all_combos %>%
  dplyr::mutate(
    Targets_A  = sapply(Drug_A, function(d) paste(get_targets(d), collapse = "|")),
    Targets_B  = sapply(Drug_B, function(d) paste(get_targets(d), collapse = "|")),
    n_targets_A = sapply(Drug_A, function(d) length(get_targets(d))),
    n_targets_B = sapply(Drug_B, function(d) length(get_targets(d))),
    n_targets_combined = n_targets_A + n_targets_B,
    Pathway_A  = sapply(Drug_A, function(d) {
      g <- get_targets(d)
      paste(unique(hub_data_cat$Category[hub_data_cat$Gene_Upper %in% g]), collapse = "|")
    }),
    Pathway_B  = sapply(Drug_B, function(d) {
      g <- get_targets(d)
      paste(unique(hub_data_cat$Category[hub_data_cat$Gene_Upper %in% g]), collapse = "|")
    }),
    Cross_Pathway = ifelse(Pathway_A != Pathway_B, "Yes", "No")
  )

# ---- complementarity (C), balance (B), network coverage (V) -----------------

all_combos <- all_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    targets_A_set  = list(if (nzchar(Targets_A)) unlist(strsplit(Targets_A, "\\|")) else character(0)),
    targets_B_set  = list(if (nzchar(Targets_B)) unlist(strsplit(Targets_B, "\\|")) else character(0)),
    n_intersect    = length(intersect(targets_A_set, targets_B_set)),
    n_union        = length(union(targets_A_set, targets_B_set)),
    Complementarity_C = if (n_union > 0) 1 - n_intersect / n_union else 0,
    Balance_B         = if ((n_targets_A + n_targets_B) > 0)
      1 - abs(n_targets_A - n_targets_B) / (n_targets_A + n_targets_B)
    else 0,
    cs = list(corr_stats(targets_A_set, targets_B_set, correlations_slim))
  ) %>%
  tidyr::unnest_wider(cs) %>%
  dplyr::ungroup()

# min-max normalize Mean_Abs_r → Network_Coverage_V
r_min <- min(all_combos$Mean_Abs_r)
r_max <- max(all_combos$Mean_Abs_r)

all_combos <- all_combos %>%
  dplyr::mutate(
    Network_Coverage_V = if (r_max > r_min) (Mean_Abs_r - r_min) / (r_max - r_min) else 0
  )

# ---- potency component (P) ---------------------------------------------------

all_combos <- all_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Mean_pIC50_A = mean_pIC50_drug(Drug_A, targets_A_set),
    Mean_pIC50_B = mean_pIC50_drug(Drug_B, targets_B_set),
    Mean_pIC50   = mean(c(Mean_pIC50_A, Mean_pIC50_B), na.rm = TRUE),
    Potency_P    = min(max(Mean_pIC50 - 4, 0), 5) / 5
  ) %>%
  dplyr::ungroup()

# ---- synergy bonus and antagonism penalty ------------------------------------

all_combos <- all_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Synergy_Bonus      = synergy_score(Drug_A, targets_A_set,
                                       Drug_B, targets_B_set, correlations_slim),
    Antagonism_Penalty = antagonism_score(Drug_A, targets_A_set,
                                          Drug_B, targets_B_set, correlations_slim)
  ) %>%
  dplyr::ungroup()

# ---- MLSS v4.0 ---------------------------------------------------------------

all_combos <- all_combos %>%
  dplyr::mutate(
    Base_Score = (0.45 * Complementarity_C +
                    0.30 * Balance_B          +
                    0.10 * Network_Coverage_V  +
                    0.15 * Potency_P) * 9,
    MLSS_v4.0  = Base_Score + Synergy_Bonus - Antagonism_Penalty
  ) %>%
  dplyr::arrange(dplyr::desc(MLSS_v4.0)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

top50 <- all_combos %>% head(50)
top15 <- all_combos %>% head(15)

cat("\nTop 15 combinations:\n")
print(top15 %>% dplyr::select(Rank, Combination, Complementarity_C, Balance_B,
                              Network_Coverage_V, Potency_P,
                              Synergy_Bonus, Antagonism_Penalty, MLSS_v4.0))

# ---- Fig 6c: component breakdown bar chart -----------------------------------

npg_mlss_colors <- c(
  "Complementarity (45%)" = npg_colors[1],
  "Balance (30%)"         = npg_colors[4],
  "Coverage (10%)"        = npg_colors[3],
  "Potency (15%)"         = npg_colors[2],
  "Synergy Bonus"         = npg_colors[5],
  "Antagonism Penalty"    = npg_colors[6]
)

plot_data <- top15 %>%
  dplyr::arrange(dplyr::desc(MLSS_v4.0)) %>%
  dplyr::mutate(
    `Complementarity (45%)` = Complementarity_C  * 0.45 * 9,
    `Balance (30%)`         = Balance_B           * 0.30 * 9,
    `Coverage (10%)`        = Network_Coverage_V  * 0.10 * 9,
    `Potency (15%)`         = Potency_P           * 0.15 * 9,
    `Synergy Bonus`         = Synergy_Bonus,
    `Antagonism Penalty`    = Antagonism_Penalty,
    Combination = factor(Combination, levels = rev(Combination))
  ) %>%
  dplyr::select(Combination, MLSS_v4.0,
                `Complementarity (45%)`, `Balance (30%)`, `Coverage (10%)`,
                `Potency (15%)`, `Synergy Bonus`, `Antagonism Penalty`) %>%
  tidyr::pivot_longer(-c(Combination, MLSS_v4.0), names_to = "Component", values_to = "Score") %>%
  dplyr::mutate(Component = factor(Component, levels = names(npg_mlss_colors)))

score_labels <- dplyr::distinct(plot_data, Combination, MLSS_v4.0)

p_6c <- ggplot(plot_data, aes(x = Combination, y = Score, fill = Component)) +
  geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.3) +
  geom_text(data = score_labels,
            aes(x = Combination, y = MLSS_v4.0, label = sprintf("%.2f", MLSS_v4.0)),
            hjust = -0.25, size = 3.5, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(name = "Score\nComponent", values = npg_mlss_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  guides(fill = guide_legend(nrow = 3)) +
  coord_flip() +
  labs(title = "MLSS Component Breakdown — Top 15 Drug Combinations",
       x = NULL, y = "Weighted Score Contribution") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(size = 14, face = "bold", hjust = 0.5,
                                      color = "black", family = "Arial"),
    axis.text          = element_text(size = 11, color = "black", family = "Arial"),
    axis.title.x       = element_text(size = 11, color = "black", family = "Arial"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 11, face = "bold", family = "Arial"),
    legend.text        = element_text(size = 11, family = "Arial"),
    legend.box         = "vertical"
  )

ggsave("Fig_6c_MLSS_top15.png", p_6c, width = 7.7, height = 6, units = "cm",
       scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_6c_MLSS_top15.svg", p_6c, width = 7.7, height = 6, units = "cm",
       scale = 2.3, dpi = 600, bg = "white")

# ---- export ------------------------------------------------------------------

export_cols <- c("Rank", "Combination", "Drug_A", "Drug_B",
                 "Pathway_A", "Pathway_B", "Cross_Pathway",
                 "n_targets_A", "n_targets_B", "n_targets_combined",
                 "Targets_A", "Targets_B",
                 "Mean_pIC50_A", "Mean_pIC50_B", "Mean_pIC50",
                 "Complementarity_C", "Balance_B", "Network_Coverage_V", "Potency_P",
                 "Mean_Pos_r", "Mean_Neg_r", "Mean_Abs_r",
                 "Base_Score", "Synergy_Bonus", "Antagonism_Penalty", "MLSS_v4.0")

write.csv(all_combos %>% dplyr::select(all_of(export_cols)),
          "Fig_6c_mlss_all_combinations.csv", row.names = FALSE)
write.csv(top15 %>% dplyr::select(all_of(export_cols)),
          "Fig_6c_mlss_top15.csv", row.names = FALSE)

cat(sprintf("Exported: %d combinations (all), top 15\n", nrow(all_combos)))
cat("Done. Output: Fig_6c_MLSS_top15 (.png/.svg)\n")

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 16: Fig_6c_MLSS_drug_combinations.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 17/21  —  Fig_6d_Fig_S2_MLSS_validation.r                             ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 6d + Fig S2: MLSS robustness — weight sensitivity, ablation analysis (8 scenarios)
# Lines: 222  |  File: Fig_6d_Fig_S2_MLSS_validation.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# MLSS v4.0 Validation - Sensitivity and Robustness Analysis

rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(reshape2)
library(gridExtra)

npg_colors <- pal_npg("nrc")(10)

# ── Load data ─────────────────────────────────────────────────────────────────
mlss_data <- read.csv("Fig_6c_mlss_all_combinations.csv", stringsAsFactors = FALSE)
top_15     <- mlss_data %>% head(15)

# ── Weight scenarios ──────────────────────────────────────────────────────────
weight_scenarios <- data.frame(
  Scenario = c("MLSS (Standard)", "C-Heavy (60:25:15)", "C-Heavy (55:20:25)",
               "B-Heavy (35:50:15)", "V-Heavy (35:30:35)", "P-Heavy (30:20:10:40)",
               "Balanced (40:40:10:10)", "Network-Focus (30:20:40:10)"),
  C_Weight = c(0.45, 0.60, 0.55, 0.35, 0.35, 0.30, 0.40, 0.30),
  B_Weight = c(0.30, 0.25, 0.20, 0.50, 0.30, 0.20, 0.40, 0.20),
  V_Weight = c(0.10, 0.10, 0.15, 0.10, 0.30, 0.10, 0.10, 0.40),
  P_Weight = c(0.15, 0.05, 0.10, 0.05, 0.05, 0.40, 0.10, 0.10),
  NPG_Color = npg_colors[1:8],
  stringsAsFactors = FALSE
)

# ── Recalculate scores for each scenario (top 15) ────────────────────────────
top_15_data    <- mlss_data %>% head(15)
scenario_results <- data.frame(Combination  = top_15_data$Combination,
                               Original_MLSS = top_15_data$MLSS_v4.0)

for (i in 1:nrow(weight_scenarios)) {
  scenario_name <- weight_scenarios$Scenario[i]
  c_w <- weight_scenarios$C_Weight[i]
  b_w <- weight_scenarios$B_Weight[i]
  v_w <- weight_scenarios$V_Weight[i]
  p_w <- weight_scenarios$P_Weight[i]
  
  new_base_score <- (c_w * top_15_data$Complementarity_C +
                       b_w * top_15_data$Balance_B +
                       v_w * top_15_data$Network_Coverage_V +
                       p_w * top_15_data$Potency_P) * 9
  
  new_mlss <- new_base_score + top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty
  scenario_results[[scenario_name]] <- new_mlss
}

# ── Ranking correlation analysis ──────────────────────────────────────────────
rank_data          <- apply(scenario_results[, -1], 2, rank)
rank_correlation   <- cor(rank_data, method = "spearman")
original_correlations <- rank_correlation[1, -1]   # correlations vs. original

# ── Top-5 stability across scenarios ─────────────────────────────────────────
top_5_variations <- list()

for (i in 2:ncol(scenario_results)) {
  scenario_name   <- colnames(scenario_results)[i]
  top5_scenario   <- which(rank(scenario_results[[i]],             ties.method = "first") <= 5)
  top5_original   <- which(rank(scenario_results$Original_MLSS,    ties.method = "first") <= 5)
  top_5_variations[[scenario_name]] <- length(intersect(top5_scenario, top5_original))
}

# ── Ablation study ────────────────────────────────────────────────────────────
ablation_results <- data.frame(Combination  = top_15_data$Combination,
                               Original_MLSS = top_15_data$MLSS_v4.0,
                               stringsAsFactors = FALSE)

ablation_results$Without_C <- ((0.30/0.55 * top_15_data$Balance_B +
                                  0.10/0.55 * top_15_data$Network_Coverage_V +
                                  0.15/0.55 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_B <- ((0.45/0.70 * top_15_data$Complementarity_C +
                                  0.10/0.70 * top_15_data$Network_Coverage_V +
                                  0.15/0.70 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_V <- ((0.45/0.90 * top_15_data$Complementarity_C +
                                  0.30/0.90 * top_15_data$Balance_B +
                                  0.15/0.90 * top_15_data$Potency_P) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Without_P <- ((0.45/0.85 * top_15_data$Complementarity_C +
                                  0.30/0.85 * top_15_data$Balance_B +
                                  0.10/0.85 * top_15_data$Network_Coverage_V) * 9) +
  top_15_data$Synergy_Bonus - top_15_data$Antagonism_Penalty

ablation_results$Base_Only  <- (0.45 * top_15_data$Complementarity_C +
                                  0.30 * top_15_data$Balance_B +
                                  0.10 * top_15_data$Network_Coverage_V +
                                  0.15 * top_15_data$Potency_P) * 9

ablation_impact <- data.frame(
  Component     = c("Complementarity (C)", "Balance (B)", "Coverage (V)",
                    "Potency (P)", "Synergy/Antagonism"),
  Avg_Rank_Shift = NA,
  Avg_Score_Drop = NA,
  Top5_Stability = NA,
  NPG_Color      = npg_colors[1:5],
  stringsAsFactors = FALSE
)

ablation_cols <- c("Without_C", "Without_B", "Without_V", "Without_P", "Base_Only")

for (i in seq_along(ablation_cols)) {
  ablation_col <- ablation_cols[i]
  rank_diff    <- abs(rank(ablation_results$Original_MLSS) - rank(ablation_results[[ablation_col]]))
  ablation_impact$Avg_Rank_Shift[i] <- mean(rank_diff, na.rm = TRUE)
  ablation_impact$Avg_Score_Drop[i] <- mean(ablation_results$Original_MLSS - ablation_results[[ablation_col]], na.rm = TRUE)
  top5_orig    <- which(rank(ablation_results$Original_MLSS,          ties.method = "first") <= 5)
  top5_abl     <- which(rank(ablation_results[[ablation_col]],         ties.method = "first") <= 5)
  ablation_impact$Top5_Stability[i] <- length(intersect(top5_orig, top5_abl)) / 5
}

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Sensitivity heatmap  →  Fig_S2a
# ═══════════════════════════════════════════════════════════════════════════════
heatmap_data <- scenario_results %>%
  select(-Original_MLSS) %>%
  arrange(match(Combination, top_15$Combination))

rownames(heatmap_data) <- heatmap_data$Combination
heatmap_data           <- heatmap_data[, -1]

heatmap_normalized <- apply(heatmap_data, 2, function(x)
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

heatmap_long <- reshape2::melt(as.matrix(heatmap_normalized))
colnames(heatmap_long) <- c("Combination", "Scenario", "Normalized_Score")

p_heatmap <- ggplot(heatmap_long, aes(x = Scenario, y = Combination, fill = Normalized_Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = npg_colors[6], high = npg_colors[1],
                      name = "Normalized\nMLSS Score") +
  labs(title = "Weight Scenario Sensitivity (Top 15)",
       x = "Weight Scenario", y = "Drug Combination") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.title   = element_text(size = 12, color = "black", face = "bold"),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    panel.grid   = element_blank()
  )

ggsave("Fig_S2a_mlss_sensitivity_top15.png", p_heatmap,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_S2a_mlss_sensitivity_top15.svg", p_heatmap,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Weight robustness bar chart  →  Fig_6d
# ═══════════════════════════════════════════════════════════════════════════════
top_5_variations_filtered <- top_5_variations[names(top_5_variations) != "Original_MLSS"]
alternative_correlations  <- original_correlations[names(original_correlations) != "MLSS (Standard)"]

robustness_data <- data.frame(
  Scenario       = names(top_5_variations_filtered),
  Top5_Stability = unlist(top_5_variations_filtered) / 5,
  stringsAsFactors = FALSE
)

robustness_data$Rank_Correlation <- c(1.0, alternative_correlations)

robustness_data$NPG_Color <- sapply(robustness_data$Scenario, function(s) {
  idx <- which(weight_scenarios$Scenario == s)
  if (length(idx) > 0) weight_scenarios$NPG_Color[idx] else npg_colors[1]
})

robustness_long <- robustness_data %>%
  pivot_longer(
    cols      = c(Top5_Stability, Rank_Correlation),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = factor(Metric,
                         levels = c("Top5_Stability", "Rank_Correlation"),
                         labels = c("Top 5 Stability", "Rank Correlation")))

p_robustness <- ggplot(robustness_long,
                       aes(x = reorder(Scenario, -Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge",
           width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Top 5 Stability" = npg_colors[1], "Rank Correlation" = npg_colors[2]),
    name   = "Robustness Metric"
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed",  color = npg_colors[6], linewidth = 0.5) +
  geom_hline(yintercept = 0.85, linetype = "dotted",  color = npg_colors[6], linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "MLSS Robustness to Weight Changes (Top 15)",
    x     = "Weight Scenario",
    y     = "Robustness Score (0–1)"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title         = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x        = element_text(size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    plot.title         = element_text(size = 16, face = "bold", hjust = 1),
    legend.position    = "bottom",
    legend.title       = element_text(size = 12, color = "black", face = "bold"),
    legend.text        = element_text(size = 12, color = "black"),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("Fig_6d_Weight_Robustness_Top15.png", plot = p_robustness,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_6d_Weight_Robustness_Top15.svg", plot = p_robustness,
       width = 7.7, height = 6, units = "cm", scale = 2.3, dpi = 600, bg = "white")

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 17: Fig_6d_Fig_S2_MLSS_validation.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 18/21  —  Fig_6e_drug_combination_network.r                           ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 6e: Drug combination network — top 15 combinations (circular layout)
# Lines: 207  |  File: Fig_6e_drug_combination_network.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Drug Combination Network - Figure 6e

rm(list = ls())
set.seed(42)

library(tidyverse)
library(igraph)
library(ggraph)
library(scales)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# ── Load MLSS data ─────────────────────────────────────────────────────────────
mlss_all <- read.csv("Fig_6c_mlss_all_combinations.csv", stringsAsFactors = FALSE)

# Normalise column names (trim whitespace, handle R mangling)
colnames(mlss_all) <- trimws(colnames(mlss_all))

# ── Select top N combinations ──────────────────────────────────────────────────
top_n <- 15

# Auto-detect cross-pathway column or derive from Pathway_A vs Pathway_B
if ("Cross_Pathway" %in% colnames(mlss_all)) {
  cross_col <- "Cross_Pathway"
} else if ("Cross.Pathway" %in% colnames(mlss_all)) {
  cross_col <- "Cross.Pathway"
} else {
  cross_col <- NULL
}

top_combos <- mlss_all %>%
  arrange(Rank) %>%
  slice(1:min(top_n, nrow(.))) %>%
  mutate(
    Rank = row_number(),
    Is_CrossPathway = if (!is.null(cross_col)) {
      ifelse(.data[[cross_col]] == "Yes", "Cross-Pathway", "Same-Pathway")
    } else {
      ifelse(Pathway_A != Pathway_B, "Cross-Pathway", "Same-Pathway")
    }
  )

# ── Build edge list ────────────────────────────────────────────────────────────
edge_list <- top_combos %>%
  transmute(
    from         = Drug_A,
    to           = Drug_B,
    weight       = MLSS_v4.0,
    relationship = Is_CrossPathway,
    rank         = Rank
  )

# ── All unique drugs ───────────────────────────────────────────────────────────
all_drugs <- sort(unique(c(top_combos$Drug_A, top_combos$Drug_B)))

# ── Assign primary pathway to each drug ───────────────────────────────────────
drug_pathways <- bind_rows(
  top_combos %>% select(Drug = Drug_A, Pathway = Pathway_A),
  top_combos %>% select(Drug = Drug_B, Pathway = Pathway_B)
) %>%
  group_by(Drug, Pathway) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Drug) %>%
  arrange(Drug, desc(n)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  group_by(Drug) %>%
  summarise(
    Primary_Pathway = if (n_distinct(Pathway) == 1) {
      first(Pathway)
    } else if (n_distinct(Pathway) > 1 && first(n) > nth(n, 2)) {
      first(Pathway)
    } else {
      "Multi-pathway"
    },
    .groups = "drop"
  )

# ── Node list ──────────────────────────────────────────────────────────────────
node_list <- tibble(name = all_drugs) %>%
  left_join(drug_pathways, by = c("name" = "Drug")) %>%
  mutate(
    Primary_Pathway  = ifelse(is.na(Primary_Pathway), "Unknown", Primary_Pathway),
    appearance_count = sapply(name, function(d) {
      sum(top_combos$Drug_A == d | top_combos$Drug_B == d)
    })
  )

# ── NPG colors for pathways ────────────────────────────────────────────────────
pathway_colors_npg <- c(
  "Calcium Signaling"   = npg_colors[4],
  "Incretin Signaling"  = npg_colors[1],
  "Wnt Signaling"       = npg_colors[2],
  "Metabolic Regulation"= npg_colors[3],
  "Multi-pathway"       = npg_colors[6],
  "Unknown"             = "grey70"
)

# ── Create igraph object ───────────────────────────────────────────────────────
g_drug <- graph_from_data_frame(
  d        = edge_list,
  vertices = node_list,
  directed = FALSE
)

V(g_drug)$appearance_count <- node_list$appearance_count
V(g_drug)$Primary_Pathway  <- node_list$Primary_Pathway
E(g_drug)$weight           <- edge_list$weight
E(g_drug)$relationship     <- edge_list$relationship

hub_threshold       <- 4
V(g_drug)$is_hub    <- ifelse(V(g_drug)$appearance_count >= hub_threshold, "Hub", "Non-hub")

# ── MLSS edge width range ──────────────────────────────────────────────────────
mlss_min <- min(E(g_drug)$weight, na.rm = TRUE)
mlss_max <- max(E(g_drug)$weight, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT — Drug Combination Network  →  Fig_6e
# ══════════════════════════════════════════════════════════════════════════════
p_Fig_6e <- ggraph(g_drug, layout = "circle") +
  geom_edge_link(
    aes(width = weight, linetype = relationship),
    color       = "grey70",
    alpha       = 0.6,
    show.legend = TRUE
  ) +
  scale_edge_linetype_manual(
    name   = "Combination Type",
    values = c("Cross-Pathway" = "solid", "Same-Pathway" = "dotted"),
    guide  = guide_legend(order = 3,
                          override.aes = list(color = "black", width = 0.8))
  ) +
  scale_edge_width_continuous(
    name   = "MLSS Score",
    range  = c(0.5, 3.5),
    breaks = c(mlss_min, (mlss_min + mlss_max) / 2, mlss_max),
    labels = c("Lower", "Intermediate", "Higher"),
    guide  = guide_legend(order = 4, override.aes = list(color = "black"))
  ) +
  geom_node_point(
    aes(size = appearance_count, color = Primary_Pathway, fill = Primary_Pathway),
    shape  = 21,
    stroke = 1.2,
    alpha  = 0.95
  ) +
  scale_color_manual(
    name  = "Functional Category",
    values = pathway_colors_npg,
    guide  = guide_legend(order = 1)
  ) +
  scale_fill_manual(
    name  = "Functional Category",
    values = pathway_colors_npg,
    guide  = "none"
  ) +
  scale_size_continuous(
    name   = "Network Appearances",
    range  = c(5, 14),
    breaks = c(1, hub_threshold, max(node_list$appearance_count)),
    labels = c("1",
               paste0(hub_threshold, " (hub)"),
               as.character(max(node_list$appearance_count))),
    guide  = guide_legend(order = 2)
  ) +
  geom_node_text(
    aes(label    = name,
        fontface = ifelse(appearance_count >= hub_threshold, "bold", "plain")),
    repel        = TRUE,
    size         = 4.5,
    max.overlaps = 30
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.45, "cm")
  ) +
  labs(title = sprintf("Drug Combination Network (Top %d)", nrow(top_combos)))

ggsave("Fig_6e_drug_combination_network.png", p_Fig_6e,
       scale = 2.3, width = 7.7, height = 6,
       units = "cm", dpi = 600, bg = "white")
ggsave("Fig_6e_drug_combination_network.svg", p_Fig_6e,
       scale = 2.3, width = 7.7, height = 6,
       units = "cm", dpi = 600, bg = "white")

# ── Export statistics ──────────────────────────────────────────────────────────
node_stats <- node_list %>%
  mutate(Hub_Status = ifelse(appearance_count >= hub_threshold, "Hub", "Non-hub")) %>%
  arrange(desc(appearance_count))

# Explicit namespace prevents tidyverse masking igraph's as_data_frame()
edge_stats <- igraph::as_data_frame(g_drug, what = "edges") %>%
  arrange(desc(weight)) %>%
  rename(
    Gene1              = from,
    Gene2              = to,
    `MLSS Score`       = weight,
    `Combination Type` = relationship
  )

write.csv(node_stats, "Fig_6e_drug_node_statistics.csv", row.names = FALSE)
write.csv(edge_stats, "Fig_6e_drug_edge_statistics.csv", row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 18: Fig_6e_drug_combination_network.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 19/21  —  Fig_6f_regulatory_mechanistic_sankey.r                      ═
# ════════════════════════════════════════════════════════════════════════════════
# Fig 6f: Six-layer regulatory-mechanistic Sankey diagram
# Lines: 182  |  File: Fig_6f_regulatory_mechanistic_sankey.r  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# Regulatory Approval & Mechanistic Pathway Sankey - Figure 6f

rm(list = ls())
set.seed(42)

library(tidyverse)
library(ggalluvial)
library(scales)
library(ggsci)

npg_colors <- pal_npg("nrc")(10)

# Load data
mlss_top15 <- read.csv("Fig_6c_mlss_top15.csv", stringsAsFactors = FALSE)
gene_centrality <- read.csv("Fig_6b_gene_centrality_summary.csv", stringsAsFactors = FALSE)

# Load regulatory databases
fda_approved <- read.csv("FDA_Approved.csv", stringsAsFactors = FALSE)
ema_approved <- read.csv("EMA_Approved.csv", stringsAsFactors = FALSE)
pmda_approved <- read.csv("PMDA_Approved.csv", stringsAsFactors = FALSE)
triple_approved <- read.csv("FDA-EMA-PMDA_Approved.csv", stringsAsFactors = FALSE)

# Extract hub drugs from top 15
top_15 <- mlss_top15 %>% arrange(Rank) %>% slice(1:15)

drug_freq_top15 <- bind_rows(
  top_15 %>% select(Drug = Drug_A),
  top_15 %>% select(Drug = Drug_B)
) %>%
  count(Drug, name = "Frequency") %>%
  arrange(desc(Frequency))

hub_threshold <- 4
hub_drugs <- drug_freq_top15 %>%
  filter(Frequency >= hub_threshold) %>%
  pull(Drug)

# Standardize drug names for regulatory lookup
fda_drugs <- tolower(fda_approved[,2])
ema_drugs <- tolower(ema_approved[,2])
pmda_drugs <- tolower(pmda_approved[,2])
triple_drugs <- tolower(triple_approved[,2])

# Assign regulatory approval tiers
drug_approval <- drug_freq_top15 %>%
  filter(Drug %in% hub_drugs) %>%
  mutate(
    Drug_lower = tolower(Drug),
    FDA = Drug_lower %in% fda_drugs,
    EMA = Drug_lower %in% ema_drugs,
    PMDA = Drug_lower %in% pmda_drugs,
    Approval_Tier = case_when(
      Drug_lower %in% triple_drugs ~ "FDA-EMA-PMDA\n(Tier 1)",
      FDA & EMA ~ "FDA-EMA\n(Tier 2)",
      FDA ~ "FDA-Only\n(Tier 3)",
      TRUE ~ "Not Approved"
    ),
    Confidence = case_when(
      Approval_Tier == "FDA-EMA-PMDA\n(Tier 1)" ~ 3.0,
      Approval_Tier == "FDA-EMA\n(Tier 2)" ~ 2.0,
      Approval_Tier == "FDA-Only\n(Tier 3)" ~ 1.0,
      TRUE ~ 0.3
    )
  ) %>%
  select(Drug, Frequency, Approval_Tier, Confidence, FDA, EMA, PMDA)

# Map drugs to pathways
drug_pathway_map <- bind_rows(
  top_15 %>% select(Drug = Drug_A, Pathway = Pathway_A),
  top_15 %>% select(Drug = Drug_B, Pathway = Pathway_B)
) %>%
  filter(Drug %in% hub_drugs) %>%
  group_by(Drug) %>%
  summarise(Primary_Pathway = first(Pathway), .groups = "drop")

# Define action types
drug_action_map <- tribble(
  ~Drug, ~Action_Type,
  "isradipine", "Antagonist",
  "(S)-nitrendipine", "Antagonist",
  "nitrendipine", "Antagonist",
  "albiglutide", "Agonist",
  "liraglutide", "Agonist"
)

# Map pathways to hub genes
pathway_gene_map <- tribble(
  ~Primary_Pathway, ~Hub_Gene, ~Connection_Strength,
  "Calcium Signaling", "Cacna1d", 1.0,
  "Calcium Signaling", "Cacna1c", 0.9,
  "Incretin Signaling", "Glp1r", 1.0,
  "Incretin Signaling", "Tcf7", 0.85,
  "Incretin Signaling", "Ptbp1", 0.75
)

# Map genes/pathways to clinical outcomes
pathway_outcome_map <- tribble(
  ~Primary_Pathway, ~Clinical_Outcome,
  "Calcium Signaling", "Enhanced insulin\nsecretion",
  "Calcium Signaling", "β-cell\npreservation",
  "Incretin Signaling", "Enhanced insulin\nsecretion",
  "Incretin Signaling", "β-cell\npreservation",
  "Incretin Signaling", "Metabolic\nhomeostasis"
)

# Build 6-layer alluvial dataset
alluvial_data <- drug_approval %>%
  left_join(drug_pathway_map, by = "Drug") %>%
  left_join(drug_action_map, by = "Drug") %>%
  left_join(pathway_gene_map, by = "Primary_Pathway") %>%
  left_join(pathway_outcome_map, by = "Primary_Pathway") %>%
  mutate(Flow = Frequency * Connection_Strength * Confidence) %>%
  select(`Approval Status` = Approval_Tier, `Hub Drugs` = Drug,
         `Action Type` = Action_Type, `Hub Genes` = Hub_Gene,
         `Pathways` = Primary_Pathway, `Clinical Outcomes` = Clinical_Outcome,
         Flow) %>%
  group_by(`Approval Status`, `Hub Drugs`, `Action Type`, `Hub Genes`, `Pathways`, `Clinical Outcomes`) %>%
  summarise(Flow = sum(Flow), .groups = "drop") %>%
  filter(Flow > 0) %>%
  arrange(desc(Flow))

# NPG colors
approval_colors <- c(
  "FDA-EMA-PMDA\n(Tier 1)" = npg_colors[1],
  "FDA-EMA\n(Tier 2)" = npg_colors[3],
  "FDA-Only\n(Tier 3)" = npg_colors[5],
  "Not Approved" = "#A6A6A6"
)

pathway_colors <- c(
  "Calcium Signaling" = npg_colors[4],
  "Incretin Signaling" = npg_colors[1]
)

action_colors <- c(
  "Agonist" = npg_colors[2],
  "Antagonist" = npg_colors[4]
)

outcome_colors <- c(
  "Enhanced insulin\nsecretion" = npg_colors[2],
  "β-cell\npreservation" = npg_colors[6],
  "Metabolic\nhomeostasis" = npg_colors[3]
)

# Create alluvial diagram
p_sankey <- ggplot(
  alluvial_data,
  aes(axis1 = `Approval Status`, axis2 = `Hub Drugs`, axis3 = `Action Type`,
      axis4 = `Hub Genes`, axis5 = `Pathways`, axis6 = `Clinical Outcomes`, y = Flow)
) +
  geom_alluvium(aes(fill = `Pathways`), width = 3/8, alpha = 0.7, curve_type = "sigmoid") +
  geom_stratum(width = 3/8, fill = "white", color = "white", size = 0.6, alpha = 0.8) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(name = "Signaling Pathway", values = pathway_colors,
                    guide = guide_legend(override.aes = list(alpha = 0.9, size = 5), order = 1)) +
  scale_x_discrete(limits = c("Approval\nStatus", "Hub\nDrugs", "Action\nType",
                               "Hub\nGenes", "Pathways", "Clinical\nOutcomes"),
                   expand = c(0.08, 0.08)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "6-Layer Regulatory-Mechanistic Network for Type 2 Diabetes Drug Repositioning")

ggsave("Fig_6f_regulatory_mechanistic_sankey.png", p_sankey,
       width = 17, height = 4.47, units = "cm", scale = 2.3, dpi = 600, bg = "white")
ggsave("Fig_6f_regulatory_mechanistic_sankey.svg", p_sankey,
       width = 17, height = 4.47, units = "cm", scale = 2.3, dpi = 600, bg = "white")

# Export supplementary tables
write.csv(drug_approval, "Fig_6f_hub_drug_approval_profile.csv", row.names = FALSE)
write.csv(alluvial_data, "Fig_6f_six_layer_complete_mapping.csv", row.names = FALSE)
write.csv(pathway_gene_map, "Fig_6f_pathway_gene_connectivity.csv", row.names = FALSE)

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 19: Fig_6f_regulatory_mechanistic_sankey.r
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 20/21  —  0_Analysis Script Combination.R                             ═
# ════════════════════════════════════════════════════════════════════════════════
# Script: 0_Analysis Script Combination.R
# Lines: 144  |  File: 0_Analysis Script Combination.R  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Combine All R Scripts into a Single Annotated File
# Filenames matched to actual directory contents
# Output: Combined_Analysis_Scripts.R
# ═══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

out_file      <- "Combined_Analysis_Scripts.R"
author_note   <- "Oontawee et al. — Notch-Wnt Pathway Inhibition in β-Cells (npj SBA, 2026)"
section_width <- 80

banner <- function(text, char = "═", width = section_width) {
  bar <- paste(rep(char, width), collapse = "")
  pad <- paste(rep(" ", max(0, width - nchar(text) - 4)), collapse = "")
  c(paste0("# ", bar),
    paste0("# ", char, " ", text, pad, " ", char),
    paste0("# ", bar))
}

# ── Actual script filenames and descriptions ──────────────────────────────────
script_order <- c(
  
  # ── Figure 1 ────────────────────────────────────────────────────────────────
  "Fig_1c_Boxplot_proliferation_revised.r"     = "Fig 1c: Total DNA content — proliferation boxplot (ANOVA/KW + Bonferroni)",
  "Fig_1d_Boxplot_Notch_Wnt_genes.r"           = "Fig 1d: Notch/Wnt gene expression boxplots — Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Tcf7, Lef1, Tcf7l2",
  
  # ── Figure 2 ────────────────────────────────────────────────────────────────
  "Fig_2b_glucose_uptake_revised.r"            = "Fig 2b: 2-NBDG glucose uptake fluorescence boxplot",
  "Fig_2cdef_functional_genes_revised.r"       = "Fig 2c–2f: Glut2, Kcnj11, Cacna1c, Cacna1d fold change boxplots",
  
  # ── Figure 3 ────────────────────────────────────────────────────────────────
  "Fig_3a_C-peptide_curve.r"                   = "Fig 3a: Glucose-stimulated C-peptide secretion LOESS curve",
  "Fig_3b_C-peptide_5.5mm_revised.r"           = "Fig 3b: C-peptide secretion at 5.5 mM glucose — boxplot",
  "Fig_3ef_Boxplot_script.r"                   = "Fig 3c–3f: Glp1r, Ins1, Ptbp1, Itpr1, Ryr1–3 gene expression boxplots",
  "Fig_3g_PCA_script.r"                        = "Fig 3g: PCA of 11 secretory-pathway genes (biplot)",
  
  # ── Figure 4 ────────────────────────────────────────────────────────────────
  "Fig_4b_correlation_matrix.r"                = "Fig 4b: Pearson correlation matrix heatmap (Hes1-normalized)",
  "Fig_4c_normalization_methods.r"             = "Fig 4c: Multi-method normalization comparison — 4 strategies",
  "Fig_4d_Venn_diagram.r"                      = "Fig 4d: Four-way Venn diagram — 28-pair consensus network",
  "Fig_4e_Hes1-normalized_correlation_matrix.R"= "Fig 4e: Hes1-normalized annotated correlation heatmap (consensus pairs marked)",
  
  # ── Figure S1 ───────────────────────────────────────────────────────────────
  "Fig_S1_validation_analysis.r"               = "Fig S1: Validation analysis — observed vs null distribution, false positive rates, sign-flip check",
  
  # ── Figure 5 ────────────────────────────────────────────────────────────────
  "Fig_5abcd_network_analysis.r"               = "Fig 5a–5d: Gene interaction network — 15 nodes, STRING classification, novel pair ranking",
  
  # ── Figure 6 ────────────────────────────────────────────────────────────────
  "Fig_6b_consensus_hub_network.r"             = "Fig 6b: Consensus hub gene network (circular layout, igraph/ggraph)",
  "Fig_6c_MLSS_drug_combinations.r"            = "Fig 6c: MLSS v4.0 scoring — all drug combinations ranked by synergy score",
  "Fig_6d_Fig_S2_MLSS_validation.r"            = "Fig 6d + Fig S2: MLSS robustness — weight sensitivity, ablation analysis (8 scenarios)",
  "Fig_6e_drug_combination_network.r"          = "Fig 6e: Drug combination network — top 15 combinations (circular layout)",
  "Fig_6f_regulatory_mechanistic_sankey.r"     = "Fig 6f: Six-layer regulatory-mechanistic Sankey diagram"
)

# ── Discover R files (exclude output file and old combine script) ─────────────
exclude <- c(out_file, "0_Combine CSV.R",
             "Fig_4c_normalization_methods_revised.r",
             "Fig_S1_validation_analysis_double.r")

all_r_files <- list.files(pattern = "\\.[Rr]$", full.names = FALSE)
all_r_files <- all_r_files[!all_r_files %in% exclude]

if (length(all_r_files) == 0) stop("No R script files found.")

known_order   <- names(script_order)[names(script_order) %in% all_r_files]
unknown_files <- sort(all_r_files[!all_r_files %in% known_order])
ordered_files <- c(known_order, unknown_files)

cat(sprintf("Found %d R script(s) — %d in order, %d additional\n\n",
            length(all_r_files), length(known_order), length(unknown_files)))

# ── Build output ──────────────────────────────────────────────────────────────
lines_out <- c(
  paste0("# ", paste(rep("═", section_width), collapse = "")),
  "#",
  "#  COMBINED ANALYSIS SCRIPTS",
  paste0("#  ", author_note),
  paste0("#  Generated: ", Sys.time()),
  paste0("#  Total scripts: ", length(ordered_files)),
  "#",
  "#  SCRIPT ORDER:",
  unlist(lapply(seq_along(ordered_files), function(i) {
    desc <- if (!is.null(script_order[[ordered_files[i]]])) {
      paste0(" — ", script_order[[ordered_files[i]]])
    } else ""
    sprintf("#    %2d. %s%s", i, ordered_files[i], desc)
  })),
  "#",
  paste0("# ", paste(rep("═", section_width), collapse = "")),
  ""
)

for (i in seq_along(ordered_files)) {
  file_name <- ordered_files[i]
  desc <- if (!is.null(script_order[[file_name]])) {
    script_order[[file_name]]
  } else paste0("Script: ", file_name)
  
  script_lines <- tryCatch(
    readLines(file_name, warn = FALSE),
    error = function(e) { warning(paste("Cannot read:", file_name)); NULL }
  )
  
  if (is.null(script_lines)) {
    cat(sprintf("⚠  Skipped: %s\n", file_name)); next
  }
  
  lines_out <- c(lines_out,
                 "",
                 banner(sprintf("SCRIPT %02d/%02d  —  %s", i, length(ordered_files), file_name)),
                 paste0("# ", desc),
                 paste0("# Lines: ", length(script_lines),
                        "  |  File: ", file_name,
                        "  |  Added: ", Sys.Date()),
                 paste0("# ", paste(rep("─", section_width), collapse = "")),
                 "",
                 script_lines,
                 "",
                 paste0("# ", paste(rep("─", section_width), collapse = "")),
                 paste0("# END OF SCRIPT ", sprintf("%02d", i), ": ", file_name),
                 paste0("# ", paste(rep("─", section_width), collapse = "")),
                 ""
  )
  
  cat(sprintf("✓  [%02d/%02d] %-42s (%d lines)\n",
              i, length(ordered_files), file_name, length(script_lines)))
}

lines_out <- c(lines_out,
               "",
               paste0("# ", paste(rep("═", section_width), collapse = "")),
               "# END OF COMBINED SCRIPT",
               paste0("# Total source files: ", length(ordered_files)),
               paste0("# Generated: ", Sys.time()),
               paste0("# ", paste(rep("═", section_width), collapse = ""))
)

writeLines(lines_out, out_file, useBytes = FALSE)
total <- length(readLines(out_file, warn = FALSE))
cat(sprintf("\n✓ Saved: %s  (%d total lines, %d scripts)\n",
            out_file, total, length(ordered_files)))

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 20: 0_Analysis Script Combination.R
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# ═ SCRIPT 21/21  —  0_Directory File Checker.R                                  ═
# ════════════════════════════════════════════════════════════════════════════════
# Script: 0_Directory File Checker.R
# Lines: 294  |  File: 0_Directory File Checker.R  |  Added: 2026-04-23
# ────────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Directory File Checker
# Scans the working directory and reports all R scripts and CSV files,
# cross-referencing against expected project outputs (Fig. 1–6).
# Exports: Directory_Check_Report.csv
# ═══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

cat("\n")
cat("████████████████████████████████████████████████████████████████████\n")
cat(" DIRECTORY FILE CHECKER\n")
cat(" Notch-Wnt Pathway Study — Fig. 1–6 Expected Outputs\n")
cat(sprintf(" Working directory: %s\n", getwd()))
cat(sprintf(" Checked: %s\n", Sys.time()))
cat("████████████████████████████████████████████████████████████████████\n\n")

# ── Expected R scripts ────────────────────────────────────────────────────────
expected_r <- c(
  "Fig_1c.R"           = "Fig 1c: DNA content proliferation analysis",
  "Fig_1d.R"           = "Fig 1d: Notch/Wnt heatmap — 10 genes",
  "Fig_2b.R"           = "Fig 2b: 2-NBDG glucose uptake",
  "Fig_2c_2f.R"        = "Fig 2c–2f: Glut2, Kcnj11, Cacna1c, Cacna1d",
  "Fig_3a_3b.R"        = "Fig 3a–3b: C-peptide secretion curves",
  "Fig_3c_3d.R"        = "Fig 3c–3d: Glp1r and Ins1 expression",
  "Fig_3e.R"           = "Fig 3e: Ptbp1, Itpr1, Ryr1–3 expression",
  "Fig_3f.R"           = "Fig 3f: PCA — 11 secretory genes",
  "Fig_4b_4c.R"        = "Fig 4b–4c: Multi-method correlation analysis",
  "Fig_4d.R"           = "Fig 4d: 28-pair consensus Venn",
  "Fig_4e.R"           = "Fig 4e: Hes1 correlation matrix heatmap",
  "FigS1a.R"           = "Fig S1a: Correlation by validation category",
  "FigS1b.R"           = "Fig S1b: False positive rate permutation",
  "Fig_5a.R"           = "Fig 5a: Gene interaction network",
  "Fig_5b_5d.R"        = "Fig 5b–5d: Interaction classification + top pairs",
  "Fig_6b.R"           = "Fig 6b: Consensus hub gene network",
  "Fig_6c_MLSS_v4.0.R" = "Fig 6c: MLSS v4.0 all combinations",
  "Fig_6d.R"           = "Fig 6d: MLSS robustness + ablation",
  "Fig_6e.R"           = "Fig 6e: Drug combination network"
)

# ── Expected CSV outputs ──────────────────────────────────────────────────────
expected_csv <- c(
  "Fig_1c_stats.csv"                               = "Fig 1c: DNA content statistics",
  "Fig_1d_Hes1_stats.csv"                          = "Fig 1d: Hes1 stats",
  "Fig_1d_Hey1_stats.csv"                          = "Fig 1d: Hey1 stats",
  "Fig_1d_Wnt2_stats.csv"                          = "Fig 1d: Wnt2 stats",
  "Fig_1d_Wnt2b_stats.csv"                         = "Fig 1d: Wnt2b stats",
  "Fig_1d_Wnt5a_stats.csv"                         = "Fig 1d: Wnt5a stats",
  "Fig_1d_Wnt5b_stats.csv"                         = "Fig 1d: Wnt5b stats",
  "Fig_1d_Wnt9a_stats.csv"                         = "Fig 1d: Wnt9a stats",
  "Fig_1d_Tcf7_stats.csv"                          = "Fig 1d: Tcf7 stats",
  "Fig_1d_Lef1_stats.csv"                          = "Fig 1d: Lef1 stats",
  "Fig_1d_Tcf7l2_stats.csv"                        = "Fig 1d: Tcf7l2 stats",
  "Fig_2b_stats.csv"                               = "Fig 2b: Glucose uptake stats",
  "Fig_2c_Glut2_stats.csv"                         = "Fig 2c: Glut2 stats",
  "Fig_2d_Kcnj11_stats.csv"                        = "Fig 2d: Kcnj11 stats",
  "Fig_2e_Cacna1c_stats.csv"                       = "Fig 2e: Cacna1c stats",
  "Fig_2f_Cacna1d_stats.csv"                       = "Fig 2f: Cacna1d stats",
  "Fig_3b_stats.csv"                               = "Fig 3b: C-peptide secretion stats",
  "Fig_3c_stats.csv"                               = "Fig 3c: Glp1r stats",
  "Fig_3d_stats.csv"                               = "Fig 3d: Ins1 stats",
  "Fig_3e_Ptbp1_stats.csv"                         = "Fig 3e: Ptbp1 stats",
  "Fig_3e_Itpr1_stats.csv"                         = "Fig 3e: Itpr1 stats",
  "Fig_3e_Ryr1_stats.csv"                          = "Fig 3e: Ryr1 stats",
  "Fig_3e_Ryr2_stats.csv"                          = "Fig 3e: Ryr2 stats",
  "Fig_3e_Ryr3_stats.csv"                          = "Fig 3e: Ryr3 stats",
  "Fig_3f_PCA_variance.csv"                        = "Fig 3f: PCA variance explained",
  "Fig_3f_PCA_coordinates.csv"                     = "Fig 3f: PCA sample/gene coordinates",
  "Fold_Change_deltadeltaCt.csv"                   = "ΔΔCt fold change — all 21 genes",
  "Fig_4c_Raw_expression_data.csv"                 = "Fig 4c: Raw 2^(-dCt) data",
  "Fig_4c_Hes1_normalized_data.csv"                = "Fig 4c: Hes1-normalized data",
  "Fig_4c_Kcnj11_normalized_data.csv"              = "Fig 4c: Kcnj11-normalized data",
  "Fig_4c_Ins1_normalized_data.csv"                = "Fig 4c: Ins1-normalized data",
  "Fig_4c_normalization_comparison_data.csv"       = "Fig 4c: Normalization comparison counts",
  "Fig_4c_permutation_null_stats.csv"              = "Fig 4c: Permutation null distribution",
  "Fig_4d_consensus_pairs.csv"                     = "Fig 4d: 28 consensus gene pairs",
  "Fig_4e_Hes1_significant_correlations_only.csv"  = "Fig 4e: Hes1 correlation matrix",
  "FigS1a_Correlation_Statistics.csv"              = "Fig S1a: Correlation statistics by category",
  "FigS1a_Correlation_Raw_Data.csv"                = "Fig S1a: All gene-pair correlations",
  "FigS1b_False_Positive_Statistics.csv"           = "Fig S1b: FP rate summary",
  "FigS1b_False_Positive_Permutation_Data.csv"     = "Fig S1b: Per-permutation FP rates",
  "Fig_5a_node_centrality.csv"                     = "Fig 5a: Node centrality metrics",
  "Fig_5a_edge_list.csv"                           = "Fig 5a: Full edge list",
  "Fig_5b_network_classification.csv"              = "Fig 5b: Gene pair classification",
  "Fig_5c_correlation_by_category.csv"             = "Fig 5c: |r| by interaction category",
  "Fig_5d_top15_novel_pairs.csv"                   = "Fig 5d: Top 15 novel gene pairs",
  "fig6b_gene_centrality_summary.csv"              = "Fig 6b: Hub network node centrality",
  "fig6b_edge_correlation_statistics.csv"          = "Fig 6b: Hub network edge statistics",
  "Fig_6c_mlss_all_combinations.csv"               = "Fig 6c: MLSS all combinations",
  "Fig_6d_ranking_correlations.csv"                = "Fig 6d: Rank correlation matrix",
  "Fig_6d_ablation_impact.csv"                     = "Fig 6d: Ablation impact summary",
  "Fig_S2a_mlss_sensitivity_data.csv"              = "Fig S2a: Sensitivity heatmap data",
  "Fig_S2b_component_importance_data.csv"          = "Fig S2b: Component importance data",
  "Fig_6e_drug_node_statistics.csv"                = "Fig 6e: Drug network node stats",
  "Fig_6e_drug_edge_statistics.csv"                = "Fig 6e: Drug network edge stats"
)

# ── Expected PNG/SVG outputs ──────────────────────────────────────────────────
expected_img <- c(
  "Fig_S2a_mlss_sensitivity_top15.png"      = "Fig S2a: Sensitivity heatmap (PNG)",
  "Fig_S2a_mlss_sensitivity_top15.svg"      = "Fig S2a: Sensitivity heatmap (SVG)",
  "Fig_S2b_component_importance_top15.png"  = "Fig S2b: Component importance (PNG)",
  "Fig_S2b_component_importance_top15.svg"  = "Fig S2b: Component importance (SVG)",
  "Fig_6d_Weight_Robustness_Top15.png"      = "Fig 6d: Weight robustness (PNG)",
  "Fig_6d_Weight_Robustness_Top15.svg"      = "Fig 6d: Weight robustness (SVG)",
  "fig6b_consensus_hub_network.png"         = "Fig 6b: Consensus hub network (PNG)",
  "fig6b_consensus_hub_network.svg"         = "Fig 6b: Consensus hub network (SVG)",
  "Fig_6e_drug_combination_network.png"     = "Fig 6e: Drug combination network (PNG)",
  "Fig_6e_drug_combination_network.svg"     = "Fig 6e: Drug combination network (SVG)"
)

# ── Scan directory ────────────────────────────────────────────────────────────
found_r   <- list.files(pattern = "\\.[Rr]$",        full.names = FALSE)
found_csv <- list.files(pattern = "\\.csv$",          full.names = FALSE, ignore.case = TRUE)
found_img <- list.files(pattern = "\\.(png|svg)$",    full.names = FALSE, ignore.case = TRUE)
found_all <- list.files(full.names = FALSE)

# ── Report helper (console) ───────────────────────────────────────────────────
check_files <- function(expected_named, found_vec, type_label) {
  expected_names <- names(expected_named)
  present <- expected_names[expected_names %in% found_vec]
  missing <- expected_names[!expected_names %in% found_vec]
  extra   <- found_vec[!found_vec %in% expected_names]
  
  cat(sprintf("── %s (%d expected, %d found, %d missing, %d unexpected) ──\n\n",
              type_label, length(expected_names),
              length(present), length(missing), length(extra)))
  
  if (length(present) > 0) {
    cat(sprintf("  ✓  PRESENT (%d):\n", length(present)))
    for (f in present) {
      info <- file.info(f)
      size <- if (!is.na(info$size)) sprintf("%.1f KB", info$size / 1024) else "?"
      cat(sprintf("       %-48s %8s  %s\n", f, size,
                  format(info$mtime, "%Y-%m-%d %H:%M")))
    }
    cat("\n")
  }
  
  if (length(missing) > 0) {
    cat(sprintf("  ✗  MISSING (%d):\n", length(missing)))
    for (f in missing)
      cat(sprintf("       %-48s  ← %s\n", f, expected_named[[f]]))
    cat("\n")
  }
  
  if (length(extra) > 0) {
    cat(sprintf("  ⚠  UNEXPECTED (%d):\n", length(extra)))
    for (f in extra) {
      info <- file.info(f)
      size <- if (!is.na(info$size)) sprintf("%.1f KB", info$size / 1024) else "?"
      cat(sprintf("       %-48s %8s\n", f, size))
    }
    cat("\n")
  }
  
  list(present = present, missing = missing, extra = extra)
}

# ── Run checks ────────────────────────────────────────────────────────────────
r_check   <- check_files(expected_r,   found_r,   "R SCRIPTS")
csv_check <- check_files(expected_csv, found_csv, "CSV OUTPUTS")
img_check <- check_files(expected_img, found_img, "PNG / SVG FIGURES")

# ── Console summary ───────────────────────────────────────────────────────────
cat("████████████████████████████████████████████████████████████████████\n")
cat(" SUMMARY\n")
cat("████████████████████████████████████████████████████████████████████\n\n")

total_expected <- length(expected_r) + length(expected_csv) + length(expected_img)
total_present  <- length(r_check$present) + length(csv_check$present) + length(img_check$present)
total_missing  <- length(r_check$missing) + length(csv_check$missing) + length(img_check$missing)
total_extra    <- length(r_check$extra)   + length(csv_check$extra)   + length(img_check$extra)

cat(sprintf("  R scripts   :  %2d / %2d present  (%d missing)\n",
            length(r_check$present),   length(expected_r),   length(r_check$missing)))
cat(sprintf("  CSV outputs :  %2d / %2d present  (%d missing)\n",
            length(csv_check$present), length(expected_csv), length(csv_check$missing)))
cat(sprintf("  PNG/SVG figs:  %2d / %2d present  (%d missing)\n",
            length(img_check$present), length(expected_img), length(img_check$missing)))
cat("  ─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  TOTAL       :  %2d / %2d present  (%d missing, %d unexpected)\n\n",
            total_present, total_expected, total_missing, total_extra))

if (total_missing == 0) {
  cat("  ✓  ALL EXPECTED FILES PRESENT — ready to build Table S3 and Combined Script.\n\n")
} else {
  cat(sprintf("  ✗  %d FILE(S) MISSING — run the corresponding analysis script(s).\n\n",
              total_missing))
}

# ── All files flat listing ────────────────────────────────────────────────────
cat("── ALL FILES IN WORKING DIRECTORY ──\n\n")
all_info <- file.info(found_all)
all_info <- all_info[order(all_info$mtime, decreasing = TRUE), ]

for (f in rownames(all_info)) {
  size  <- if (!is.na(all_info[f, "size"])) sprintf("%8.1f KB", all_info[f, "size"] / 1024) else "       ? KB"
  mtime <- format(all_info[f, "mtime"], "%Y-%m-%d %H:%M")
  ext   <- tolower(tools::file_ext(f))
  tag   <- switch(ext, "r" = "[R]  ", "csv" = "[CSV]", "png" = "[PNG]",
                  "svg" = "[SVG]", "xlsx" = "[XLS]", "txt" = "[TXT]", "[   ]")
  cat(sprintf("  %s  %-48s %s  %s\n", tag, f, size, mtime))
}

# ═══════════════════════════════════════════════════════════════════════════════
# CSV EXPORT — Directory_Check_Report.csv
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── EXPORTING CHECK REPORT ──\n\n")

# ── Helper: build rows for one file type ─────────────────────────────────────
build_rows <- function(expected_named, found_vec, file_type) {
  all_expected <- names(expected_named)
  extra        <- found_vec[!found_vec %in% all_expected]
  
  rows <- lapply(all_expected, function(f) {
    present  <- f %in% found_vec
    info     <- if (present) file.info(f) else NULL
    size_kb  <- if (present && !is.na(info$size)) round(info$size / 1024, 2) else NA
    mod_time <- if (present && !is.na(info$mtime)) format(info$mtime, "%Y-%m-%d %H:%M:%S") else NA
    
    data.frame(
      File            = f,
      File_Type       = file_type,
      Status          = ifelse(present, "Present", "Missing"),
      In_Manifest     = "Yes",
      Description     = expected_named[[f]],
      Size_KB         = size_kb,
      Last_Modified   = mod_time,
      stringsAsFactors = FALSE
    )
  })
  
  extra_rows <- lapply(extra, function(f) {
    info     <- file.info(f)
    size_kb  <- if (!is.na(info$size)) round(info$size / 1024, 2) else NA
    mod_time <- if (!is.na(info$mtime)) format(info$mtime, "%Y-%m-%d %H:%M:%S") else NA
    
    data.frame(
      File            = f,
      File_Type       = file_type,
      Status          = "Unexpected",
      In_Manifest     = "No",
      Description     = "",
      Size_KB         = size_kb,
      Last_Modified   = mod_time,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, c(rows, extra_rows))
}

# ── Build full report data frame ──────────────────────────────────────────────
report_df <- rbind(
  build_rows(expected_r,   found_r,   "R Script"),
  build_rows(expected_csv, found_csv, "CSV Output"),
  build_rows(expected_img, found_img, "PNG/SVG Figure")
)

# ── Add summary row ───────────────────────────────────────────────────────────
summary_row <- data.frame(
  File          = "── SUMMARY ──",
  File_Type     = "",
  Status        = sprintf("%d Present | %d Missing | %d Unexpected",
                          total_present, total_missing, total_extra),
  In_Manifest   = "",
  Description   = sprintf("Checked: %s  |  Working dir: %s", Sys.time(), getwd()),
  Size_KB       = NA,
  Last_Modified = "",
  stringsAsFactors = FALSE
)

report_df <- rbind(summary_row, report_df)

# ── Sort: Missing first, then Present, then Unexpected ───────────────────────
report_df[-1, ] <- report_df[-1, ][order(
  factor(report_df$Status[-1],
         levels = c("Missing", "Present", "Unexpected")),
  report_df$File_Type[-1],
  report_df$File[-1]
), ]

# ── Write CSV ─────────────────────────────────────────────────────────────────
out_csv <- "Directory_Check_Report.csv"
write.csv(report_df, out_csv, row.names = FALSE, na = "")

cat(sprintf("  ✓ Saved: %s  (%d rows)\n", out_csv, nrow(report_df) - 1))

cat("\n")
cat("████████████████████████████████████████████████████████████████████\n")
cat(sprintf(" Check complete: %s\n", Sys.time()))
cat("████████████████████████████████████████████████████████████████████\n\n")

# ────────────────────────────────────────────────────────────────────────────────
# END OF SCRIPT 21: 0_Directory File Checker.R
# ────────────────────────────────────────────────────────────────────────────────


# ════════════════════════════════════════════════════════════════════════════════
# END OF COMBINED SCRIPT
# Total source files: 21
# Generated: 2026-04-23 18:58:41.720304
# ════════════════════════════════════════════════════════════════════════════════

