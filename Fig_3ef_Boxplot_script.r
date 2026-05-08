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

