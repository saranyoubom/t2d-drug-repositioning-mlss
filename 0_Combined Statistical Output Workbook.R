# ═══════════════════════════════════════════════════════════════════════════════
# Table S3 — Combined Statistical Output Workbook (Fig. 1–6)
# Filenames matched to actual directory contents
# Font: Arial 10 pt throughout
# ═══════════════════════════════════════════════════════════════════════════════

rm(list = ls())
library(openxlsx)

# ── Known descriptions: actual filenames from directory ───────────────────────
known_descriptions <- list(
  
  # ── Figure 1 ────────────────────────────────────────────────────────────────
  "Fig_1c_stats.csv"                               = "Fig 1c: Total DNA content (proliferation) — Shapiro-Wilk normality, one-way ANOVA/KW, Bonferroni pairwise comparisons (CTRL vs DAPT vs DKK-1)",
  "Fig_1d_stats.csv"                               = "Fig 1d: Notch/Wnt gene expression — combined statistical results for all 10 genes (Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Tcf7, Lef1, Tcf7l2)",
  "Fig_1d_Statistical_Results.csv"                 = "Fig 1d: Per-gene statistical results table — test statistic, p-value, adjusted p-value, significance label",
  "Fig_1d_Test_Type_Per_Gene.csv"                  = "Fig 1d: Test type selected per gene (ANOVA or Kruskal-Wallis) based on Shapiro-Wilk normality results",
  
  # ── Figure 2 ────────────────────────────────────────────────────────────────
  "Fig_2b_stats.csv"                               = "Fig 2b: 2-NBDG fluorescence intensity (glucose uptake) — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_2c_stats.csv"                               = "Fig 2c: Glut2 mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_2d_stats.csv"                               = "Fig 2d: Kcnj11 mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_2e_stats.csv"                               = "Fig 2e: Cacna1c mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_2f_stats.csv"                               = "Fig 2f: Cacna1d mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  
  # ── Figure 3 ────────────────────────────────────────────────────────────────
  "Fig_3a_LOESS_smooth_CI.csv"                     = "Fig 3a: LOESS-smoothed C-peptide secretion curves with 95% confidence intervals across glucose concentrations (CTRL, DAPT, DKK-1)",
  "Fig_3b_stats.csv"                               = "Fig 3b: C-peptide secretion at 5.5 mM glucose (ng/g DNA) — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_3c_stats.csv"                               = "Fig 3c: Glp1r mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_3d_stats.csv"                               = "Fig 3d: Ins1 mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_3e_stats.csv"                               = "Fig 3e: Calcium/RNA-binding genes — combined statistics for Ptbp1, Itpr1, Ryr1, Ryr2, Ryr3",
  "Fig_3f_Ptbp1_stats.csv"                         = "Fig 3f: Ptbp1 mRNA fold change — Shapiro-Wilk, one-way ANOVA/KW, Bonferroni pairwise comparisons",
  "Fig_3g_PCA_variance_explained.csv"              = "Fig 3f/3g: PCA variance explained per dimension — 11 secretory-pathway genes across CTRL, DAPT, DKK-1",
  "Fig_3g_PCA_scores.csv"                          = "Fig 3f/3g: PCA sample scores (coordinates) — treatment group separation across 11 selected genes",
  "Fig_3g_PCA_loadings.csv"                        = "Fig 3f/3g: PCA gene loadings (eigenvectors) — contribution of each gene to principal components",
  
  # ── Figure 4 ────────────────────────────────────────────────────────────────
  "Fig_4b_full_correlation_matrix.csv"             = "Fig 4b: Full Pearson correlation matrix (21 × 21 genes) — Hes1-normalized expression data, all pairwise r values",
  "Fig_4b_pvalue_matrix.csv"                       = "Fig 4b: Full p-value matrix (21 × 21 genes) — corresponding p-values for all Pearson correlations",
  "Fig_4b_significant_correlations.csv"            = "Fig 4b: Significant gene pairs only (p < 0.05) — r value, p-value, adjusted p-value",
  "Fig_4c_Raw_expression_data.csv"                 = "Fig 4c: Raw 2^(−dCt) expression data (21 genes, n = 12) used as input for all four normalization methods",
  "Fig_4c_Hes1_normalized_data.csv"                = "Fig 4c: Expression data normalized to Hes1 — 21 genes, n = 12",
  "Fig_4c_Kcnj11_normalized_data.csv"              = "Fig 4c: Expression data normalized to Kcnj11 — 21 genes, n = 12",
  "Fig_4c_Ins1_normalized_data.csv"                = "Fig 4c: Expression data normalized to Ins1 — 21 genes, n = 12",
  "Fig_4c_normalization_comparison_data.csv"       = "Fig 4c: Significant / Strong (r ≥ 0.7) / Very Strong (r ≥ 0.8) gene pair counts across four normalization methods",
  "Fig_4c_permutation_null_stats.csv"              = "Fig 4c: Empirical null distribution — mean ± SD significant pairs from 1,000 column-wise permutations (alpha = 0.05)",
  "Fig_4d_consensus_pairs.csv"                     = "Fig 4d: 28 consensus gene pairs identified across all four normalization methods (four-way Venn intersection)",
  "Fig_4d_Hes1_normalized_pairs.csv"               = "Fig 4d: Significant pairs from Hes1-normalized method only",
  "Fig_4d_Kcnj11_normalized_pairs.csv"             = "Fig 4d: Significant pairs from Kcnj11-normalized method only",
  "Fig_4d_Ins1_normalized_pairs.csv"               = "Fig 4d: Significant pairs from Ins1-normalized method only",
  "Fig_4d_Raw_expression_pairs.csv"                = "Fig 4d: Significant pairs from raw expression method only",
  "Fig_4e_Hes1_significant_correlations_only.csv"  = "Fig 4e: Pearson correlation matrix (Hes1-normalized) — significant pairs only (p < 0.05); 28 consensus pairs annotated",
  "Fig_4e_consensus_annotations.csv"               = "Fig 4e: Consensus annotation key — gene pair labels, r values, consensus status, STRING classification",
  
  # ── Figure S1 ───────────────────────────────────────────────────────────────
  "Fig_S1a_Observed_vs_Null.csv"                   = "Fig S1a: Observed significant pair counts vs permutation null distribution — all four normalization methods",
  "Fig_S1b_CI_width_pooled_vs_single.csv"          = "Fig S1b: False positive rate comparison — 4-Method Consensus (0.5%) vs Hes1-Only (4.2%); 1,000 permutation iterations",
  "Fig_S1c_Hey1_signflip.csv"                      = "Fig S1c: Sign-flip validation for Hey1 — correlation direction consistency across normalization methods",
  
  # ── Figure 5 ────────────────────────────────────────────────────────────────
  "Fig_5_All_Classified_Pairs.csv"                 = "Fig 5: Complete classified gene pair list — all 181 high-confidence pairs (r ≥ 0.7) with STRING annotation, category (Known Validated / Novel Discovery / STRING Only / Not Validated), and r value",
  "Fig_5_Summary_Statistics.csv"                   = "Fig 5: Summary statistics — pair counts and mean |r| by interaction category",
  "Fig_5b_Coverage_Summary.csv"                    = "Fig 5b: Network coverage summary — proportion of pairs in each validation category",
  "Fig_5d_Top_Novel_Discoveries.csv"               = "Fig 5d: Top novel gene pair discoveries — ranked by absolute Pearson r, absent from STRING (r ≥ 0.93)",
  
  # ── Figure 6b ───────────────────────────────────────────────────────────────
  "Fig_6b_gene_centrality_summary.csv"             = "Fig 6b: Consensus hub gene network — node centrality summary (degree, betweenness, closeness, eigenvector) and hub classification",
  "Fig_6b_edge_correlation_statistics.csv"         = "Fig 6b: Consensus hub gene network — edge list with Pearson r, |r|, edge type (positive/negative), and edge width",
  "fig6b_gene_centrality_summary.csv"              = "Fig 6b: Consensus hub gene network — node centrality summary (degree, betweenness, closeness, eigenvector) and hub classification",
  "fig6b_edge_correlation_statistics.csv"          = "Fig 6b: Consensus hub gene network — edge list with Pearson r, |r|, edge type (positive/negative), and edge width",
  
  # ── Figure 6c ───────────────────────────────────────────────────────────────
  "Fig_6c_mlss_all_combinations.csv"               = "Fig 6c: MLSS v4.0 scores for all drug combinations — Complementarity (C), Balance (B), Coverage (V), Potency (P), Synergy Bonus, Antagonism Penalty, final MLSS, and Rank",
  "Fig_6c_mlss_top15.csv"                          = "Fig 6c: Top 15 drug combinations by MLSS v4.0 — final ranked list with all scoring components",
  
  # ── Figure 6d ───────────────────────────────────────────────────────────────
  "Fig_6d_ranking_correlations.csv"                = "Fig 6d: Spearman rank correlation matrix — pairwise ranking stability across 8 MLSS weight scenarios (top 15 combinations)",
  "Fig_6d_ablation_impact.csv"                     = "Fig 6d: Ablation study impact — average rank shift and score drop per component (C, B, V, P, Synergy/Antagonism); top-5 stability fraction",
  "Fig_S2a_mlss_sensitivity_data.csv"              = "Fig S2a: Normalized MLSS scores for top 15 combinations across 8 weight scenarios (sensitivity heatmap raw data)",
  "Fig_S2b_component_importance_data.csv"          = "Fig S2b: Component importance ablation results — Avg Rank Shift and Avg Score Drop for each MLSS component",
  
  # ── Figure 6e ───────────────────────────────────────────────────────────────
  "Fig_6e_drug_node_statistics.csv"                = "Fig 6e: Drug combination network — node statistics (drug name, primary pathway, appearance count in top 15, hub status)",
  "Fig_6e_drug_edge_statistics.csv"                = "Fig 6e: Drug combination network — edge list (Drug A, Drug B, MLSS score, combination type: Cross-Pathway/Same-Pathway, rank)",
  
  # ── Figure 6f ───────────────────────────────────────────────────────────────
  "Fig_6f_hub_drug_approval_profile.csv"           = "Fig 6f: Hub drug regulatory approval profile — FDA/EMA/PMDA approval tier, approval year, indication, and approval status per drug",
  "Fig_6f_pathway_gene_connectivity.csv"           = "Fig 6f: Pathway-gene connectivity — hub gene to drug target mapping across calcium signaling, incretin, and metabolic pathways",
  "Fig_6f_six_layer_complete_mapping.csv"          = "Fig 6f: Six-layer complete mechanistic mapping — disease → pathway → gene → drug target → drug → approval status",
  
  # ── Reference / database files ──────────────────────────────────────────────
  "drug.target.interaction.csv"                    = "Reference: Drug-target interaction database — all drug-gene/protein interaction pairs used for network pharmacology scoring",
  "string_interactions.csv"                        = "Reference: STRING protein-protein interaction database export — interaction scores and evidence channels for 21 hub genes",
  "FDA_Approved.csv"                               = "Reference: FDA-approved drug list — drug names, targets, indications, and approval dates",
  "EMA_Approved.csv"                               = "Reference: EMA-approved drug list — drug names, targets, indications, and approval dates",
  "PMDA_Approved.csv"                              = "Reference: PMDA-approved drug list — drug names, targets, indications, and approval dates",
  "FDA-EMA-PMDA_Approved.csv"                      = "Reference: Triple-agency (FDA + EMA + PMDA) approved drugs — Tier 1 regulatory maturity drugs used in MLSS prioritization"
)

# ── Sheet name order (Missing/Reference sheets last) ──────────────────────────
priority_order <- c(
  "Fig_1c_stats.csv",
  "Fig_1d_stats.csv",
  "Fig_1d_Statistical_Results.csv",
  "Fig_1d_Test_Type_Per_Gene.csv",
  "Fig_2b_stats.csv",
  "Fig_2c_stats.csv",
  "Fig_2d_stats.csv",
  "Fig_2e_stats.csv",
  "Fig_2f_stats.csv",
  "Fig_3a_LOESS_smooth_CI.csv",
  "Fig_3b_stats.csv",
  "Fig_3c_stats.csv",
  "Fig_3d_stats.csv",
  "Fig_3e_stats.csv",
  "Fig_3f_Ptbp1_stats.csv",
  "Fig_3g_PCA_variance_explained.csv",
  "Fig_3g_PCA_scores.csv",
  "Fig_3g_PCA_loadings.csv",
  "Fig_4b_full_correlation_matrix.csv",
  "Fig_4b_pvalue_matrix.csv",
  "Fig_4b_significant_correlations.csv",
  "Fig_4c_Raw_expression_data.csv",
  "Fig_4c_Hes1_normalized_data.csv",
  "Fig_4c_Kcnj11_normalized_data.csv",
  "Fig_4c_Ins1_normalized_data.csv",
  "Fig_4c_normalization_comparison_data.csv",
  "Fig_4c_permutation_null_stats.csv",
  "Fig_4d_consensus_pairs.csv",
  "Fig_4d_Hes1_normalized_pairs.csv",
  "Fig_4d_Kcnj11_normalized_pairs.csv",
  "Fig_4d_Ins1_normalized_pairs.csv",
  "Fig_4d_Raw_expression_pairs.csv",
  "Fig_4e_Hes1_significant_correlations_only.csv",
  "Fig_4e_consensus_annotations.csv",
  "Fig_S1a_Observed_vs_Null.csv",
  "Fig_S1b_CI_width_pooled_vs_single.csv",
  "Fig_S1c_Hey1_signflip.csv",
  "Fig_5_All_Classified_Pairs.csv",
  "Fig_5_Summary_Statistics.csv",
  "Fig_5b_Coverage_Summary.csv",
  "Fig_5d_Top_Novel_Discoveries.csv",
  "Fig_6b_gene_centrality_summary.csv",
  "Fig_6b_edge_correlation_statistics.csv",
  "fig6b_gene_centrality_summary.csv",
  "fig6b_edge_correlation_statistics.csv",
  "Fig_6c_mlss_all_combinations.csv",
  "Fig_6c_mlss_top15.csv",
  "Fig_6d_ranking_correlations.csv",
  "Fig_6d_ablation_impact.csv",
  "Fig_S2a_mlss_sensitivity_data.csv",
  "Fig_S2b_component_importance_data.csv",
  "Fig_6e_drug_node_statistics.csv",
  "Fig_6e_drug_edge_statistics.csv",
  "Fig_6f_hub_drug_approval_profile.csv",
  "Fig_6f_pathway_gene_connectivity.csv",
  "Fig_6f_six_layer_complete_mapping.csv",
  "drug.target.interaction.csv",
  "string_interactions.csv",
  "FDA_Approved.csv",
  "EMA_Approved.csv",
  "PMDA_Approved.csv",
  "FDA-EMA-PMDA_Approved.csv"
)

# ── Helper: sanitize sheet names ──────────────────────────────────────────────
make_sheet_name <- function(filename) {
  nm <- tools::file_path_sans_ext(basename(filename))
  nm <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm)
  substr(nm, 1, 31)
}

# ── Auto-discover CSV files ───────────────────────────────────────────────────
csv_found <- list.files(pattern = "\\.csv$", full.names = FALSE, ignore.case = TRUE)
if (length(csv_found) == 0) stop("No CSV files found in the working directory.")

# Sort: priority order first, then remaining alphabetically
in_priority  <- priority_order[priority_order %in% csv_found]
not_priority <- sort(csv_found[!csv_found %in% priority_order])
csv_ordered  <- c(in_priority, not_priority)

cat(sprintf("Found %d CSV file(s) — %d in priority order, %d additional\n\n",
            length(csv_found), length(in_priority), length(not_priority)))

# ── Styles (Arial 10 pt) ──────────────────────────────────────────────────────
HDR_STYLE <- createStyle(
  fontName = "Arial", fontSize = 10, fontColour = "#FFFFFF",
  fgFill = "#1F3864", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", wrapText = TRUE,
  border = "TopBottomLeftRight", borderColour = "#FFFFFF"
)
DATA_STYLE <- createStyle(
  fontName = "Arial", fontSize = 10, halign = "LEFT",
  border = "TopBottomLeftRight", borderColour = "#D4D1CA"
)
NUM_STYLE <- createStyle(
  fontName = "Arial", fontSize = 10, halign = "RIGHT",
  numFmt = "0.000000",
  border = "TopBottomLeftRight", borderColour = "#D4D1CA"
)
META_TITLE <- createStyle(
  fontName = "Arial", fontSize = 11,
  fontColour = "#1F3864", textDecoration = "bold"
)
META_DESC <- createStyle(
  fontName = "Arial", fontSize = 10,
  fontColour = "#7A7974", wrapText = TRUE
)
ALT_FILL <- createStyle(fgFill = "#F3F0EC")

# ── Build workbook ────────────────────────────────────────────────────────────
wb         <- createWorkbook()
addWorksheet(wb, "Index")

index_rows <- list()
loaded     <- 0
skipped    <- character(0)
used_names <- character(0)

for (file_name in csv_ordered) {
  
  df <- tryCatch(
    read.csv(file_name, stringsAsFactors = FALSE),
    error = function(e) { warning(paste("Cannot read:", file_name)); NULL }
  )
  if (is.null(df) || nrow(df) == 0) {
    skipped <- c(skipped, file_name); next
  }
  
  base_name  <- make_sheet_name(file_name)
  sheet_name <- base_name
  i <- 2
  while (sheet_name %in% used_names) {
    sheet_name <- substr(paste0(base_name, "_", i), 1, 31); i <- i + 1
  }
  used_names <- c(used_names, sheet_name)
  
  description <- if (!is.null(known_descriptions[[file_name]])) {
    known_descriptions[[file_name]]
  } else {
    paste0("CSV export: ", file_name)
  }
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, x = sheet_name,  startRow = 1, startCol = 2)
  addStyle(wb, sheet_name, META_TITLE, rows = 1, cols = 2)
  writeData(wb, sheet_name, x = description, startRow = 2, startCol = 2)
  addStyle(wb, sheet_name, META_DESC,  rows = 2, cols = 2)
  mergeCells(wb, sheet_name, cols = 2:(2 + ncol(df) - 1), rows = 2)
  setRowHeights(wb, sheet_name, rows = 2, heights = 30)
  writeData(wb, sheet_name,
            x = paste0("Source: ", file_name, "  |  Generated: ", Sys.Date()),
            startRow = 3, startCol = 2)
  addStyle(wb, sheet_name, META_DESC, rows = 3, cols = 2)
  mergeCells(wb, sheet_name, cols = 2:(2 + ncol(df) - 1), rows = 3)
  
  data_start <- 5
  writeData(wb, sheet_name, df,
            startRow = data_start, startCol = 2,
            headerStyle = HDR_STYLE)
  
  last_row <- data_start + nrow(df)
  n_cols   <- ncol(df)
  num_cols <- which(sapply(df, is.numeric))
  
  # ── FAST (new) — column-range styling: O(cols) addStyle calls ────────────────
  data_rows <- (data_start + 1):last_row
  
  # Apply DATA_STYLE to all text columns at once
  text_cols <- setdiff(seq_len(n_cols), num_cols)
  if (length(text_cols) > 0) {
    for (ci in text_cols) {
      addStyle(wb, sheet_name, DATA_STYLE,
               rows = data_rows, cols = ci + 1,
               gridExpand = FALSE, stack = TRUE)
    }
  }
  
  # Apply NUM_STYLE to all numeric columns at once
  if (length(num_cols) > 0) {
    for (ci in num_cols) {
      addStyle(wb, sheet_name, NUM_STYLE,
               rows = data_rows, cols = ci + 1,
               gridExpand = FALSE, stack = TRUE)
    }
  }
  
  # Apply alternating fill to even rows only — vectorised
  alt_rows <- data_rows[(data_rows - data_start) %% 2 == 0]
  if (length(alt_rows) > 0) {
    addStyle(wb, sheet_name, ALT_FILL,
             rows = alt_rows, cols = 2:(n_cols + 1),
             gridExpand = TRUE, stack = TRUE)
  }
  
  setColWidths(wb, sheet_name, cols = 1,              widths = 3)
  setColWidths(wb, sheet_name, cols = 2:(n_cols + 1), widths = "auto")
  freezePane(wb, sheet_name,
             firstActiveRow = data_start + 1, firstActiveCol = 2)
  
  index_rows[[length(index_rows) + 1]] <- data.frame(
    Sheet       = sheet_name,
    File        = file_name,
    Rows        = nrow(df),
    Columns     = n_cols,
    Description = description,
    stringsAsFactors = FALSE
  )
  
  loaded <- loaded + 1
  cat(sprintf("✓  %-38s  (%d rows)\n", file_name, nrow(df)))
}

# ── Index sheet ───────────────────────────────────────────────────────────────
index_df <- do.call(rbind, index_rows)

writeData(wb, "Index",
          x = "Table S3: Combined Statistical Output — Notch-Wnt Pathway Inhibition in β-Cells (Fig. 1–6)",
          startRow = 1, startCol = 2)
addStyle(wb, "Index",
         createStyle(fontName = "Arial", fontSize = 14,
                     fontColour = "#1F3864", textDecoration = "bold"),
         rows = 1, cols = 2)
mergeCells(wb, "Index", cols = 2:6, rows = 1)

writeData(wb, "Index",
          x = paste0("Generated: ", Sys.Date(),
                     "  |  Sheets: ", loaded,
                     if (length(skipped) > 0)
                       paste0("  |  Skipped: ", paste(skipped, collapse = ", "))
                     else ""),
          startRow = 2, startCol = 2)
addStyle(wb, "Index",
         createStyle(fontName = "Arial", fontSize = 10,
                     fontColour = "#7A7974", wrapText = TRUE),
         rows = 2, cols = 2)
mergeCells(wb, "Index", cols = 2:6, rows = 2)

writeData(wb, "Index", index_df,
          startRow = 4, startCol = 2, headerStyle = HDR_STYLE)

if (!is.null(index_df) && nrow(index_df) > 0) {
  for (r in 5:(4 + nrow(index_df))) {
    alt <- (r - 4) %% 2 == 0
    for (ci in 1:ncol(index_df)) {
      style <- if (ci %in% c(3, 4)) NUM_STYLE else DATA_STYLE
      addStyle(wb, "Index", style, rows = r, cols = ci + 1, stack = TRUE)
    }
    if (alt) addStyle(wb, "Index", ALT_FILL,
                      rows = r, cols = 2:(ncol(index_df) + 1), stack = TRUE)
  }
}

setColWidths(wb, "Index", cols = 1,   widths = 3)
setColWidths(wb, "Index", cols = 2:6, widths = c(30, 45, 8, 8, 70))
freezePane(wb, "Index", firstActiveRow = 5, firstActiveCol = 2)

out_file <- "Table_S3_Statistical_Output.xlsx"
saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("\n✓ Saved: %s  (%d sheets + Index)\n", out_file, loaded))

