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

