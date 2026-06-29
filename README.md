# Analysis Scripts — Revision 1 & Revision 2

**Repository:** [saranyoubom/t2d-drug-repositioning-mlss](https://github.com/saranyoubom/t2d-drug-repositioning-mlss)  
**Subfolder:** `R_scripts_Rev1/`  
**Generated:** 2026-05-09  
**Compiled source:** `Combined_Analysis_Scripts_Fig_1-6.txt`

---

## Overview

This repository contains the complete set of analysis scripts (R and Python) used to generate all figures and supplementary analyses for:

> **Systems Biology Analysis for Type 2 Diabetes Drug Repositioning via Multi-Layered Network Pharmacology and MLSS Algorithm**  
> *npj Systems Biology and Applications*

Scripts were updated across two revision rounds. Key revision-specific changes are noted per script below.

---

## Revision 2 Additions (June 2026)

One new R script was added to address Reviewer 2 Minor Comment 6:

| File | Figure | Description |
|------|--------|-------------|
| `Fig_S4_linregpcr_plot.R` | Supplementary Fig. 4 | Reads per-gene amplification efficiency CSVs and produces a 4 × 6 panel figure of representative LinRegPCR amplification curves with per-condition E-value annotations |

---

## Revision 1 Changes Summary

| Change | Scripts Affected |
|--------|-----------------|
| `Ryr1` excluded from calcium gene panels (near-undetectable expression in CTRL) | `Fig_3ef_Boxplot_script.r`, `Fig_4d_Venn_diagram.r`, `Fig_4e_Hes1-normalized_correlation_matrix.R` |
| `Hmisc` package replaced with base R `cor.test()` for portability | `Fig_4d_Venn_diagram.r`, `Fig_4e_Hes1-normalized_correlation_matrix.R` |
| Bonferroni-corrected p-values applied to C-peptide comparisons | `Fig_3b_C-peptide_5.5mm_revised.r` |
| Gene interaction category labels updated per reviewer response (R2.2) | `Fig_5abcd_network_analysis.r` |
| MLSS formula updated to v4.0 (January 2026) with potency multiplier | `Fig_6c_MLSS_drug_combinations.r` |
| Fig. S1 split into two validated panel scripts | `Fig_S1_validation_analysis.r`, `Fig_S1_validation_analysis_double.r` |

---

## Folder Contents

### Utility / Helper Scripts

| File | Lines | Description |
|------|-------|-------------|
| `0_Analysis Script Combination.R` | 77 | Compiles all `.R` scripts in the working directory into a single output file |
| `0_Combined Statistical Output Workbook.R` | 360 | Builds `Table_S3_Statistical_Output.xlsx` — all CSV statistical outputs assembled into a formatted Excel workbook (Arial 10 pt, alternating rows, Index sheet) |
| `0_Directory File Checker.R` | 294 | Scans the working directory and cross-references all R scripts, CSV outputs, and PNG/SVG figures against the expected project manifest; exports `Directory_Check_Report.csv` |

---

### Figure Scripts

#### Figure 1 — Notch/Wnt Pathway Inhibition & Gene Expression

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_1c_Boxplot_proliferation_revised.r` | 147 | Fig. 1c | Total DNA content (cell proliferation) — Shapiro-Wilk normality test, one-way ANOVA or Kruskal-Wallis, Bonferroni pairwise comparisons (CTRL vs DAPT vs DKK-1); outputs `Fig_1c_stats.csv` |
| `Fig_1d_Boxplot_Notch_Wnt_genes.r` | 375 | Fig. 1d | Notch target and Wnt pathway gene expression boxplots — 10 genes (Hes1, Hey1, Wnt2, Wnt2b, Wnt5a, Wnt5b, Wnt9a, Tcf7, Lef1, Tcf7l2); 5-column × 2-row layout (16.2 × 12.8 cm); outputs `Fig_1d_stats.csv`, `Fig_1d_Statistical_Results.csv`, `Fig_1d_Test_Type_Per_Gene.csv` |

#### Figure 2 — Glucose Uptake & Functional Gene Expression

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_2b_glucose_uptake_revised.r` | 166 | Fig. 2b | 2-NBDG fluorescence intensity (glucose uptake) — Shapiro-Wilk, ANOVA/KW, Bonferroni comparisons; outputs `Fig_2b_stats.csv` |
| `Fig_2cdef_functional_genes_revised.r` | 305 | Fig. 2c–f | Glut2, Kcnj11, Cacna1c, Cacna1d mRNA fold-change boxplots — unified theme matching Fig. 1d; outputs `Fig_2c_stats.csv`, `Fig_2d_stats.csv`, `Fig_2e_stats.csv`, `Fig_2f_stats.csv` |

#### Figure 3 — C-Peptide Secretion, Insulin Pathway & PCA

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_3a_C-peptide_curve.r` | 83 | Fig. 3a | LOESS-smoothed C-peptide secretion curves with 95% CI across glucose concentrations (CTRL, DAPT, DKK-1); outputs `Fig_3a_LOESS_smooth_CI.csv` |
| `Fig_3b_C-peptide_5.5mm_revised.r` | 165 | Fig. 3b | C-peptide secretion at 5.5 mM glucose — **Revision 1**: Bonferroni-corrected p-values applied; outputs `Fig_3b_stats.csv` |
| `Fig_3ef_Boxplot_script.r` | 360 | Fig. 3e–f | Calcium signalling & RNA-binding gene expression — **Revision 1**: Ryr1 excluded (near-undetectable CTRL); Fig. 3e: Ryr2, Ryr3, Itpr1 (16.2 × 5.3 cm); Fig. 3f: Ptbp1 (8 × 5.3 cm); outputs `Fig_3e_stats.csv`, `Fig_3f_Ptbp1_stats.csv` |
| `Fig_3g_PCA_script.r` | 111 | Fig. 3g | PCA of 11 secretory-pathway genes across CTRL, DAPT, DKK-1; outputs `Fig_3g_PCA_variance_explained.csv`, `Fig_3g_PCA_scores.csv`, `Fig_3g_PCA_loadings.csv` |

#### Figure 4 — Multi-Method Correlation Analysis

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_4b_correlation_matrix.r` | 60 | Fig. 4b | Full 21 × 21 Pearson correlation matrix (Hes1-normalized); outputs `Fig_4b_full_correlation_matrix.csv`, `Fig_4b_pvalue_matrix.csv`, `Fig_4b_significant_correlations.csv` |
| `Fig_4c_normalization_methods.r` | 119 | Fig. 4c | Comparison of four normalization methods (Raw, Hes1, Kcnj11, Ins1) — significant/strong/very-strong pair counts; outputs `Fig_4c_normalization_comparison_data.csv`, `Fig_4c_permutation_null_stats.csv` |
| `Fig_4c_normalization_methods_revised.r` | 146 | Fig. 4c | **Revision 1 version** — updated normalization comparison with revised thresholds and permutation null; all four expression datasets exported |
| `Fig_4d_Venn_diagram.r` | 99 | Fig. 4d | Four-way Venn diagram of 28 consensus gene pairs — **Revision 1**: `Hmisc` replaced with base R `cor.test()`; Ryr1 excluded; outputs `Fig_4d_consensus_pairs.csv`, normalization-specific pair CSVs |
| `Fig_4e_Hes1-normalized_correlation_matrix.R` | 214 | Fig. 4e | Hes1-normalized significant-only correlation heatmap with 28 consensus pair annotations — **Revision 1**: `Hmisc` removed; Ryr1 excluded; outputs `Fig_4e_Hes1_significant_correlations_only.csv`, `Fig_4e_consensus_annotations.csv` |

#### Figure 5 — Gene Interaction Network

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_5abcd_network_analysis.r` | 351 | Fig. 5a–d | Full gene interaction network with STRING validation — **Revision 1**: category labels updated per reviewer R2.2; 181 high-confidence pairs (|r| ≥ 0.7); outputs `Fig_5_All_Classified_Pairs.csv`, `Fig_5_Summary_Statistics.csv`, `Fig_5b_Coverage_Summary.csv`, `Fig_5d_Top_Novel_Discoveries.csv` |

#### Figure 6 — Hub Network, MLSS Drug Repositioning & Regulatory Approval

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_6b_consensus_hub_network.r` | 165 | Fig. 6b | Consensus hub gene network — node centrality (degree, betweenness, closeness, eigenvector); outputs `fig6b_gene_centrality_summary.csv`, `fig6b_edge_correlation_statistics.csv`, `fig6b_consensus_hub_network.png/.svg` |
| `Fig_6c_MLSS_drug_combinations.r` | 341 | Fig. 6c | **MLSS v4.0** — Multi-Layer Synergy Score for all drug combinations; formula: `MLSS = (0.45·C + 0.30·B + 0.10·V + 0.15·P) × 9 ± potency-weighted r`; outputs `Fig_6c_mlss_all_combinations.csv`, `Fig_6c_mlss_top15.csv` |
| `Fig_6d_Fig_S2_MLSS_validation.r` | 223 | Fig. 6d, S2 | MLSS robustness and ablation study — 8 weight scenarios, Spearman rank correlation matrix, component importance; outputs `Fig_6d_ranking_correlations.csv`, `Fig_6d_ablation_impact.csv`, `Fig_S2a_mlss_sensitivity_data.csv`, `Fig_S2b_component_importance_data.csv` |
| `Fig_6e_drug_combination_network.r` | 208 | Fig. 6e | Drug combination network — node/edge statistics for top-15 combinations; outputs `Fig_6e_drug_node_statistics.csv`, `Fig_6e_drug_edge_statistics.csv`, `Fig_6e_drug_combination_network.png/.svg` |
| `Fig_6f_regulatory_mechanistic_sankey.r` | 183 | Fig. 6f | Six-layer mechanistic Sankey diagram: Disease → Pathway → Gene → Target → Drug → Regulatory approval; outputs `Fig_6f_hub_drug_approval_profile.csv`, `Fig_6f_pathway_gene_connectivity.csv`, `Fig_6f_six_layer_complete_mapping.csv` |

#### Supplementary Figure S1 — Consensus Validation

| File | Lines | Figure Panel | Description |
|------|-------|-------------|-------------|
| `Fig_S1_validation_analysis.r` | 492 | Fig. S1a–c | Multi-normalization consensus validation — observed vs permutation null, false positive rate comparison, Hey1 sign-flip check; outputs `Fig_S1a_Observed_vs_Null.csv`, `Fig_S1b_CI_width_pooled_vs_single.csv`, `Fig_S1c_Hey1_signflip.csv` |
| `Fig_S1_validation_analysis_double.r` | 411 | Fig. S1a–b | **Double-validation version** — Observed vs permutation null bar chart + Pearson r CI validation plot; exports `FigS1a_Validation_Summary.csv`, `FigS1a_Validation_BarChart.png`, `FigS1b_Pearson_r_CI_Plot.png/.svg`; console reports CI (0.211–0.909) and r_min (0.576) |

---

### Master Combined Script

| File | Lines | Description |
|------|-------|-------------|
| `Combined_Analysis_Scripts_Fig_1-6.R` | 4,893 | Full concatenation of all 21 figure scripts in execution order — used for reproducible single-session runs |

---

## Required R Packages

Install all dependencies with:

```r
install.packages(c(
  "ggplot2", "ggpubr", "ggrepel", "ggcorrplot", "ggraph", "ggfortify",
  "ggalluvial", "ggsci",
  "dplyr", "tidyr", "tidyverse", "data.table", "reshape2",
  "readxl", "readr",
  "patchwork", "gridExtra", "grid",
  "igraph", "VennDiagram",
  "factoextra",
  "rstatix", "outliers", "scales", "rlang",
  "openxlsx"
))
```

> **Note:** `Hmisc` was used in pre-revision scripts but has been replaced with base R `cor.test()` in all Revision 1 scripts.

---

## Execution Order

For a full reproducible run, execute scripts in the following order:

```
1.  Fig_1c_Boxplot_proliferation_revised.r
2.  Fig_1d_Boxplot_Notch_Wnt_genes.r
3.  Fig_2b_glucose_uptake_revised.r
4.  Fig_2cdef_functional_genes_revised.r
5.  Fig_3a_C-peptide_curve.r
6.  Fig_3b_C-peptide_5.5mm_revised.r
7.  Fig_3ef_Boxplot_script.r
8.  Fig_3g_PCA_script.r
9.  Fig_4b_correlation_matrix.r
10. Fig_4c_normalization_methods_revised.r
11. Fig_4d_Venn_diagram.r
12. Fig_4e_Hes1-normalized_correlation_matrix.R
13. Fig_5abcd_network_analysis.r
14. Fig_6b_consensus_hub_network.r
15. Fig_6c_MLSS_drug_combinations.r
16. Fig_6d_Fig_S2_MLSS_validation.r
17. Fig_6e_drug_combination_network.r
18. Fig_6f_regulatory_mechanistic_sankey.r
19. Fig_S1_validation_analysis.r
20. Fig_S1_validation_analysis_double.r
21. 0_Combined Statistical Output Workbook.R   ← run last (reads all CSVs)
```

Or run the master script directly:
```r
source("Combined_Analysis_Scripts_Fig_1-6.R")
```

---

## Input Data Files Required

The following reference/input CSV files must be present in the working directory before running:

| File | Used by |
|------|---------|
| `drug.target.interaction.csv` | Fig_6c, Fig_6e, Fig_6f |
| `string_interactions.csv` | Fig_5abcd, Fig_4e |
| `FDA_Approved.csv` | Fig_6f |
| `EMA_Approved.csv` | Fig_6f |
| `PMDA_Approved.csv` | Fig_6f |
| `FDA-EMA-PMDA_Approved.csv` | Fig_6c, Fig_6f |

Raw experimental data (qPCR ΔCt, DNA content, C-peptide, 2-NBDG) must be provided as `.xlsx` files — see individual script headers for expected filenames.

---

## Citation

If you use these scripts, please cite:

> [Manuscript citation to be added upon acceptance]

---

*Scripts generated: 2026-05-04 | Revision 1 compiled: 2026-05-09 | Revision 2 additions: 2026-06-29*
