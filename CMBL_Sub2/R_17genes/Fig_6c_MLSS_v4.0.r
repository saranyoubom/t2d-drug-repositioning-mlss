###############################################################################
#
# MLSS v4.0: NETWORK PHARMACOLOGY WITH POTENCY_MULTIPLIER WEIGHTING
#
# Full integration of drug-binding potency (pIC50) into synergy and antagonism
# calculations via Potency_Multiplier = max(pIC50_i, pIC50_j) / 9
#
# REVISED FORMULA:
# MLSS_v4.0 = [C×0.45 + B×0.30 + V×0.10 + P×0.15] × 9
#           + Σ(r_pos_ij × Potency_Multiplier_ij) × 1
#           - Σ(|r_neg_ij| × Potency_Multiplier_ij) × 1
#
# Key Changes:
# - Base Score multiplier: 10 → 9
# - Synergy scale factor: 5 → 1
# - Antagonism scale factor: 5 → 1
# - Result: More conservative scoring with larger base score contribution
#
# Uses: Fig_6b_Gene_Centrality_Summary.csv
#       Fig_4e_Hes1_significant_correlations_only.csv
#       drug.target.interaction.csv
#
# Output: Publication-ready figures, CSVs with potency weighting
#
# v4.0 REVISED - January 2026
#
###############################################################################

rm(list = ls())

# ============================================================================
# 0. CONFIGURATION & SETUP
# ============================================================================

setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")
set.seed(42)

required_packages <- c("dplyr", "tidyr", "ggplot2", "ggsci", "tibble")

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

npg_colors <- pal_npg("nrc")(10)

cat("\n")
cat("████████████████████████████████████████████████████████████████████\n")
cat(" MLSS v4.0 REVISED - NETWORK + POTENCY_MULTIPLIER WEIGHTING\n")
cat(" Base × 9 | Synergy/Antagonism × 1\n")
cat("████████████████████████████████████████████████████████████████████\n\n")

# ============================================================================
# 1. LOAD DATA FILES
# ============================================================================

cat("STEP 1: Loading data files...\n")

hub_data <- read.csv("Fig_6b_Gene_Centrality_Summary.csv", stringsAsFactors = FALSE)
drug_data <- read.csv("drug.target.interaction.csv", stringsAsFactors = FALSE)
corr_mat_raw <- read.csv(
  "Fig_4e_Hes1_significant_correlations_only.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_names <- corr_mat_raw[[1]]
corr_mat <- corr_mat_raw[, -1]
rownames(corr_mat) <- toupper(gene_names)
colnames(corr_mat) <- toupper(gene_names)

cat(sprintf(" ✓ Genes loaded: %d\n", nrow(hub_data)))
cat(sprintf(" ✓ Drug interactions loaded: %d\n", nrow(drug_data)))
cat(sprintf(" ✓ Correlation matrix: %d × %d genes\n", nrow(corr_mat), ncol(corr_mat)))
cat(" ✓ Gene names standardized to UPPERCASE\n\n")

# ============================================================================
# 2. EXTRACT AND PREPARE POTENCY DATA (pIC50)
# ============================================================================

cat("STEP 2: Extracting potency data (pIC50) from drug.target.interaction...\n")

potency_raw <- drug_data %>%
  dplyr::filter(
    ACT_TYPE %in% c("IC50", "Ki", "EC50"),
    ORGANISM == "Homo sapiens",
    RELATION == "=",
    !is.na(ACT_VALUE)
  ) %>%
  dplyr::select(DRUG_NAME, GENE, ACT_VALUE, ACT_TYPE) %>%
  dplyr::mutate(Gene_Standardized = toupper(GENE))

potency_lookup <- potency_raw %>%
  dplyr::group_by(DRUG_NAME, Gene_Standardized) %>%
  dplyr::summarise(
    pIC50 = max(ACT_VALUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::distinct(DRUG_NAME, Gene_Standardized, pIC50)

cat(sprintf(" ✓ Potency pairs extracted: %d (Drug, Gene) → pIC50\n", nrow(potency_lookup)))
cat(sprintf(" ✓ pIC50 range: %.2f - %.2f\n", 
            min(potency_lookup$pIC50, na.rm = TRUE),
            max(potency_lookup$pIC50, na.rm = TRUE)))
cat(sprintf(" ✓ Unique drugs with potency data: %d\n",
            length(unique(potency_lookup$DRUG_NAME))))
cat(sprintf(" ✓ Unique genes with potency data: %d\n\n",
            length(unique(potency_lookup$Gene_Standardized))))

# ============================================================================
# 3. FILTER FOR HUB GENES (Degree >= 5; was >= 6 before Hey1 exclusion)
# ============================================================================

cat("STEP 3: Identifying hub genes (Degree >= 5)...\n")

hub_genes_df <- hub_data %>%
  dplyr::filter(Degree >= 5)

cat(sprintf(" ✓ Hub genes identified: %d\n", nrow(hub_genes_df)))
cat("\n Hub Genes:\n")
for (i in 1:nrow(hub_genes_df)) {
  cat(sprintf("   %d. %-10s (Degree: %2d, %s)\n",
              i,
              hub_genes_df$Gene[i],
              hub_genes_df$Degree[i],
              hub_genes_df$Category[i]))
}

hub_genes <- unique(toupper(hub_genes_df$Gene))
cat(sprintf("\n ✓ Using %d hub genes for analysis\n\n", length(hub_genes)))

# ============================================================================
# 4. CONVERT CORRELATION MATRIX TO PAIRWISE FORMAT
# ============================================================================

cat("STEP 4: Converting correlation matrix to pairwise format...\n")

correlations <- as.data.frame(as.matrix(corr_mat)) %>%
  tibble::rownames_to_column(var = "Gene1") %>%
  tidyr::pivot_longer(
    cols = -Gene1,
    names_to = "Gene2",
    values_to = "r"
  ) %>%
  dplyr::filter(Gene1 != Gene2, !is.na(r)) %>%
  dplyr::mutate(
    pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")
  ) %>%
  dplyr::distinct(pair_key, .keep_all = TRUE)

correlations_slim <- correlations %>% dplyr::select(pair_key, r)

cat(sprintf(" ✓ Correlation pairs: %d\n", nrow(correlations)))
cat(sprintf(" ✓ Correlation range: %.4f to %.4f\n\n",
            min(correlations$r), max(correlations$r)))

# ============================================================================
# 5. FILTER DRUG-TARGET INTERACTIONS FOR HUB GENES
# ============================================================================

cat("STEP 5: Filtering drug-target interactions (hub genes only)...\n")

drug_targets <- drug_data %>%
  dplyr::filter(toupper(GENE) %in% hub_genes, ORGANISM == "Homo sapiens") %>%
  dplyr::select(DRUG_NAME, GENE, TARGET_NAME, ACTION_TYPE) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Gene_Standardized = toupper(GENE))

unique_drugs <- sort(unique(drug_targets$DRUG_NAME))

cat(sprintf(" ✓ Drug-target pairs (hub genes): %d\n", nrow(drug_targets)))
cat(sprintf(" ✓ Unique drugs: %d\n", length(unique_drugs)))
cat(sprintf(" ✓ Unique hub genes targeted: %d\n\n", 
            length(unique(drug_targets$Gene_Standardized))))

# ============================================================================
# 6. DEFINE FUNCTIONAL CATEGORIES
# ============================================================================

cat("STEP 6: Defining functional categories...\n\n")

hub_data_cat <- hub_genes_df %>%
  dplyr::mutate(Gene_Upper = toupper(Gene))

for (cat_name in unique(hub_data_cat$Category)) {
  genes_in_cat <- hub_data_cat %>%
    dplyr::filter(Category == cat_name) %>%
    dplyr::pull(Gene)
  cat(sprintf(" %-20s: %s\n", cat_name, paste(genes_in_cat, collapse = ", ")))
}

cat("\n")

# ============================================================================
# 7. GENERATE ALL POSSIBLE DRUG COMBINATIONS
# ============================================================================

cat("STEP 7: Generating all drug combinations...\n")

all_drug_combos <- expand.grid(
  Drug_A = unique_drugs,
  Drug_B = unique_drugs,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(Drug_A < Drug_B) %>%
  dplyr::mutate(
    Combination = paste(Drug_A, "+", Drug_B),
    Combo_ID = dplyr::row_number()
  )

total_combos <- nrow(all_drug_combos)

cat(sprintf(" ✓ Total combinations: %d\n", total_combos))
cat(sprintf(" ✓ Expected C(%d, 2): %d\n\n",
            length(unique_drugs),
            length(unique_drugs) * (length(unique_drugs) - 1) / 2))

# ============================================================================
# 8. ENRICH COMBINATIONS WITH TARGETS AND PATHWAYS
# ============================================================================

cat("STEP 8: Enriching combinations with targets and pathways...\n")

get_drug_targets <- function(drug_name) {
  targets <- drug_targets %>%
    dplyr::filter(DRUG_NAME == drug_name) %>%
    dplyr::pull(Gene_Standardized) %>%
    unique()
  return(targets)
}

all_drug_combos <- all_drug_combos %>%
  dplyr::mutate(
    Targets_A = sapply(Drug_A, function(d) paste(get_drug_targets(d), collapse = "|")),
    Targets_B = sapply(Drug_B, function(d) paste(get_drug_targets(d), collapse = "|")),
    n_targets_A = sapply(Drug_A, function(d) length(get_drug_targets(d))),
    n_targets_B = sapply(Drug_B, function(d) length(get_drug_targets(d))),
    n_targets_combined = n_targets_A + n_targets_B,
    Pathway_A = sapply(Drug_A, function(d) {
      genes <- get_drug_targets(d)
      pathways <- unique(hub_data_cat$Category[hub_data_cat$Gene_Upper %in% genes])
      paste(pathways, collapse = "|")
    }),
    Pathway_B = sapply(Drug_B, function(d) {
      genes <- get_drug_targets(d)
      pathways <- unique(hub_data_cat$Category[hub_data_cat$Gene_Upper %in% genes])
      paste(pathways, collapse = "|")
    }),
    Cross_Pathway = ifelse(Pathway_A != Pathway_B, "Yes", "No")
  )

cat(sprintf(" ✓ Enrichment complete: %d combinations\n\n", nrow(all_drug_combos)))

# ============================================================================
# 9. HELPER FUNCTIONS: GET pIC50 AND POTENCY_MULTIPLIER
# ============================================================================

cat("STEP 9: Defining potency lookup and Potency_Multiplier functions...\n")

get_pIC50 <- function(drug_name, gene) {
  gene_upper <- toupper(gene)
  pIC50 <- potency_lookup %>%
    dplyr::filter(DRUG_NAME == drug_name, Gene_Standardized == gene_upper) %>%
    dplyr::pull(pIC50)
  
  if (length(pIC50) == 0) return(NA)
  return(pIC50)
}

calculate_potency_multiplier <- function(drug_a, gene_i, drug_b, gene_j) {
  pIC50_i <- get_pIC50(drug_a, gene_i)
  pIC50_j <- get_pIC50(drug_b, gene_j)
  
  if (is.na(pIC50_i) || is.na(pIC50_j)) return(0)
  
  max_pIC50 <- max(pIC50_i, pIC50_j)
  potency_mult <- max_pIC50 / 9
  
  return(min(potency_mult, 1.0))
}

get_mean_pIC50 <- function(drug_name, targets) {
  targets_vec <- unique(toupper(targets))
  pIC50_values <- sapply(targets_vec, function(g) get_pIC50(drug_name, g))
  mean_pIC50 <- mean(pIC50_values[!is.na(pIC50_values)], na.rm = TRUE)
  return(if (is.nan(mean_pIC50)) 0 else mean_pIC50)
}

cat(" ✓ get_pIC50() function defined\n")
cat(" ✓ calculate_potency_multiplier() function defined\n")
cat(" ✓ get_mean_pIC50() function defined\n\n")

# ============================================================================
# 10. CORRELATION STATS AND COMPLEMENTARITY/BALANCE/COVERAGE
# ============================================================================

cat("STEP 10: Calculating correlation statistics (C, B, V)...\n")

get_corr_stats <- function(genes_A, genes_B, correlations_df) {
  genes_A <- unique(toupper(genes_A))
  genes_B <- unique(toupper(genes_B))
  
  if (length(genes_A) == 0 || length(genes_B) == 0) {
    return(data.frame(
      Mean_Pos_r = 0,
      Max_Pos_r = 0,
      Mean_Neg_r = 0,
      Sum_Neg_Abs_r = 0,
      Mean_Abs_r = 0
    ))
  }
  
  pair_df <- expand.grid(
    Gene1 = genes_A,
    Gene2 = genes_B,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")
    ) %>%
    dplyr::left_join(correlations_df, by = "pair_key") %>%
    dplyr::filter(!is.na(r))
  
  if (nrow(pair_df) == 0) {
    return(data.frame(
      Mean_Pos_r = 0,
      Max_Pos_r = 0,
      Mean_Neg_r = 0,
      Sum_Neg_Abs_r = 0,
      Mean_Abs_r = 0
    ))
  }
  
  pos_r <- pair_df$r[pair_df$r > 0]
  neg_r <- pair_df$r[pair_df$r < 0]
  
  mean_pos_r <- if (length(pos_r) > 0) mean(pos_r) else 0
  mean_neg_r <- if (length(neg_r) > 0) mean(neg_r) else 0
  sum_neg_abs_r <- if (length(neg_r) > 0) sum(abs(neg_r)) else 0
  mean_abs_r <- mean(abs(pair_df$r))
  
  return(data.frame(
    Mean_Pos_r = mean_pos_r,
    Max_Pos_r = if (length(pos_r) > 0) max(pos_r) else 0,
    Mean_Neg_r = mean_neg_r,
    Sum_Neg_Abs_r = sum_neg_abs_r,
    Mean_Abs_r = mean_abs_r
  ))
}

all_drug_combos <- all_drug_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    targets_A_set = list(if (nzchar(Targets_A)) unlist(strsplit(Targets_A, "\\|")) else character(0)),
    targets_B_set = list(if (nzchar(Targets_B)) unlist(strsplit(Targets_B, "\\|")) else character(0)),
    
    intersection_size = length(intersect(targets_A_set, targets_B_set)),
    union_size = length(union(targets_A_set, targets_B_set)),
    Complementarity_C = if (union_size > 0) 1 - (intersection_size / union_size) else 0,
    
    Balance_B = if ((n_targets_A + n_targets_B) > 0) {
      1 - abs(n_targets_A - n_targets_B) / (n_targets_A + n_targets_B)
    } else {
      0
    },
    
    corr_stats = list(get_corr_stats(targets_A_set, targets_B_set, correlations_slim))
  ) %>%
  tidyr::unnest_wider(corr_stats) %>%
  dplyr::ungroup()

mean_abs_r_min <- min(all_drug_combos$Mean_Abs_r)
mean_abs_r_max <- max(all_drug_combos$Mean_Abs_r)

all_drug_combos <- all_drug_combos %>%
  dplyr::mutate(
    Network_Coverage_V = if (mean_abs_r_max > mean_abs_r_min) {
      (Mean_Abs_r - mean_abs_r_min) / (mean_abs_r_max - mean_abs_r_min)
    } else {
      0
    }
  )

cat(" ✓ Complementarity (C) calculated\n")
cat(" ✓ Balance (B) calculated\n")
cat(" ✓ Network Coverage (V) normalized\n")
cat(sprintf(" ✓ Network Coverage range: %.4f - %.4f\n\n",
            min(all_drug_combos$Network_Coverage_V),
            max(all_drug_combos$Network_Coverage_V)))

# ============================================================================
# 11. CALCULATE POTENCY COMPONENT (P)
# ============================================================================

cat("STEP 11: Calculating Potency component (P)...\n")

all_drug_combos <- all_drug_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Mean_pIC50_A = get_mean_pIC50(Drug_A, targets_A_set),
    Mean_pIC50_B = get_mean_pIC50(Drug_B, targets_B_set),
    Mean_pIC50 = mean(c(Mean_pIC50_A, Mean_pIC50_B), na.rm = TRUE),
    Potency_P = min(max(Mean_pIC50 - 4, 0), 5) / 5
  ) %>%
  dplyr::ungroup()

cat(sprintf(" ✓ Mean pIC50 calculated for all combinations\n"))
cat(sprintf(" ✓ Potency component (P) range: %.3f - %.3f\n",
            min(all_drug_combos$Potency_P, na.rm = TRUE),
            max(all_drug_combos$Potency_P, na.rm = TRUE)))
cat(sprintf(" ✓ Mean Potency score: %.3f\n\n",
            mean(all_drug_combos$Potency_P, na.rm = TRUE)))

# ============================================================================
# 12. SYNERGY WITH POTENCY_MULTIPLIER (SCALE = 1)
# ============================================================================

cat("STEP 12: Calculating Synergy bonus with Potency_Multiplier (scale × 1)...\n")

synergy_bonus_with_potency <- function(drug_a, targets_a, drug_b, targets_b, 
                                       correlations_df) {
  targets_a <- unique(toupper(targets_a))
  targets_b <- unique(toupper(targets_b))
  
  synergy_total <- 0
  
  for (i in seq_along(targets_a)) {
    for (j in seq_along(targets_b)) {
      gene_i <- targets_a[i]
      gene_j <- targets_b[j]
      
      pair_key <- paste(pmin(gene_i, gene_j), pmax(gene_i, gene_j), sep = "_")
      r_ij <- correlations_df %>%
        dplyr::filter(pair_key == !!pair_key) %>%
        dplyr::pull(r)
      
      if (length(r_ij) == 0 || is.na(r_ij)) next
      
      if (r_ij > 0) {
        potency_mult <- calculate_potency_multiplier(drug_a, gene_i, drug_b, gene_j)
        synergy_ij <- r_ij * potency_mult * 1  # SCALE = 1 (changed from 5)
        synergy_total <- synergy_total + synergy_ij
      }
    }
  }
  
  return(synergy_total)
}

all_drug_combos <- all_drug_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Synergy_Bonus = synergy_bonus_with_potency(
      Drug_A, targets_A_set,
      Drug_B, targets_B_set,
      correlations_slim
    )
  ) %>%
  dplyr::ungroup()

cat(sprintf(" ✓ Synergy bonus (with Potency_Multiplier × 1) calculated\n"))
cat(sprintf(" ✓ Synergy range: %.2f - %.2f\n",
            min(all_drug_combos$Synergy_Bonus, na.rm = TRUE),
            max(all_drug_combos$Synergy_Bonus, na.rm = TRUE)))
cat(sprintf(" ✓ Combinations with synergy > 0: %d\n\n",
            sum(all_drug_combos$Synergy_Bonus > 0)))

# ============================================================================
# 13. ANTAGONISM WITH POTENCY_MULTIPLIER (SCALE = 1)
# ============================================================================

cat("STEP 13: Calculating Antagonism penalty with Potency_Multiplier (scale × 1)...\n")

antagonism_penalty_with_potency <- function(drug_a, targets_a, drug_b, targets_b, 
                                            correlations_df) {
  targets_a <- unique(toupper(targets_a))
  targets_b <- unique(toupper(targets_b))
  
  penalty_total <- 0
  
  for (i in seq_along(targets_a)) {
    for (j in seq_along(targets_b)) {
      gene_i <- targets_a[i]
      gene_j <- targets_b[j]
      
      pair_key <- paste(pmin(gene_i, gene_j), pmax(gene_i, gene_j), sep = "_")
      r_ij <- correlations_df %>%
        dplyr::filter(pair_key == !!pair_key) %>%
        dplyr::pull(r)
      
      if (length(r_ij) == 0 || is.na(r_ij)) next
      
      if (r_ij < 0) {
        potency_mult <- calculate_potency_multiplier(drug_a, gene_i, drug_b, gene_j)
        penalty_ij <- abs(r_ij) * potency_mult * 1  # SCALE = 1 (changed from 5)
        penalty_total <- penalty_total + penalty_ij
      }
    }
  }
  
  return(penalty_total)
}

all_drug_combos <- all_drug_combos %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Antagonism_Penalty = antagonism_penalty_with_potency(
      Drug_A, targets_A_set,
      Drug_B, targets_B_set,
      correlations_slim
    )
  ) %>%
  dplyr::ungroup()

cat(sprintf(" ✓ Antagonism penalty (with Potency_Multiplier × 1) calculated\n"))
cat(sprintf(" ✓ Antagonism range: %.2f - %.2f\n",
            min(all_drug_combos$Antagonism_Penalty, na.rm = TRUE),
            max(all_drug_combos$Antagonism_Penalty, na.rm = TRUE)))
cat(sprintf(" ✓ Combinations with antagonism > 0: %d\n\n",
            sum(all_drug_combos$Antagonism_Penalty > 0)))

# ============================================================================
# 14. FINAL MLSS v4.0 CALCULATION (BASE × 9)
# ============================================================================

cat("STEP 14: Calculating final MLSS v4.0 scores (Base × 9)...\n")

all_drug_combos <- all_drug_combos %>%
  dplyr::mutate(
    # Base Score = (0.45·C + 0.30·B + 0.10·V + 0.15·P) × 9 (changed from × 10)
    Base_Score = (0.45 * Complementarity_C + 
                  0.30 * Balance_B + 
                  0.10 * Network_Coverage_V + 
                  0.15 * Potency_P) * 9,
    
    # Final MLSS v4.0 = Base_Score + Synergy - Antagonism
    MLSS_v4.0 = Base_Score + Synergy_Bonus - Antagonism_Penalty
  ) %>%
  dplyr::arrange(dplyr::desc(MLSS_v4.0)) %>%
  dplyr::mutate(Rank = dplyr::row_number())

cat(sprintf(" ✓ Base Score range: %.2f - %.2f\n",
            min(all_drug_combos$Base_Score, na.rm = TRUE),
            max(all_drug_combos$Base_Score, na.rm = TRUE)))
cat(sprintf(" ✓ MLSS v4.0 range: %.2f - %.2f\n",
            min(all_drug_combos$MLSS_v4.0, na.rm = TRUE),
            max(all_drug_combos$MLSS_v4.0, na.rm = TRUE)))
cat(sprintf(" ✓ Mean MLSS v4.0: %.2f\n\n",
            mean(all_drug_combos$MLSS_v4.0, na.rm = TRUE)))

# ============================================================================
# 15. SELECT Top 15 AND TOP 50
# ============================================================================

cat("STEP 15: Selecting top combinations...\n")

top_50_combos <- all_drug_combos %>% head(min(50, nrow(all_drug_combos)))
top_15_combos <- all_drug_combos %>% head(min(15, nrow(all_drug_combos)))

cat(sprintf(" ✓ Top 15 selected\n"))
cat(sprintf(" ✓ Top 50 selected\n\n"))

cat("═══════════════════════════════════════════════════════════════\n")
cat("Top 15 DRUG COMBINATIONS (MLSS v4.0 REVISED)\n")
cat("Base × 9 | Synergy/Antagonism × 1\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

print(
  top_15_combos %>%
    dplyr::select(
      Rank, Combination, Mean_pIC50,
      Complementarity_C, Balance_B, Network_Coverage_V, Potency_P,
      Synergy_Bonus, Antagonism_Penalty, MLSS_v4.0
    )
)

cat("\n")

# ============================================================================
# 16. NPG COLOR PALETTE
# ============================================================================

cat("STEP 16: Applying NPG color palette...\n\n")

npg_mlss_colors <- c(
  "Complementarity (45%)" = npg_colors[1],
  "Balance (30%)" = npg_colors[4],
  "Coverage (10%)" = npg_colors[3],
  "Potency (15%)" = npg_colors[2],
  "Synergy Bonus" = npg_colors[5],
  "Antagonism Penalty" = npg_colors[6]
)

cat("NPG Color Assignments:\n")
cat(sprintf(" Complementarity (45%%): %s\n", npg_colors[1]))
cat(sprintf(" Balance (30%%): %s\n", npg_colors[4]))
cat(sprintf(" Coverage (10%%): %s\n", npg_colors[3]))
cat(sprintf(" Potency (15%%): %s\n", npg_colors[2]))
cat(sprintf(" Synergy Bonus: %s\n", npg_colors[5]))
cat(sprintf(" Antagonism Penalty: %s\n\n", npg_colors[6]))

# ============================================================================
# 17. VISUALIZATION: Top 15
# ============================================================================

cat("STEP 17: Creating visualization (Top 15)...\n")

component_data <- top_15_combos %>%
  dplyr::arrange(dplyr::desc(MLSS_v4.0)) %>%
  dplyr::mutate(
    C_Contribution = Complementarity_C * 0.45 * 9,
    B_Contribution = Balance_B * 0.30 * 9,
    V_Contribution = Network_Coverage_V * 0.10 * 9,
    P_Contribution = Potency_P * 0.15 * 9,
    Synergy_Component = Synergy_Bonus,
    Penalty_Component = Antagonism_Penalty,
    Combination = factor(Combination, levels = rev(Combination))
  ) %>%
  dplyr::select(
    Combination, C_Contribution, B_Contribution, V_Contribution, P_Contribution,
    Synergy_Component, Penalty_Component, MLSS_v4.0
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      C_Contribution, B_Contribution, V_Contribution, P_Contribution,
      Synergy_Component, Penalty_Component
    ),
    names_to = "Component",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Component = factor(
      Component,
      levels = c(
        "C_Contribution", "B_Contribution", "V_Contribution", "P_Contribution",
        "Synergy_Component", "Penalty_Component"
      ),
      labels = c(
        "Complementarity (45%)",
        "Balance (30%)",
        "Coverage (10%)",
        "Potency (15%)",
        "Synergy Bonus",
        "Antagonism Penalty"
      )
    )
  )

mlss_labels <- component_data %>%
  dplyr::group_by(Combination, MLSS_v4.0) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

p <- ggplot(component_data, aes(x = Combination, y = Score, fill = Component)) +
  geom_bar(
    stat = "identity",
    width = 0.85,
    color = "white",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey50",
    linewidth = 0.5,
    linetype = "solid"
  ) +
  geom_text(
    data = mlss_labels,
    aes(
      x = Combination,
      y = MLSS_v4.0,
      label = sprintf("%.2f", MLSS_v4.0)
    ),
    hjust = -0.3,
    size = 3.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Score\nComponent",
    values = npg_mlss_colors
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  guides(fill = guide_legend(nrow = 3)) +
  coord_flip() +
  labs(
    title = "MLSS Component Breakdown of Top 15 Drug Combinations",
    x = NULL,
    y = "Weighted Score Contribution"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 1,
      color = "black",
      family = "Arial"
    ),
    axis.text.y = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    axis.title.x = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12, family = "Arial", hjust = 1, face = "bold"),
    legend.text = element_text(size = 12, family = "Arial"),
    legend.box = "vertical"
  )

print(p)

ggsave(
  "Fig_6c_MLSS_top15.png",
  plot = p,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Fig_6c_MLSS_top15.svg",
  plot = p,
  width = 7.7,
  height = 6,
  units = "cm",
  scale = 2.3,
  dpi = 600,
  bg = "white"
)

cat(" ✓ Saved: Fig_6c_MLSS_top15.png (600 DPI)\n")
cat(" ✓ Saved: Fig_6c_MLSS_top15.pdf (600 DPI)\n\n")

# ============================================================================
# 18. EXPORT RESULTS
# ============================================================================

cat("STEP 18: Exporting results to CSV...\n")

write.csv(
  all_drug_combos %>%
    dplyr::select(
      Rank, Combination, Drug_A, Drug_B,
      Pathway_A, Pathway_B, Cross_Pathway,
      n_targets_A, n_targets_B, n_targets_combined,
      Targets_A, Targets_B,
      Mean_pIC50_A, Mean_pIC50_B, Mean_pIC50,
      Complementarity_C, Balance_B, Network_Coverage_V, Potency_P,
      Mean_Pos_r, Mean_Neg_r, Mean_Abs_r,
      Base_Score, Synergy_Bonus, Antagonism_Penalty,
      MLSS_v4.0
    ),
  "Fig_6c_MLSS_ALL_Combinations.csv",
  row.names = FALSE
)

cat(sprintf(" ✓ Saved: Fig_6c_MLSS_ALL_Combinations.csv (%d rows)\n",
            nrow(all_drug_combos)))

write.csv(
  top_50_combos %>%
    dplyr::select(
      Rank, Combination, Drug_A, Drug_B,
      Pathway_A, Pathway_B, Cross_Pathway,
      n_targets_combined,
      Mean_pIC50,
      Complementarity_C, Balance_B, Network_Coverage_V, Potency_P,
      Mean_Pos_r, Mean_Neg_r,
      Base_Score, Synergy_Bonus, Antagonism_Penalty,
      MLSS_v4.0
    ),
  "Fig_6c_MLSS_Top50.csv",
  row.names = FALSE
)

cat(sprintf(" ✓ Saved: Fig_6c_MLSS_Top50.csv (%d rows)\n\n",
            nrow(top_50_combos)))

# ============================================================================
# 19. SUMMARY STATISTICS
# ============================================================================

cat("█████████████████████████████████████████████████████████████████\n")
cat(" MLSS v4.0 REVISED - COMPLETE\n")
cat(" Base × 9 | Synergy/Antagonism × 1\n")
cat("█████████████████████████████████████████████████████████████████\n\n")

cat("POTENCY DATA SUMMARY:\n")
cat(sprintf(" • Potency pairs loaded: %d\n", nrow(potency_lookup)))
cat(sprintf(" • pIC50 range: %.2f - %.2f\n",
            min(potency_lookup$pIC50, na.rm = TRUE),
            max(potency_lookup$pIC50, na.rm = TRUE)))
cat(sprintf(" • Drugs with potency data: %d\n", length(unique(potency_lookup$DRUG_NAME))))
cat(sprintf(" • Hub genes with potency data: %d/%d\n",
            length(intersect(unique(potency_lookup$Gene_Standardized), hub_genes)),
            length(hub_genes)))

cat("\nCOMBINATION STATISTICS:\n")
cat(sprintf(" • Total combinations: %d\n", total_combos))
cat(sprintf(" • Hub genes used: %d\n", length(hub_genes)))
cat(sprintf(" • Cross-pathway: %d (%.1f%%)\n",
            sum(all_drug_combos$Cross_Pathway == "Yes"),
            100 * sum(all_drug_combos$Cross_Pathway == "Yes") / total_combos))
cat(sprintf(" • Same-pathway: %d (%.1f%%)\n",
            sum(all_drug_combos$Cross_Pathway == "No"),
            100 * sum(all_drug_combos$Cross_Pathway == "No") / total_combos))

cat("\nSCORE COMPONENT RANGES:\n")
cat(sprintf(" • Complementarity (C): %.3f - %.3f\n",
            min(all_drug_combos$Complementarity_C),
            max(all_drug_combos$Complementarity_C)))
cat(sprintf(" • Balance (B): %.3f - %.3f\n",
            min(all_drug_combos$Balance_B),
            max(all_drug_combos$Balance_B)))
cat(sprintf(" • Coverage (V): %.3f - %.3f\n",
            min(all_drug_combos$Network_Coverage_V),
            max(all_drug_combos$Network_Coverage_V)))
cat(sprintf(" • Potency (P): %.3f - %.3f\n",
            min(all_drug_combos$Potency_P, na.rm = TRUE),
            max(all_drug_combos$Potency_P, na.rm = TRUE)))
cat(sprintf(" • Mean pIC50: %.2f - %.2f\n",
            min(all_drug_combos$Mean_pIC50, na.rm = TRUE),
            max(all_drug_combos$Mean_pIC50, na.rm = TRUE)))
cat(sprintf(" • Base Score: %.2f - %.2f\n",
            min(all_drug_combos$Base_Score, na.rm = TRUE),
            max(all_drug_combos$Base_Score, na.rm = TRUE)))
cat(sprintf(" • Synergy Bonus: %.2f - %.2f\n",
            min(all_drug_combos$Synergy_Bonus, na.rm = TRUE),
            max(all_drug_combos$Synergy_Bonus, na.rm = TRUE)))
cat(sprintf(" • Antagonism Penalty: %.2f - %.2f\n",
            min(all_drug_combos$Antagonism_Penalty, na.rm = TRUE),
            max(all_drug_combos$Antagonism_Penalty, na.rm = TRUE)))
cat(sprintf(" • MLSS v4.0: %.2f - %.2f\n\n",
            min(all_drug_combos$MLSS_v4.0, na.rm = TRUE),
            max(all_drug_combos$MLSS_v4.0, na.rm = TRUE)))

cat("Top 15 CHARACTERISTICS:\n")
top15_summary <- top_15_combos %>%
  dplyr::summarise(
    Mean_C = mean(Complementarity_C),
    Mean_B = mean(Balance_B),
    Mean_V = mean(Network_Coverage_V),
    Mean_P = mean(Potency_P),
    Mean_pIC50 = mean(Mean_pIC50, na.rm = TRUE),
    Mean_r_pos = mean(Mean_Pos_r),
    Mean_r_neg = mean(Mean_Neg_r),
    Cross_Pct = 100 * sum(Cross_Pathway == "Yes") / dplyr::n(),
    With_Synergy = sum(Synergy_Bonus > 0),
    With_Antagonism = sum(Antagonism_Penalty > 0)
  )

cat(sprintf(" • Mean Complementarity: %.3f\n", top15_summary$Mean_C))
cat(sprintf(" • Mean Balance: %.3f\n", top15_summary$Mean_B))
cat(sprintf(" • Mean Coverage: %.3f\n", top15_summary$Mean_V))
cat(sprintf(" • Mean Potency: %.3f\n", top15_summary$Mean_P))
cat(sprintf(" • Mean pIC50: %.2f\n", top15_summary$Mean_pIC50))
cat(sprintf(" • Mean r_pos: %.3f\n", top15_summary$Mean_r_pos))
cat(sprintf(" • Mean r_neg: %.3f\n", top15_summary$Mean_r_neg))
cat(sprintf(" • Cross-pathway %%: %.1f%%\n", top15_summary$Cross_Pct))
cat(sprintf(" • With synergy: %d/%d\n", top15_summary$With_Synergy, nrow(top_15_combos)))
cat(sprintf(" • With antagonism: %d/%d\n\n", top15_summary$With_Antagonism, nrow(top_15_combos)))

cat("OUTPUT FILES CREATED:\n")
cat(" 1. Fig_6c_MLSS_top15.png (600 DPI)\n")
cat(" 2. Fig_6c_MLSS_top15.pdf (600 DPI)\n")
cat(" 3. Fig_6c_MLSS_ALL_Combinations.csv (ALL combos)\n")
cat(" 4. Fig_6c_MLSS_Top50.csv (TOP 50)\n\n")

cat("═════════════════════════════════════════════════════════════════\n")
cat(" ✅ MLSS v4.0 REVISED WITH POTENCY_MULTIPLIER COMPLETE\n")
cat(" ✅ Formula: Base × 9 + Synergy × 1 - Antagonism × 1\n")
cat(" ✅ Base score is now dominant contributor\n")
cat(" ✅ PUBLICATION READY\n")
cat("═════════════════════════════════════════════════════════════════\n\n")

# End of script

