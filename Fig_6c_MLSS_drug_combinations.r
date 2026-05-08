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

