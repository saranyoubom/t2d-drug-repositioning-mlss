###############################################################################
#
# FIGURE 6b – CONSENSUS HUB GENE NETWORK WITH CORRELATION STRENGTH VISUALIZATION
#
# REVISED: Edge Width = Correlation Strength | Linetype = Sign (Solid/Dotted)
#
# NPG COLOR PALETTE (Nature Publishing Group)
#
# Uses Fig_4d & Fig_4e output CSV files from SBA submission folder
#
# Publication-Ready for Research in Veterinary Science
#
###############################################################################

rm(list = ls())

# Run with the working directory set to a folder containing this round's
# corrected Fig_4d_consensus_pairs.csv plus Fig_4e_Hes1_significant_correlations_only.csv
# from the R_scripts_Rev1 pipeline (see README.md, Input Data Files Required).

set.seed(42)

#==============================================================================
# 0. PACKAGE INSTALLATION & LOADING
#==============================================================================

required_packages <- c(
  "tidyverse",      # Data manipulation
  "igraph",         # Network analysis
  "ggraph",         # Network visualization with ggplot2
  "scales",         # Scaling functions
  "ggsci"           # Nature Publishing Group (NPG) color palettes
)

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}))

# Extract NPG color palette
npg_colors <- pal_npg("nrc")(10)

cat("\n✓ NPG Color Palette loaded\n")
cat(sprintf("  Colors: %s\n\n", paste(npg_colors[1:6], collapse = " | ")))

#==============================================================================
# 1. DATA LOADING - FROM FIGURE 4 CSV OUTPUTS
#==============================================================================

# Load consensus gene pairs (one column: Consensus_Pairs like 'Wnt5a_Hey1')
common_pairs_raw <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)

# Split Consensus_Pairs into Gene1 and Gene2
common_pairs <- common_pairs_raw %>%
  tidyr::separate(Consensus_Pairs, into = c("Gene1", "Gene2"), sep = "_", remove = FALSE)

# Load Hes1-normalized significant correlation matrix (wide format)
corr_mat <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv",
                     stringsAsFactors = FALSE,
                     check.names = FALSE)

# First column is gene names (row labels)
gene_names <- corr_mat[[1]]
corr_mat <- corr_mat[, -1]
rownames(corr_mat) <- gene_names

cat("═══════════════════════════════════════════════════════════════════\n")
cat("FIGURE 6b: CONSENSUS HUB GENE NETWORK WITH CORRELATION STRENGTH\n")
cat("NPG Color Palette - Publication Ready\n")
cat("Data Source: Fig_4d & Fig_4e\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Dataset loaded:\n")
cat(sprintf(" Consensus pairs: %d\n", nrow(common_pairs)))
cat(sprintf(" Correlation genes: %d\n\n", ncol(corr_mat)))
cat(sprintf(" Unique genes in consensus list: %d\n\n",
            length(unique(c(common_pairs$Gene1, common_pairs$Gene2)))))

#==============================================================================
# 2. MERGE CORRELATION DATA INTO EDGE LIST
#==============================================================================

# Convert wide correlation matrix to long edge list
corr_long <- as.data.frame(as.matrix(corr_mat)) %>%
  tibble::rownames_to_column(var = "Gene1") %>%
  tidyr::pivot_longer(
    cols = -Gene1,
    names_to = "Gene2",
    values_to = "r"
  )

# Remove self-correlations and NA values
corr_long <- corr_long %>% 
  dplyr::filter(Gene1 != Gene2, !is.na(r))

# Build symmetric pair key (order-independent)
corr_long <- corr_long %>%
  dplyr::mutate(
    pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")
  ) %>%
  # Remove duplicates (keep only unique pairs)
  dplyr::distinct(pair_key, .keep_all = TRUE)

# Add pair_key to consensus pairs
common_pairs <- common_pairs %>%
  dplyr::mutate(
    pair_key = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")
  )

# Merge: keep only consensus pairs, bring in r
common_pairs <- common_pairs %>%
  dplyr::left_join(
    corr_long %>% dplyr::select(pair_key, r),
    by = "pair_key"
  )

# Optional: assign p_value = NA (not available directly from this CSV)
common_pairs$p_value <- NA_real_

# Classify edges: Positive vs Negative
common_pairs$edge_type <- ifelse(common_pairs$r > 0, "Positive", "Negative")

cat("Correlation Summary:\n")
cat(sprintf(" Positive consensus correlations: %d (r > 0)\n", sum(common_pairs$r > 0, na.rm = TRUE)))
cat(sprintf(" Negative consensus correlations: %d (r < 0)\n", sum(common_pairs$r < 0, na.rm = TRUE)))
cat(sprintf(" Mean positive r: %.4f\n", mean(common_pairs$r[common_pairs$r > 0], na.rm = TRUE)))
cat(sprintf(" Mean negative r: %.4f\n", mean(common_pairs$r[common_pairs$r < 0], na.rm = TRUE)))
cat(sprintf(" Range |r|: %.4f to %.4f\n\n",
            min(abs(common_pairs$r), na.rm = TRUE),
            max(abs(common_pairs$r), na.rm = TRUE)))

#==============================================================================
# 3. GENE FUNCTIONAL CLASSIFICATION - NPG COLOR MAPPED
#==============================================================================

# Get all unique genes from consensus pairs
all_genes <- unique(c(common_pairs$Gene1, common_pairs$Gene2))

cat("All genes in consensus network:\n")
print(sort(all_genes))
cat("\n")

# Define functional categories
calcium_signaling <- c("Cacna1c", "Cacna1d", "Itpr1", "Ryr1", "Ryr2", "Ryr3")
incretin <- c("Glp1r")
wnt_signaling <- c("Hey1", "Tcf7", "Tcf7l2", "Wnt5a", "Wnt5b", "Wnt9a")
metabolic <- c("Glut2", "Ptbp1")

# Create metadata for ONLY the genes that appear in consensus pairs
gene_metadata <- data.frame(
  Gene = all_genes,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    Category = dplyr::case_when(
      Gene %in% calcium_signaling ~ "Calcium Signaling",
      Gene %in% incretin ~ "Incretin Signaling",
      Gene %in% wnt_signaling ~ "Wnt Signaling",
      Gene %in% metabolic ~ "Metabolic Regulation",
      TRUE ~ "Other"
    )
  )

cat("Gene metadata (consensus network genes):\n")
print(gene_metadata)
cat("\n")

#==============================================================================
# 4. EXPORT FUNCTION (PDF + PNG, 600 DPI)
#==============================================================================

export_figure <- function(filename_base, plot_object, width_cm = 7.7, height_cm = 6) {

  ggsave(
    paste0(filename_base, ".svg"),
    plot_object,
    scale = 2.3,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    paste0(filename_base, ".png"),
    plot_object,
    scale = 2.3,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 600,
    bg = "white"
  )

  cat(sprintf("✓ Exported: %s (.svg & .png)\n", filename_base))

}

#==============================================================================
# 5. NETWORK CONSTRUCTION WITH CORRELATION ATTRIBUTES
#==============================================================================

# Create igraph object from edge list
g_consensus <- graph_from_data_frame(
  d = common_pairs %>% dplyr::select(Gene1, Gene2, r, edge_type, p_value),
  vertices = gene_metadata,
  directed = FALSE
)

# Add edge attributes
E(g_consensus)$correlation <- common_pairs$r
E(g_consensus)$abs_correlation <- abs(common_pairs$r)
E(g_consensus)$edge_type <- common_pairs$edge_type

# Calculate node metrics
V(g_consensus)$degree <- degree(g_consensus)
V(g_consensus)$betweenness <- betweenness(g_consensus, normalized = TRUE)
V(g_consensus)$closeness <- closeness(g_consensus, normalized = TRUE)
V(g_consensus)$eigen_centrality <- eigen_centrality(g_consensus)$vector

# Identify hub genes (degree >= 6 connections)
hub_threshold <- 6
V(g_consensus)$is_hub <- ifelse(V(g_consensus)$degree >= hub_threshold, "Hub", "Non-hub")

# Print network statistics
cat("\nNetwork Statistics:\n")
cat(sprintf(" Nodes: %d\n", vcount(g_consensus)))
cat(sprintf(" Edges: %d\n", ecount(g_consensus)))
cat(sprintf(" Network density: %.3f\n", edge_density(g_consensus)))
cat(sprintf(" Average degree: %.2f\n", mean(V(g_consensus)$degree)))
cat(sprintf(" Hub genes (degree ≥%d): %d\n\n", hub_threshold, sum(V(g_consensus)$is_hub == "Hub")))

# Print gene statistics
gene_stats <- data.frame(
  Gene = V(g_consensus)$name,
  Degree = V(g_consensus)$degree,
  Betweenness = round(V(g_consensus)$betweenness, 4),
  Closeness = round(V(g_consensus)$closeness, 4),
  Eigenvector = round(V(g_consensus)$eigen_centrality, 4),
  Category = V(g_consensus)$Category,
  Hub_Status = V(g_consensus)$is_hub
) %>%
  arrange(desc(Degree))

cat("Top Hub Genes (ranked by degree centrality):\n")
print(gene_stats[1:nrow(gene_stats), ], row.names = FALSE)
cat("\n")

#==============================================================================
# 6. EDGE WIDTH SCALING BASED ON ABSOLUTE CORRELATION STRENGTH
#==============================================================================

# Scale edge widths: min |r| = 0.5, max |r| = 3.0
edge_width_range <- c(0.5, 3.0)
min_corr <- min(E(g_consensus)$abs_correlation)
max_corr <- max(E(g_consensus)$abs_correlation)

E(g_consensus)$width <- edge_width_range[1] +
  (E(g_consensus)$abs_correlation - min_corr) / (max_corr - min_corr) *
  (edge_width_range[2] - edge_width_range[1])

cat(sprintf("Edge Width Scaling:\n"))
cat(sprintf(" |r| range: %.4f to %.4f\n", min_corr, max_corr))
cat(sprintf(" Width range: %.2f to %.2f (pt)\n\n", edge_width_range[1], edge_width_range[2]))

#==============================================================================
# 7. EDGE LINETYPE ASSIGNMENT: SOLID (Positive) vs DOTTED (Negative)
#==============================================================================

# Define edge linetypes based on correlation sign
edge_linetypes <- ifelse(E(g_consensus)$edge_type == "Positive", "solid", "dotted")

cat(sprintf("Edge Linetypes:\n"))
cat(sprintf(" Positive correlations (r > 0): SOLID [%d edges]\n",
            sum(E(g_consensus)$edge_type == "Positive")))
cat(sprintf(" Negative correlations (r < 0): DOTTED [%d edges]\n\n",
            sum(E(g_consensus)$edge_type == "Negative")))

#==============================================================================
# 8. NPG COLOR PALETTE ASSIGNMENT - FUNCTIONAL CATEGORIES
#==============================================================================

# NPG Palette Assignments (6 distinct colors)
# npg_colors[1] = Red/Orange
# npg_colors[2] = Green
# npg_colors[3] = Blue (light)
# npg_colors[4] = Blue (dark)
# npg_colors[5] = Yellow
# npg_colors[6] = Purple

npg_category_colors <- c(
  "Calcium Signaling" = npg_colors[4],        # Dark Blue
  "Incretin Signaling" = npg_colors[1],       # Red/Orange
  "Wnt Signaling" = npg_colors[2],            # Green
  "Metabolic Regulation" = npg_colors[3]      # Light Blue
)

cat("NPG Color Assignments:\n")
cat(sprintf(" Calcium Signaling: %s (npg_colors[4])\n", npg_colors[4]))
cat(sprintf(" Incretin Signaling: %s (npg_colors[1])\n", npg_colors[1]))
cat(sprintf(" Wnt Signaling: %s (npg_colors[2])\n", npg_colors[2]))
cat(sprintf(" Metabolic Regulation: %s (npg_colors[3])\n\n", npg_colors[3]))

#==============================================================================
# 9. CIRCULAR LAYOUT NETWORK VISUALIZATION WITH NPG COLORS
#==============================================================================

# Create the network plot with NPG palette
p_network <- ggraph(g_consensus, layout = "circle") +

  # ENHANCED EDGES: Width = |correlation|, Linetype = sign, Color = grey70
  geom_edge_link(
    aes(
      width = width,
      linetype = edge_type
    ),
    color = "grey70", # All edges grey70
    alpha = 0.6,
    show.legend = TRUE
  ) +

  # Edge linetype scale: Solid for positive, Dotted for negative
  scale_edge_linetype_manual(
    name = "Correlation",
    values = c("Positive" = "solid", "Negative" = "dotted"),
    guide = guide_legend(order = 3, override.aes = list(color = "black", width = 1))
  ) +

  # Edge width scale
  scale_edge_width_continuous(
    name = "Correlation Strength",
    range = edge_width_range,
    breaks = c(0.5, 1.5, 2.5, 3.0),
    labels = c("Weak", "Moderate", "Strong", "Very Strong"),
    guide = guide_legend(order = 4, override.aes = list(color = "black"))
  ) +

  # Node coloring - NPG PALETTE
  scale_color_manual(
    name = "Functional Category",
    values = npg_category_colors,
    guide = guide_legend(order = 1)
  ) +

  scale_fill_manual(
    name = "Functional Category",
    values = npg_category_colors,
    guide = "none"
  ) +

  # Node size = degree centrality
  geom_node_point(
    aes(
      size = degree,
      color = Category,
      fill = Category
    ),
    shape = 21,
    stroke = 1.2,
    alpha = 0.85
  ) +

  # Node labels
  geom_node_text(
    aes(
      label = name,
      fontface = ifelse(degree >= hub_threshold, "bold.italic", "italic")
    ),
    repel = TRUE,
    size = 5,
    max.overlaps = 25,
    family = "Arial"
  ) +

  # Node size scale
  scale_size_continuous(
    name = "Degree Centrality",
    range = c(5, 14),
    breaks = c(2, 4, 6),
    limits = c(min(gene_stats$Degree), max(gene_stats$Degree)),
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

  labs(
    title = "Consensus Hub Gene Network"
  )

# print(p_network)  # disabled for headless Rscript run (14-Sep-2026 Rev3 regeneration) --
                     # default device lacks Arial in its PostScript font-metric database;
                     # actual outputs come from ggsave/export_figure below, unaffected.

export_figure("Fig_6b_ConsensusHubNetwork_NPG",
              p_network, width_cm = 7.7, height_cm = 6)

#==============================================================================
# 10. EDGE STATISTICS EXPORT
#==============================================================================

edge_stats <- as_data_frame(g_consensus, what = "edges") %>%
  arrange(desc(abs(correlation))) %>%
  select(from, to, correlation, abs_correlation, edge_type, width) %>%
  rename(
    Gene1 = from,
    Gene2 = to,
    `Correlation (r)` = correlation,
    `|Correlation|` = abs_correlation,
    `Edge Type` = edge_type,
    `Edge Width (pt)` = width
  )

write.csv(edge_stats, "Fig_6b_Edge_Correlation_Statistics.csv", row.names = FALSE)

cat("\nTop 20 Edges by Correlation Strength:\n")
print(edge_stats[1:min(20, nrow(edge_stats)), ], row.names = FALSE)
cat("\n")

#==============================================================================
# 11. EXPORT GENE STATISTICS
#==============================================================================

write.csv(gene_stats, "Fig_6b_Gene_Centrality_Summary.csv", row.names = FALSE)

#==============================================================================
# 12. LEGEND EXPLANATION - NPG PALETTE
#==============================================================================

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat(" ✅ FIGURE 6b COMPLETE - NPG COLOR PALETTE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("VISUALIZATION GUIDE:\n")
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("NODE SIZE (Degree Centrality):\n")
cat(" • Larger circles = More gene connections\n")
cat(" • Ranges: Low (2) → Medium (4) → High (6+)\n")
cat(" • Bold text = Hub genes (≥5 connections)\n\n")

cat("EDGE WIDTH (Correlation Strength):\n")
cat(" • Thin edges: |r| ≈ 0.6-0.7 (weak correlation)\n")
cat(" • Medium edges: |r| ≈ 0.8-0.9 (moderate correlation)\n")
cat(" • Thick edges: |r| ≈ 0.95-1.0 (very strong correlation)\n\n")

cat("EDGE LINETYPE (Correlation Sign):\n")
cat(" • SOLID line: Positive correlation (r > 0)\n")
cat(" → Genes move together; co-regulated\n")
cat(" • DOTTED line: Negative correlation (r < 0)\n")
cat(" → Genes move oppositely; antagonistic\n\n")

cat("EDGE COLOR:\n")
cat(" • All edges: GREY70 (uniform color)\n\n")

cat("NODE COLOR (Functional Category - NPG Palette):\n")
cat(sprintf(" • DARK BLUE (%s): Calcium Signaling\n", npg_colors[4]))
cat("   (Cacna1c, Cacna1d, Itpr1, Ryr1-3)\n\n")

cat(sprintf(" • RED/ORANGE (%s): Incretin Signaling\n", npg_colors[1]))
cat("   (Glp1r)\n\n")

cat(sprintf(" • GREEN (%s): Wnt Signaling\n", npg_colors[2]))
cat("   (Tcf7, Tcf7l2, Wnt*, Hey1)\n\n")

cat(sprintf(" • LIGHT BLUE (%s): Metabolic Regulation\n", npg_colors[3]))
cat("   (Glut2, Ptbp1)\n\n")

cat("═══════════════════════════════════════════════════════════════════\n")
cat(" DATA LINEAGE:\n")
cat(" ────────────────────────────────────────────────────────────────\n")
cat(" Input files:\n")
cat(sprintf("   • Fig_4d_consensus_pairs.csv (%d consensus pairs)\n", nrow(common_pairs)))
cat("   • Fig_4e_Hes1_significant_correlations_only.csv (Hes1-norm correlations)\n\n")
cat(" Output files generated:\n")
cat("   • Fig_6b_ConsensusHubNetwork_NPG.svg\n")
cat("   • Fig_6b_ConsensusHubNetwork_NPG.png\n")
cat("   • Fig_6b_Edge_Correlation_Statistics.csv\n")
cat("   • Fig_6b_Gene_Centrality_Summary.csv\n\n")
cat(" Color Palette:\n")
cat("   • NPG (Nature Publishing Group) via ggsci::pal_npg()\n\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# End of script