################################################################################
# FIGURE 5 - GENE INTERACTION NETWORK & STRING DATABASE NOVELTY ANALYSIS
# 4-Panel Comprehensive Analysis:
#   Panel A: Gene Interaction Network (STRING-validated + novel discoveries)
#   Panel B: STRING Database Coverage Distribution
#   Panel C: Correlation Strength Distribution by Category
#   Panel D: Top 15 Novel Gene Pair Discoveries
# Author: Revised for publication
# Date: January 2026
# Color Palette: NPG (Nature Publishing Group)
# Font: Arial
# Output: Individual panels saved as PNG (600 DPI) and SVG
################################################################################

rm(list = ls())
setwd("C:/Users/saran/OneDrive/VSCBIC/C2F-2025/Publication/Submission/CMBL_Sub2_291026/Analysis_outputs/R_17genes")

cat("\n", rep("=", 80), "\n", sep="")
cat("FIGURE 5: GENE INTERACTION NETWORK & STRING DATABASE NOVELTY DISCOVERY\n")
cat("Publication-Ready with NPG Color Palette\n")
cat("4-Panel Comprehensive Analysis\n")
cat(rep("=", 80), "\n\n", sep="")

################################################################################
# LOAD REQUIRED LIBRARIES
################################################################################

cat("Loading packages...\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(igraph)
  library(ggraph)
  library(patchwork)
  library(ggsci)
})

# Extract NPG color palette
npg_colors <- pal_npg("nrc")(10)

cat("✓ Packages loaded\n\n")

################################################################################
# CREATE OUTPUT DIRECTORIES
################################################################################

dir.create("./figures_publication", showWarnings = FALSE, recursive = TRUE)
dir.create("./figures_publication/data_export", showWarnings = FALSE, recursive = TRUE)

cat("✓ Output directories created\n")
cat("  → Figures: ./figures_publication/\n")
cat("  → Data export: ./figures_publication/data_export/\n\n")

################################################################################
# HELPER FUNCTION: CREATE CANONICAL PAIR NAME
################################################################################

# Always sort gene names alphabetically to create consistent pair identifiers
create_pair_name <- function(gene1, gene2) {
  pairs <- mapply(function(g1, g2) {
    sorted <- sort(c(g1, g2))
    paste(sorted, collapse = "-")
  }, gene1, gene2, SIMPLIFY = TRUE)
  return(pairs)
}

################################################################################
# LOAD DATA
################################################################################

cat("Loading data files...\n")

# Load STRING interactions
string_data <- read.csv("string_interactions.csv", stringsAsFactors = FALSE) %>%
  mutate(Pair = create_pair_name(node1, node2)) %>%
  select(Pair, node1, node2, combined_score) %>%
  rename(STRING_Score = combined_score) %>%
  distinct(Pair, .keep_all = TRUE)

cat(sprintf("✓ STRING interactions: %d pairs\n", nrow(string_data)))

# Load 4-method validated pairs (consensus from all normalization methods)
# These use underscore separator: Gene1_Gene2
validated_pairs_raw <- read.csv("Fig_4d_consensus_pairs.csv", stringsAsFactors = FALSE)

validated_pairs <- validated_pairs_raw %>%
  mutate(
    # Split the underscore-separated pair names
    Gene1 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 1),
    Gene2 = sapply(strsplit(Consensus_Pairs, "_"), `[`, 2),
    # Create canonical pair name (alphabetically sorted with hyphen)
    Pair = create_pair_name(Gene1, Gene2)
  ) %>%
  select(Pair, Gene1, Gene2) %>%
  distinct(Pair, .keep_all = TRUE)

cat(sprintf("✓ Validated pairs (4-method consensus): %d pairs\n", nrow(validated_pairs)))

# Load Hes1-normalized significant correlations (correlation matrix format)
correlation_data <- read.csv("Fig_4e_Hes1_significant_correlations_only.csv", 
                             stringsAsFactors = FALSE, row.names = 1) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene1") %>%
  pivot_longer(-Gene1, names_to = "Gene2", values_to = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  filter(Gene1 < Gene2) %>%  # Keep only upper triangle to avoid duplicates
  mutate(
    Pair = create_pair_name(Gene1, Gene2),
    r = Correlation
  ) %>%
  select(Pair, Gene1, Gene2, r)

cat(sprintf("✓ Correlation data: %d significant pairs (p < 0.05)\n\n", nrow(correlation_data)))

################################################################################
# CLASSIFY GENE PAIRS BY STRING STATUS & VALIDATION
################################################################################

cat("Classifying gene pairs...\n")

# Create sets for fast lookup
string_set <- string_data$Pair
validated_set <- validated_pairs$Pair

# Filter strong correlations (|r| ≥ 0.7) for network analysis
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

# Get novel discoveries (top candidates for new biological knowledge)
novel_discoveries <- all_pairs %>%
  filter(Category == "Novel Discovery") %>%
  arrange(desc(abs(r)))

# Summary statistics by category
summary_stats <- all_pairs %>%
  group_by(Category) %>%
  summarise(
    Count = n(),
    Mean_r = mean(abs(r)),
    Median_r = median(abs(r)),
    Min_r = min(abs(r)),
    Max_r = max(abs(r)),
    Strong_r_0.9 = sum(abs(r) >= 0.9),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

cat(rep("-", 80), "\n", sep="")
cat("PAIR CLASSIFICATION SUMMARY\n")
cat(rep("-", 80), "\n", sep="")
print(summary_stats)
cat("\n")

cat(sprintf("✓ Total pairs analyzed (|r| ≥ 0.7): %d\n", nrow(all_pairs)))
cat(sprintf("✓ Novel discoveries: %d (%.1f%%)\n", 
            nrow(novel_discoveries), 
            100 * nrow(novel_discoveries) / nrow(all_pairs)))
cat(sprintf("✓ Known validated (STRING + 4-method): %d\n",
            sum(all_pairs$Category == "Known Validated")))
cat(sprintf("✓ STRING only (not in 4-method): %d\n",
            sum(all_pairs$Category == "STRING Only")))
cat(sprintf("✓ Not validated (not in STRING, not in 4-method): %d\n",
            sum(all_pairs$Category == "Not Validated")))
cat("\n")

################################################################################
# CUSTOM THEME FOR FIGURE 5
################################################################################

theme_fig5 <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = "Arial"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", family = "Arial"),
    axis.title = element_text(size = 10, face = "bold", color = "black", family = "Arial"),
    axis.text = element_text(size = 10, color = "black", family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 10, face = "bold", family = "Arial"),
    panel.grid.major = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )
}

################################################################################
# PANEL A: GENE INTERACTION NETWORK
################################################################################

cat("Creating Panel A: Gene Interaction Network...\n")

# Select top pairs for network visualization
# Prioritize STRING-validated and novel discoveries
string_validated <- all_pairs %>%
  filter(Category == "Known Validated") %>%
  arrange(desc(abs(r))) %>%
  head(30)

novel_top <- novel_discoveries %>%
  head(20)

# If we don't have enough pairs, supplement with STRING Only and Not Validated
remaining_needed <- 50 - nrow(string_validated) - nrow(novel_top)

if (remaining_needed > 0) {
  supplement <- all_pairs %>%
    filter(Category %in% c("STRING Only", "Not Validated")) %>%
    arrange(desc(abs(r))) %>%
    head(remaining_needed)
} else {
  supplement <- data.frame()
}

cat(sprintf("Network composition:\n"))
cat(sprintf("  STRING-validated: %d pairs\n", nrow(string_validated)))
cat(sprintf("  Novel discoveries: %d pairs\n", nrow(novel_top)))
cat(sprintf("  Supplemental (STRING Only/Not Validated): %d pairs\n", nrow(supplement)))

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

# Create graph object
g <- graph_from_data_frame(edges, directed = FALSE)

# Calculate centrality metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, weights = E(g)$weight, directed = FALSE)

cat(sprintf("Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

# Identify novel genes
novel_genes <- unique(c(novel_top$Gene1, novel_top$Gene2))
V(g)$node_type <- ifelse(V(g)$name %in% novel_genes, "Novel", "Known")

# Set edge colors - NPG palette
edge_colors <- case_when(
  E(g)$EdgeType == "Novel Discovery" ~ npg_colors[1],      
  E(g)$EdgeType == "STRING Validated" ~ npg_colors[4],     
  E(g)$EdgeType == "STRING Only" ~ npg_colors[2],          
  TRUE ~ "gray70"                                          
)
E(g)$color <- edge_colors
E(g)$width <- E(g)$weight * 1.5

# Set layout
set.seed(42)
layout_net <- layout_with_fr(g, weights = E(g)$weight)

# Create network plot
# ---- display text for the categories (manuscript wording); the category codes themselves are unchanged ----
LAB_CAT <- c("Known Validated" = "Known,
analytically consistent",
             "Novel Discovery" = "Analytically consistent,\ndatabase-unregistered",
             "STRING Only"     = "STRING-only",
             "Not Validated"   = "Not analytically consistent")
LAB_EDGE <- c("STRING Validated" = "Known, analytically consistent", "Novel Discovery" = LAB_CAT[["Novel Discovery"]],
              "STRING Only" = LAB_CAT[["STRING Only"]], "Not Validated" = LAB_CAT[["Not Validated"]])

panel_a <- ggraph(g, layout = layout_net) +
  geom_edge_link(aes(edge_color = EdgeType, edge_width = weight),
                 alpha = 0.6,show.legend = TRUE) +
  geom_node_point(aes(size = degree, color = node_type, fill = node_type),
                  alpha = 0.8, shape = 21) +
  geom_node_text(aes(label = name, color = node_type),
                 size = 2.5, repel = TRUE, max.overlaps = 20,
                 fontface = "italic", family = "Arial") +
  scale_edge_color_manual(
    name = "Pair category",
    labels = LAB_EDGE,
    values = c("STRING Validated" = npg_colors[4], 
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  scale_edge_width_continuous(name = expression("Abs. Correlation (|" * italic(r) * "|)"), range = c(0.3, 2)) +
  scale_color_manual(
    name = "Gene status",
    labels = c("Novel" = "Database-unregistered", "Known" = "Known"),
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_fill_manual(
    name = "Gene status",
    labels = c("Novel" = "Database-unregistered", "Known" = "Known"),
    values = c("Novel" = npg_colors[1], "Known" = npg_colors[4])
  ) +
  scale_size_continuous(name = "Degree", range = c(3, 10)) +
  labs(title = "Co-expression network",
       #subtitle = sprintf("%d STRING-validated + %d novel discoveries",
                          #nrow(string_validated), nrow(novel_top))
) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "right",
    legend.key.size = unit(0.42, "cm"),     # the category names are longer than before: keep the legend inside the canvas
    legend.key.height = unit(0.42, "cm"),
    legend.key.width = unit(0.42, "cm"),
    legend.text = element_text(size = 7, family = "Arial"),
    legend.title = element_text(size = 7.5, face = "bold", family = "Arial"),
    legend.spacing.y = unit(0.12, "cm"),
    legend.margin = margin(2, 0, 2, 0),
    plot.margin = margin(5, 5, 5, 5)
  )

panel_a

# Panel A - Network
ggsave("./figures_publication/Fig5a_Network.png",
       plot = panel_a,
       scale = 1.5,
       width = 8, height = 7.5, 
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5a_Network.svg",
       plot = panel_a,
       scale = 1.5,
       width = 8, height = 7.5, 
       dpi = 600, units = "cm", bg = "white")
cat("✓ Panel A created\n\n")

################################################################################
# PANEL B: STRING COVERAGE DISTRIBUTION
################################################################################

cat("Creating Panel B: STRING Database Coverage...\n")

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
    labels = LAB_CAT,
    values = c("Known Validated" = npg_colors[4],      # Green
               "Novel Discovery" = npg_colors[1],      # Red
               "STRING Only" = npg_colors[2],          # Yellow
               "Not Validated" = "gray70")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +   # two rows: the category names are long
  labs(title = "STRING Database Coverage",
       #subtitle = sprintf("n = %d gene pairs (|r| ≥ 0.7)", nrow(all_pairs))
       ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = "Arial"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, family = "Arial"),
    legend.position = "bottom",
    legend.text = element_text(size = 8, family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    plot.margin = margin(5, 5, 5, 5)
  )

panel_b

# Panel B - Coverage
ggsave("./figures_publication/Fig5b_Coverage.png",
       plot = panel_b,
       scale = 1.5,
       width = 8, height = 7.5, 
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5b_Coverage.svg",
       plot = panel_b,
       scale = 1.5,
       width = 8, height = 7.5, 
       dpi = 600, units = "cm", bg = "white")


# Export coverage summary
write.csv(panel_b_data, 
          "./figures_publication/data_export/Fig5B_Coverage_Summary.csv", 
          row.names = FALSE)

cat("✓ Panel B created\n")
cat("  → Exported: Fig5B_Coverage_Summary.csv\n\n")

################################################################################
# PANEL C: CORRELATION STRENGTH DISTRIBUTION
################################################################################

cat("Creating Panel C: Correlation Distribution...\n")

panel_c <- ggplot(all_pairs, aes(x = abs(r), fill = Category)) +
  geom_histogram(bins = 15, alpha = 0.7, position = "identity", 
                 color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black", 
             linewidth = 0.4, alpha = 0.5) +
  facet_wrap(~factor(Category, 
                     levels = c("Known Validated", "Novel Discovery", 
                               "STRING Only", "Not Validated")),
             ncol = 2, scales = "free_y", labeller = as_labeller(LAB_CAT)) +
  scale_fill_manual(
    values = c("Known Validated" = npg_colors[4],
               "Novel Discovery" = npg_colors[1],
               "STRING Only" = npg_colors[2],
               "Not Validated" = "gray70")
  ) +
  labs(title = "Correlation Strength Distribution",
       x = expression("Absolute Correlation (|" * italic(r) * "|)"),
       y = "Frequency") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    axis.title = element_text(size = 10, family = "Arial", color = "black")
  )

panel_c

# Panel C - Distribution
ggsave("./figures_publication/Fig5c_Distribution.png",
       plot = panel_c,
       scale = 1.5,
       width = 8, height = 8.5, 
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5c_Distribution.svg",
       plot = panel_c,
       scale = 1.5,
       width = 8, height = 8.5, 
       dpi = 600, units = "cm", bg = "white")

# Export all classified pairs
write.csv(all_pairs %>% select(Pair, Gene1, Gene2, r, Category) %>% 
            arrange(desc(abs(r))), 
          "./figures_publication/data_export/Fig5C_All_Pairs_Classified.csv", 
          row.names = FALSE)

cat("✓ Panel C created\n")
cat("  → Exported: Fig5C_All_Pairs_Classified.csv\n\n")

################################################################################
# PANEL D: TOP NOVEL DISCOVERIES
################################################################################

cat("Creating Panel D: Top Novel Discoveries...\n")

# Determine how many novel discoveries we have
n_novel <- min(15, nrow(novel_discoveries))

if (n_novel > 0) {
  novel_top_n <- novel_discoveries %>%
    head(n_novel) %>%
    mutate(
      PairLabel = paste0(Gene1, "–", Gene2),
      PairLabel = factor(PairLabel, levels = rev(PairLabel))
    )
  
  panel_d <- ggplot(novel_top_n, aes(x = PairLabel, y = abs(r), fill = abs(r))) +
    geom_col(color = "white", linewidth = 0.3) +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray40", 
               linewidth = 0.3, alpha = 0.6) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray20", 
               linewidth = 0.3, alpha = 0.6) +
    geom_text(aes(label = sprintf("%.3f", abs(r))), 
              hjust = -0.1, size = 2.5, fontface = "bold", family = "Arial") +
    coord_flip() +
    scale_fill_gradient(low = npg_colors[5], high = npg_colors[1],
                       name = expression("|" * italic(r) * "|"),
                       limits = c(0.930951812154159, 0.996076286491127)) +
    scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
    labs(title = sprintf("Top %d analytically consistent,\ndatabase-unregistered pairs", n_novel),
         #subtitle = "Not in STRING but validated by 4-method consensus",
         x = "",
         y = expression("Absolute Correlation (|" * italic(r) * "|)")) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 12, family = "Arial"),
      axis.text.y = element_text(size = 8, family = "Arial", face = "italic", color = "black"),
      axis.text.x = element_text(size = 8, family = "Arial", color = "black"),
      legend.position = "none"
    )
    
  # Export top novel discoveries
  write.csv(novel_top_n %>% select(Gene1, Gene2, r) %>%
              arrange(desc(abs(r))),
            "./figures_publication/data_export/Fig5D_Novel_Discoveries.csv",
            row.names = FALSE)
  
  cat(sprintf("✓ Panel D created (showing %d novel discoveries)\n", n_novel))
  cat("  → Exported: Fig5D_Novel_Discoveries.csv\n\n")
  
} else {
  # If no novel discoveries, create placeholder
  panel_d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No novel discoveries found\n(all 4-method validated pairs are in STRING)",
             size = 5, family = "Arial") +
    labs(title = "D. Novel Gene Pair Discoveries",
         subtitle = "Not in STRING but validated by 4-method consensus") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10, family = "Arial"),
      plot.subtitle = element_text(hjust = 0.5, size = 7, family = "Arial")
    )
  
  cat("⚠ Panel D: No novel discoveries to display\n")
  cat("  → All 4-method validated pairs are already in STRING database\n\n")
}

panel_d

# Panel D - Novel Discoveries
ggsave("./figures_publication/Fig5d_Novel_Discoveries.png",
       plot = panel_d,
       scale = 1.5,
       width = 8, height = 8.5, 
       dpi = 600, units = "cm", bg = "white")

ggsave("./figures_publication/Fig5d_Novel_Discoveries.svg",
       plot = panel_d,
       scale = 1.5,
       width = 8, height = 8.5, 
       dpi = 600, units = "cm", bg = "white")

################################################################################
# NETWORK ANALYSIS SUMMARY
################################################################################

cat(rep("=", 80), "\n", sep="")
cat("NETWORK ANALYSIS SUMMARY\n")
cat(rep("=", 80), "\n\n", sep="")

is_connected <- is_connected(g)
n_components <- components(g)$no

cat("Network Connectivity:\n")
cat(sprintf("  Connected: %s\n", ifelse(is_connected, "Yes", "No")))
cat(sprintf("  Components: %d\n", n_components))
cat(sprintf("  Nodes: %d genes\n", vcount(g)))
cat(sprintf("  Edges: %d interactions\n", ecount(g)))
cat(sprintf("  Density: %.3f\n", edge_density(g)))
cat(sprintf("  Average degree: %.2f\n", mean(degree(g))))
cat("\n")

cat("Gene Pair Classification:\n")
print(summary_stats)
cat("\n")

if (vcount(g) > 0) {
  cat("Hub Genes (Top 5 by connectivity):\n")
  node_data <- data.frame(
    Gene = V(g)$name,
    Degree = V(g)$degree,
    Betweenness = V(g)$betweenness,
    NodeType = V(g)$node_type,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Degree)) %>%
    head(5)
  
  print(node_data)
  cat("\n")
}

if (nrow(novel_discoveries) > 0) {
  cat("Novel Discovery Statistics:\n")
  novel_stats <- novel_discoveries %>%
    summarise(
      Count = n(),
      Mean_r = mean(abs(r)),
      Median_r = median(abs(r)),
      Strong_r_0.9 = sum(abs(r) >= 0.9),
      Very_Strong_r_0.95 = sum(abs(r) >= 0.95)
    )
  
  cat(sprintf("  Total novel pairs: %d\n", novel_stats$Count))
  cat(sprintf("  Mean |r|: %.3f\n", novel_stats$Mean_r))
  cat(sprintf("  Median |r|: %.3f\n", novel_stats$Median_r))
  cat(sprintf("  Pairs with |r| ≥ 0.9: %d (%.1f%%)\n", 
              novel_stats$Strong_r_0.9, 
              100 * novel_stats$Strong_r_0.9 / novel_stats$Count))
  cat(sprintf("  Pairs with |r| ≥ 0.95: %d (%.1f%%)\n",
              novel_stats$Very_Strong_r_0.95,
              100 * novel_stats$Very_Strong_r_0.95 / novel_stats$Count))
} else {
  cat("Novel Discovery Statistics:\n")
  cat("  No novel discoveries found\n")
  cat("  → All 4-method validated pairs are already documented in STRING\n")
}
cat("\n")

################################################################################
# EXPORT NETWORK DATA
################################################################################

cat("Exporting network analysis data...\n")

if (vcount(g) > 0) {
  # Node properties
  all_node_data <- data.frame(
    Gene = V(g)$name,
    Degree = V(g)$degree,
    Betweenness = V(g)$betweenness,
    NodeType = V(g)$node_type,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Degree))
  
  write.csv(all_node_data,
            "./figures_publication/data_export/Fig5_Network_Nodes.csv",
            row.names = FALSE)
  
  # Edge list
  edge_data <- edges %>%
    mutate(abs_corr = weight) %>%
    arrange(desc(abs_corr))
  
  write.csv(edge_data %>% select(from, to, weight, EdgeType),
            "./figures_publication/data_export/Fig5_Network_Edges.csv",
            row.names = FALSE)
  
  cat("✓ Fig5_Network_Nodes.csv\n")
  cat("✓ Fig5_Network_Edges.csv\n")
}

# Summary statistics
write.csv(summary_stats,
          "./figures_publication/data_export/Fig5_Summary_Statistics.csv",
          row.names = FALSE)

cat("✓ Fig5_Summary_Statistics.csv\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat(rep("=", 80), "\n", sep="")
cat("FIGURE 5 COMPLETE - PUBLICATION READY\n")
cat(rep("=", 80), "\n\n", sep="")

cat("INDIVIDUAL PANEL FIGURES SAVED:\n")
cat("  Format: PNG (600 DPI) + SVG\n")
cat("  Dimensions: 9.16 × 9.16 cm each\n\n")
cat("  ✓ Fig5a_Network.png / .svg\n")
cat(sprintf("      → %d nodes, %d edges\n", vcount(g), ecount(g)))
cat(sprintf("      → %d STRING-validated, %d novel discoveries\n\n", 
            nrow(string_validated), nrow(novel_top)))
cat("  ✓ Fig5b_Coverage.png / .svg\n")
cat(sprintf("      → %d pairs classified into 4 categories\n\n", nrow(all_pairs)))
cat("  ✓ Fig5c_Distribution.png / .svg\n")
cat("      → Correlation strength histograms by category\n\n")
cat("  ✓ Fig5d_Novel_Discoveries.png / .svg\n")
cat(sprintf("      → Top %d novel gene pairs ranked by |r|\n\n", n_novel))

cat("DATA EXPORTS:\n")
cat("  ✓ Fig5B_Coverage_Summary.csv\n")
cat("  ✓ Fig5C_All_Pairs_Classified.csv\n")
if (n_novel > 0) {
  cat("  ✓ Fig5D_Novel_Discoveries.csv\n")
}
cat("  ✓ Fig5_Network_Nodes.csv\n")
cat("  ✓ Fig5_Network_Edges.csv\n")
cat("  ✓ Fig5_Summary_Statistics.csv\n\n")

cat("All panels saved individually. Ready for manuscript submission.\n\n")

################################################################################
# END OF SCRIPT
################################################################################
