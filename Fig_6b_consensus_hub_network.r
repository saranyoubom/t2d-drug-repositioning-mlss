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

