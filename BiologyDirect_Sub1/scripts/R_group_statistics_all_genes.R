# Group statistics for every RT-qPCR gene (fold change, 2^-ddCt), total DNA, and C-peptide at 5.5 mM,
# recomputed from the replicate values with the procedure stated in the Methods:
#   Shapiro-Wilk within each treatment group; all p > 0.05 -> one-way ANOVA with Bonferroni pairwise
#   comparisons (pooled-SD t-tests, pairwise.t.test); otherwise Kruskal-Wallis with Bonferroni pairwise
#   Wilcoxon rank-sum (pairwise.wilcox.test; exact when there are no ties, normal approximation with
#   continuity correction when ties exist).
# Input : input_data/replicates_long.csv (gene, group, rep, value; n = 4 technical replicates per group)
# Output: R_group_statistics_all_genes.csv
# Usage : Rscript R_group_statistics_all_genes.R [input_csv] [output_csv]
args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "input_data/replicates_long.csv"
outp <- if (length(args) >= 2) args[2] else "R_group_statistics_all_genes.csv"
d <- read.csv(inp, stringsAsFactors = FALSE)
G <- c("CTRL", "DAPT", "DKK-1")
d$group <- factor(d$group, levels = G)
sem <- function(x) sd(x) / sqrt(length(x))
pairs <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))
res <- list()
for (g in unique(d$gene)) {
  s <- d[d$gene == g, ]
  sh <- sapply(G, function(x) shapiro.test(s$value[s$group == x])$p.value)
  param <- all(sh > 0.05)
  if (param) {
    omn <- summary(aov(value ~ group, data = s))[[1]][["Pr(>F)"]][1]
    pw  <- pairwise.t.test(s$value, s$group, p.adjust.method = "none")$p.value      # pooled SD, raw
    test <- "One-way ANOVA + Bonferroni pooled-SD t"
  } else {
    omn <- kruskal.test(value ~ group, data = s)$p.value
    pw  <- suppressWarnings(pairwise.wilcox.test(s$value, s$group, p.adjust.method = "none"))$p.value
    test <- "Kruskal-Wallis + Bonferroni Wilcoxon"
  }
  getp <- function(a, b) { v <- if (!is.na(pw[b, a])) pw[b, a] else pw[a, b]; as.numeric(v) }
  row <- data.frame(dataset = g, test = test, omnibus_p = omn,
                    shapiro_p_CTRL = sh[1], shapiro_p_DAPT = sh[2], shapiro_p_DKK1 = sh[3],
                    mean_CTRL = mean(s$value[s$group == "CTRL"]), sem_CTRL = sem(s$value[s$group == "CTRL"]),
                    mean_DAPT = mean(s$value[s$group == "DAPT"]), sem_DAPT = sem(s$value[s$group == "DAPT"]),
                    mean_DKK1 = mean(s$value[s$group == "DKK-1"]), sem_DKK1 = sem(s$value[s$group == "DKK-1"]),
                    stringsAsFactors = FALSE)
  for (p in pairs) {
    raw <- getp(p[1], p[2]); nm <- paste0(p[1], "_vs_", gsub("-", "", p[2]))
    row[[paste0("p_raw_", nm)]] <- raw
    row[[paste0("p_bonf_", nm)]] <- min(1, raw * 3)
  }
  res[[g]] <- row
}
out <- do.call(rbind, res)
write.csv(out, outp, row.names = FALSE)
cat("Wrote", nrow(out), "datasets to", outp, "\n")
cat(R.version.string, "\n")
