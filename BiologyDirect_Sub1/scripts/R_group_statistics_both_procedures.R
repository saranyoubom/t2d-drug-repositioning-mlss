# Supplementary Table 4: for every RT-qPCR dataset (plus total DNA and C-peptide at 5.5 mM) report BOTH
#   (A) one-way ANOVA with Bonferroni-corrected pooled-variance t-tests, and
#   (B) Kruskal-Wallis with Bonferroni-corrected Wilcoxon rank-sum tests,
# and mark the procedure that was reported (Supplementary Data 2) and any pair whose significance (adjusted p < 0.05) differs.
# Input: input_data/replicates_long.csv (n = 4 technical replicates per group)   Output: SuppTable4_both_procedures.csv
d <- read.csv("input_data/replicates_long.csv", stringsAsFactors = FALSE)
G <- c("CTRL", "DAPT", "DKK-1"); d$group <- factor(d$group, levels = G)
reported_KW <- c("Hey1", "Wnt5a", "Lef1", "Cpeptide_5.5mM")       # rank-based in the reported analysis (Supplementary Data 2)
sem <- function(x) sd(x) / sqrt(length(x))
pairs <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))
res <- list()
for (g in unique(d$gene)) {
  s <- d[d$gene == g, ]
  m <- tapply(s$value, s$group, mean); e <- tapply(s$value, s$group, sem)
  an  <- summary(aov(value ~ group, data = s))[[1]][["Pr(>F)"]][1]
  ptt <- pairwise.t.test(s$value, s$group, p.adjust.method = "bonferroni")$p.value
  kw  <- kruskal.test(value ~ group, data = s)$p.value
  pww <- suppressWarnings(pairwise.wilcox.test(s$value, s$group, p.adjust.method = "bonferroni"))$p.value
  get <- function(M, a, b) { v <- if (!is.na(M[b, a])) M[b, a] else M[a, b]; as.numeric(v) }
  row <- data.frame(dataset = g, reported = ifelse(g %in% reported_KW, "Kruskal-Wallis", "ANOVA"),
                    CTRL = sprintf("%.3g +- %.2g", m["CTRL"], e["CTRL"]), DAPT = sprintf("%.3g +- %.2g", m["DAPT"], e["DAPT"]), DKK1 = sprintf("%.3g +- %.2g", m["DKK-1"], e["DKK-1"]),
                    ANOVA_p = an, KW_p = kw, stringsAsFactors = FALSE)
  differs <- c()
  for (p in pairs) {
    nm <- paste0(p[1], "_vs_", gsub("-", "", p[2]))
    a <- get(ptt, p[1], p[2]); b <- get(pww, p[1], p[2])
    row[[paste0("ANOVA_padj_", nm)]] <- a; row[[paste0("KW_padj_", nm)]] <- b
    if ((a < 0.05) != (b < 0.05)) differs <- c(differs, nm)
  }
  row$classification_differs <- if (length(differs)) paste(differs, collapse = "; ") else ""
  res[[g]] <- row
}
out <- do.call(rbind, res)
write.csv(out, "SuppTable4_both_procedures.csv", row.names = FALSE)
cat("datasets:", nrow(out), "| with a differing classification:", sum(out$classification_differs != ""), "\n")
print(out[out$classification_differs != "", c("dataset", "reported", "classification_differs")], row.names = FALSE)
cat("min attainable KW (Wilcoxon, Bonferroni) padj with n = 4 vs 4:", round(min(out$KW_padj_CTRL_vs_DAPT[out$dataset == "Hes1"]), 4), "\n")
