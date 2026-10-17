# Extra row of Supplementary Table S4: C-peptide at 5.5 mM glucose as the RAW ELISA signal (dOD = OD450 - OD590) of the four wells per group,
# with both procedures exactly as in R_group_statistics_both_procedures.R (one-way ANOVA + Bonferroni pooled-variance t-tests; Kruskal-Wallis + Bonferroni Wilcoxon).
# The per-ug-DNA row (workbook values, group-mean DNA) is already in SuppTable4_both_procedures.csv. Reported procedure: Kruskal-Wallis (Shapiro-Wilk p < 0.05).
# Input: Cpeptide_raw_elisa_tidy.csv (elisa_raw_audit.py)   Output: SuppTable4_cpeptide_dOD_row.csv   Run from Analysis_outputs/.
tidy <- read.csv("Cpeptide_raw_elisa_tidy.csv", stringsAsFactors = FALSE)
s <- tidy[tidy$condition_mM == 5.5, ]; s$group <- factor(s$group, levels = c("CTRL", "DAPT", "DKK-1")); s$value <- s$dOD
sem <- function(x) sd(x) / sqrt(length(x))
pairs <- list(c("CTRL", "DAPT"), c("CTRL", "DKK-1"), c("DAPT", "DKK-1"))
m <- tapply(s$value, s$group, mean); e <- tapply(s$value, s$group, sem)
sw <- shapiro.test(s$value)$p.value
an  <- summary(aov(value ~ group, data = s))[[1]][["Pr(>F)"]][1]
ptt <- pairwise.t.test(s$value, s$group, p.adjust.method = "bonferroni")$p.value
kw  <- kruskal.test(value ~ group, data = s)$p.value
pww <- suppressWarnings(pairwise.wilcox.test(s$value, s$group, p.adjust.method = "bonferroni"))$p.value
get <- function(M, a, b) { v <- if (!is.na(M[b, a])) M[b, a] else M[a, b]; as.numeric(v) }
row <- data.frame(dataset = "Cpeptide_dOD_5.5mM", reported = ifelse(sw < 0.05, "Kruskal-Wallis", "ANOVA"),
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
cols <- names(read.csv("SuppTable4_both_procedures.csv", nrows = 1))
stopifnot(identical(cols, names(row)))
write.csv(row, "SuppTable4_cpeptide_dOD_row.csv", row.names = FALSE)
cat("Shapiro-Wilk p =", signif(sw, 3), "| ANOVA p =", signif(an, 3), "| KW p =", signif(kw, 3), "\n"); print(t(row))
