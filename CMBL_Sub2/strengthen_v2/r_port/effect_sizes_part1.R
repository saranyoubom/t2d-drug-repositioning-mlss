## effect_sizes_part1.R - R port of wp5_reporting/effect_sizes.py, PART 1 (lines 1-173 of the Python file): paths, constants, loaders and statistics helpers.
## Parts 2 and 3 (rt_panel ... main, from line 174) follow; the three parts are sourced in order by effect_sizes.R. Nothing random is used.
suppressMessages(library(readxl))
if (!exists("ES_HERE")) ES_HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
ROOT <- Sys.getenv("CMBL_ROOT", normalizePath(file.path(ES_HERE, "..", "..", "..")))
PKG <- file.path(ROOT, "Submission_package_JTM"); AF3 <- file.path(PKG, "Additional file 3.xlsx"); AF4 <- file.path(PKG, "Additional file 4.xlsx")
GAPDH_JSON <- file.path(ROOT, "CMBL_Sub2", "reanalysis_17genes", "gapdh_exact.json")
FIG4B <- file.path(dirname(ROOT), "SBA", "Figures", "Fig_4b_Correlation_matrix_Relative_mRNA_expression.xlsx")
LABEL <- "within-experiment, one preparation, wells are not independent biological replicates"
G5 <- c("Tcf7l2", "Kcnj11", "Cacna1d", "Cacna1c", "Ptbp1")      # DAPT values depend on the 26 August Gapdh wells
CONTRASTS <- c("DAPT", "DKK-1"); GI <- list(CTRL = 1:4, DAPT = 5:8, `DKK-1` = 9:12); GATE <- list()

## ---------------------------------------------------------------- loaders (openpyxl data_only=True: every cell keeps its own type, empty cells are NULL)
read_rows <- function(path, sheet) {
  df <- read_excel(path, sheet = sheet, col_names = FALSE, col_types = "list", .name_repair = "minimal")
  lapply(seq_len(nrow(df)), function(i) unname(lapply(seq_len(ncol(df)), function(j) { x <- df[[j]][[i]]; if (is.null(x) || (length(x) == 1 && is.na(x))) NULL else x })))
}
strip_left <- function(row) { while (length(row) && is.null(row[[1]])) row <- row[-1]; row }
cell <- function(r, i) if (i <= length(r)) r[[i]] else NULL          # 1-based cell access, NULL beyond the row
num <- function(x) { if (is.null(x)) return(NA_real_); v <- suppressWarnings(as.numeric(x)); if (length(v) != 1) NA_real_ else v }
is_str <- function(x) !is.null(x) && is.character(x)
hdr_names <- function(r) vapply(r, function(x) if (is.null(x)) NA_character_ else as.character(x), "")   # header cells as text, empty cells NA
load_af3 <- function() {
  rows <- read_rows(AF3, "Analyzed genes (17)"); h <- which(vapply(rows, function(r) is_str(cell(r, 1)) && cell(r, 1) == "Gene", TRUE))[1]; out <- list()
  for (r in rows[(h + 1):length(rows)]) {
    if (is.null(cell(r, 1)) || (is_str(cell(r, 1)) && startsWith(cell(r, 1), "Abbreviations"))) break
    g <- function(i) cell(r, i); gn <- function(a, b) vapply(a:b, function(i) num(g(i)), 1)
    out[[length(out) + 1]] <- list(gene = g(1), category = g(2), test = g(3), overall_p = as.character(g(4)),
      mean = c(CTRL = num(g(5)), DAPT = num(g(11)), `DKK-1` = num(g(17))), sem = c(CTRL = num(g(6)), DAPT = num(g(12)), `DKK-1` = num(g(18))),
      rep = list(CTRL = gn(7, 10), DAPT = gn(13, 16), `DKK-1` = gn(19, 22)), p_raw = list(DAPT = g(23), `DKK-1` = g(27), `D-K` = g(31)), p_adj = list(DAPT = g(24), `DKK-1` = g(28), `D-K` = g(32)))
  }
  out
}
af4_rows <- function(path, name) lapply(read_rows(path, name), strip_left)
## descriptive block of an Additional file 4 'stats' sheet (numbers are stored as text there): first row per group -> mean, SEM, n; with start_label only rows after the first row starting with it
stat_means <- function(rows, start_label = NULL) {
  out <- list(); active <- is.null(start_label)
  for (r in rows) {
    if (length(r) && !is.null(start_label) && is_str(r[[1]]) && startsWith(r[[1]], start_label)) active <- TRUE
    if (active && length(r) >= 4 && is_str(r[[1]]) && r[[1]] %in% c("CTRL", "DAPT", "DKK-1") && !is.na(num(r[[2]])) && !is.na(num(r[[3]])) && !is.na(num(r[[4]])) && is.null(out[[r[[1]]]]))
      out[[r[[1]]]] <- list(mean = num(r[[2]]), sem = num(r[[3]]), n = as.integer(num(r[[4]])))
  }
  out
}
stat_pairwise <- function(rows, start_label = NULL) {       # CTRL vs DAPT and CTRL vs DKK-1 rows (p_raw, p_adjusted)
  out <- list(); active <- is.null(start_label)
  for (r in rows) {
    if (length(r) && !is.null(start_label) && is_str(r[[1]]) && startsWith(r[[1]], start_label)) active <- TRUE
    if (active && length(r) >= 4 && is_str(r[[1]]) && r[[1]] == "CTRL" && is_str(r[[2]]) && r[[2]] %in% c("DAPT", "DKK-1") && !is.na(num(r[[3]])) && !is.na(num(r[[4]])) && is.null(out[[r[[2]]]]))
      out[[r[[2]]]] <- c(num(r[[3]]), num(r[[4]]))
  }
  out
}

## ---------------------------------------------------------------- statistics
sem_ <- function(a) sd(a) / sqrt(length(a))
welch <- function(a, b, alpha = 0.05) {           # treated a minus control b on the log scale: difference, 95% t-interval (Welch-Satterthwaite df), two-sided p
  va <- var(a) / length(a); vb <- var(b) / length(b); se <- sqrt(va + vb); df <- (va + vb)^2 / (va^2 / (length(a) - 1) + vb^2 / (length(b) - 1)); d <- mean(a) - mean(b); t <- qt(1 - alpha / 2, df)
  c(d, d - t * se, d + t * se, df, 2 * pt(abs(d / se), df, lower.tail = FALSE))
}
pooled3 <- function(groups, i, j, alpha = 0.05) {   # difference of means groups[j] - groups[i] with the within-group variance pooled over all groups (df = N - k)
  n <- vapply(groups, length, 1); ss <- sum(vapply(groups, function(g) sum((g - mean(g))^2), 1)); df <- sum(n) - length(groups)
  se <- sqrt(ss / df * (1 / n[i] + 1 / n[j])); d <- mean(groups[[j]]) - mean(groups[[i]]); t <- qt(1 - alpha / 2, df)
  c(d, d - t * se, d + t * se, df, 2 * pt(abs(d / se), df, lower.tail = FALSE))
}
pairwise_p_anova <- function(groups, i, j) pooled3(groups, i, j)[5]
## Wilcoxon rank-sum, two-sided; exact when there are no ties (sizes here are 4 and 4), normal approximation with continuity correction otherwise (scipy method='auto' = R wilcox.test)
pairwise_p_wilcox <- function(a, b) suppressWarnings(wilcox.test(a, b, alternative = "two.sided", correct = TRUE)$p.value)
printed_match <- function(val, s) {                 # TRUE if val agrees with the printed string s to the printed precision
  s <- trimws(as.character(s)); if (startsWith(s, "<")) return(val < as.numeric(substring(s, 2))); f <- as.numeric(s)
  if (grepl("e", tolower(s), fixed = TRUE)) { pr <- strsplit(tolower(s), "e", fixed = TRUE)[[1]]; nd <- if (grepl(".", pr[1], fixed = TRUE)) nchar(strsplit(pr[1], ".", fixed = TRUE)[[1]][2]) else 0; return(abs(val - f) <= 0.5 * 10^(as.integer(pr[2]) - nd) + 1e-300) }
  nd <- if (grepl(".", s, fixed = TRUE)) nchar(strsplit(s, ".", fixed = TRUE)[[1]][2]) else 0; abs(val - f) <= 0.5 * 10^(-nd) + 1e-12
}
bh_ <- function(p) p.adjust(p, "BH")
fmt <- function(x, d = 2) if (is.null(x) || is.na(x)) "" else sprintf(paste0("%.", d, "f"), x)
fmtp <- function(p) if (is.null(p) || is.na(p)) "" else if (p >= 0.001) sprintf("%.3f", p) else sprintf("%.2e", p)
