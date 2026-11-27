## check_readme_numbers.R - R port of wp4_islets/check_readme_numbers.py: every headline number of the README of the islet cross-check is rebuilt from the result files and searched for, in the format used in the README (so a typing
## slip is caught). Prints one line per check and ends with the number of failures (exit status 1 if any). No analysis: it only reads the results of the R run (or, with CHECK_SOURCE=archive, the archived Python results) and README.md.
## The README itself is the archived text (wp4_islets/README.md): a number of the R run that does not appear in it as written is reported as FAIL, which is the agreement test of the port.
## Usage: Rscript check_readme_numbers.R          CHECK_SOURCE = R (default; results_R/wp4/...) | archive (wp4_islets/...)
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(R_PORT, "pycompat.R")); suppressMessages(library(jsonlite))
ARCH <- normalizePath(file.path(R_PORT, "..", "wp4_islets")); W4 <- file.path(R_PORT, "results_R", "wp4"); SRC <- Sys.getenv("CHECK_SOURCE", "R")
at <- function(rdir, adir, file) if (SRC == "archive") file.path(ARCH, adir, file) else file.path(W4, rdir, file)      # R folder of the run | archive folder
read_text <- function(p) { b <- readBin(p, "raw", file.size(p)); s <- rawToChar(b); Encoding(s) <- "UTF-8"; s }
RD <- read_text(file.path(ARCH, "README.md")); fails <- character(0)
has <- function(label, s) { ok <- grepl(enc2utf8(s), RD, fixed = TRUE, useBytes = TRUE); cat(if (ok) "ok   " else "FAIL ", label, "   [", substr(s, 1, 90), "]\n", sep = ""); if (!ok) fails <<- c(fails, label) }
sg <- function(x, d = 3) { s <- py_fmt_f(x, d); ifelse(startsWith(s, "-"), s, paste0("+", s)) }      # '%+.*f'
f3 <- function(x) py_fmt_f(x, 3L); f2 <- function(x) py_fmt_f(x, 2L)
rcsv <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
if (SRC == "archive") I <- rcsv(file.path(ARCH, "results", "index.csv")) else { fs <- list.files(file.path(W4, "results"), "^index_.*[.]csv$", full.names = TRUE); I <- do.call(rbind, lapply(fs, rcsv)) }      # the per-dataset indexes of the R run
ms <- fromJSON(at("meta", "meta", "meta_summary.json"), simplifyVector = FALSE)
NAME <- c(gse81608 = "GSE81608", gse84133_human = "GSE84133 human", gse86469 = "GSE86469", gse81547 = "GSE81547", emtab5061 = "E-MTAB-5061", gse84133_mouse = "GSE84133 mouse")
for (t in names(NAME)) {
  n <- NAME[[t]]; u <- I[I$tag == t & I$variant == "primary" & I$mode == "unadjusted", ]; a <- I[I$tag == t & I$variant == "primary" & I$mode == "adjusted", ]
  if (nrow(u) == 0 || nrow(a) == 0) { cat("FAIL ", n, ": no row in the index\n", sep = ""); fails <- c(fails, paste(n, "index")); next }
  has(paste0(n, ": section 4 row"), sprintf("| %s | %s; %s | %s; %s |", n, sg(u$median_rho_consensus), sg(u$E), sg(a$median_rho_consensus), sg(a$E)))
  if (t != "gse84133_mouse") has(paste0(n, ": section 5 row"), sprintf("| %s | %d | %d | %s | %s |", n, as.integer(a$beta_cells), as.integer(a$consensus_pairs_evaluable), sg(a$E), f3(a$SE_E)))
  has(paste0(n, ": section 3/1 evaluable pairs count in Table-1 list"), as.character(as.integer(a$consensus_pairs_evaluable)))
}
V <- list(); for (v in ms$variants) V[[paste(v$variant, v$mode)]] <- v
est <- function(key, pre) Filter(function(x) startsWith(x$estimator, pre), V[[key]]$estimators)[[1]]; ci <- function(e) sprintf("%s [%s, %s]", sg(e$mu), sg(e$lo), sg(e$hi))
e <- est("primary adjusted", "REML+HK (primary)"); has("pooled primary", ci(e)); has("tau2 and I2", sprintf("tau2 %s, I2 %s%%", py_fmt_f(e$tau2, 4L), py_fmt_f(100 * e$I2, 0L)))
for (k in list(c("nondiabetic adjusted", "non-diabetic adjusted"), c("primary unadjusted", "unadjusted"), c("nondiabetic unadjusted", "non-diabetic unadjusted"))) has(k[2], ci(est(k[1], "REML+HK (primary)")))
for (y in V[["primary adjusted"]]$estimators) if (y$estimator %in% c("DL+HK", "DL, normal CI", "fixed effect")) has(paste("sensitivity", y$estimator), sprintf("[%s, %s]", sg(y$lo), sg(y$hi)))
for (y in V[["primary adjusted"]]$leave_one_out) has(paste("leave out", y$left_out), ci(y))
P <- rcsv(at("meta", "meta", "pair_level_pooled_rho.csv"))
for (i in seq_len(nrow(P))) { g <- c(P$gene1[i], P$gene2[i])
  if (identical(g, c("Tcf7l2", "Itpr1"))) has(paste0("pair ", g[1], "-", g[2]), f3(P$pooled_rho[i]))
  else if (any(vapply(list(c("Cacna1c", "Cacna1d"), c("Cacna1d", "Ptbp1"), c("Cacna1d", "Glp1r")), function(z) identical(g, z), TRUE))) has(paste0("pair ", g[1], "-", g[2]), sprintf("%s [%s, %s]", sg(P$pooled_rho[i]), sg(P$lo[i]), sg(P$hi[i]))) }
gg <- rcsv(at("meta", "meta", "gene_graph_meta.csv")); g0 <- gg[gg$estimator == "REML+HK", ][1, ]; has("gene-graph pooled", sprintf("%s [%s, %s]", sg(g0$mu), sg(g0$lo), sg(g0$hi)))
G <- rcsv(at("meta", "meta", "gene_graph_by_dataset.csv")); for (i in seq_len(nrow(G))) has(paste("gene-graph SE ratio", G$tag[i]), f2(G$SE_ratio[i]))
Ms <- rcsv(at("metacell", "metacell", "metacell_summary.csv"))
for (i in seq_len(nrow(Ms))) { r <- Ms[i, ]; if (is.na(r$consensus_pairs_evaluable_med)) next
  has(sprintf("metacell %s size %d", r$tag, r$size), if (r$size == 20) sprintf("%d; %d (%d); %s; %s", as.integer(r$metacells_median), as.integer(r$consensus_pairs_evaluable_med), as.integer(r$cell_level_pairs), sg(r$median_rho_med), sg(r$E_med))
                                                      else sprintf("%d; %d; %s; %s", as.integer(r$metacells_median), as.integer(r$consensus_pairs_evaluable_med), sg(r$median_rho_med), sg(r$E_med))) }
pc <- fromJSON(at("perturbation", "perturbation", "gse60803_concordance.json"), simplifyVector = FALSE)
has("perturbation primary", sprintf("| 17 | 7 | %s | %s | %s |", f2(pc$all_donors$binomial_p), sg(pc$all_donors$spearman_rho), f2(pc$all_donors$spearman_p)))
has("perturbation S114", sprintf("| 17 | 10 | %s | %s | %s |", f2(pc$S114_only$binomial_p), sg(pc$S114_only$spearman_rho), f2(pc$S114_only$spearman_p)))
Cd <- rcsv(at("perturbation_decisions", "perturbation", "candidates_decisions.csv")); has("candidates", sprintf("%d unique series", nrow(Cd))); has("dropped", sprintf("Dropped %d", sum(Cd$decision == "drop")))
cat("FAILURES: ", length(fails), "\n", sep = ""); if (sys.nframe() == 0L) quit(status = if (length(fails)) 1L else 0L)
