## gene_graph_null.R - R port of wp4_islets/gene_graph_null.py: ADD-ON, POST HOC (written after the primary meta-analysis result was known; not in the protocol): does the verdict of 4.4 survive when the null keeps the gene graph of
## the consensus pairs? The pre-specified null draws the two genes of every consensus pair independently, but the evaluable consensus pairs share genes, so the mean of their correlations varies more than a mean over independent random
## pairs and the pre-specified SE of E is probably too small. For each dataset (depth-adjusted analysis, all frozen beta-cells, same eligible genes, same detection deciles, same partial-Spearman ranks as wp4_core analyse) every
## distinct gene of the evaluable consensus pairs is replaced by ONE random gene of its own detection decile (not one of the 17, no gene used twice in a draw) and the mean Fisher z of the same pairs is computed on the replaced genes;
## 2000 draws, seed 2028. E_graph = mean z of the observed pairs minus the mean of the null; SE_graph = SD of the null means. The meta-analysis (same estimators and margin) is repeated with E_graph and SE_graph.
## Output (WP4_OUT, default results_R/wp4/meta): gene_graph_by_dataset.csv, gene_graph_meta.csv, gene_graph_report.md. Needs the per-dataset summaries of run_datasets.R (WP4_RES; default the archived results).
## Usage: Rscript gene_graph_null.R [tag ...]
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_parse.R")); source(file.path(HERE, "wp4_meta_lib.R")); options(py_json_nan = TRUE)
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); RES <- Sys.getenv("WP4_RES", normalizePath(file.path(R_PORT, "..", "wp4_islets", "results"))); OUT <- Sys.getenv("WP4_OUT", file.path(R_PORT, "results_R", "wp4", "meta")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SEED <- 2028L; TAGS <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse"); HUMAN <- TAGS[1:5]
graph_null <- function(cells, thr = 0.10, nn = 2000L, seed = SEED) {
  P <- prep_adjusted(cells, thr); sym <- P$sym; pairs <- P$pairs; Rk <- P$Rk; idx <- P$idx; genes <- GENES[vapply(GENES, function(g) any(vapply(pairs, function(p) g %in% p, TRUE)), TRUE)]; gi <- setNames(seq_along(genes), genes)
  rho <- vapply(pairs, function(p) sum(Rk[idx[[sym(p[1])]], ] * Rk[idx[[sym(p[2])]], ]), 1); z_obs <- mean(atanh(pmin(pmax(rho, -CLIP), CLIP)))
  rng <- np_rng(seed); M <- matrix(0L, nn, length(genes))
  for (s in seq_len(nn)) { used <- integer(0); for (k in seq_along(genes)) { arr <- P$pool[[P$bins[idx[[sym(genes[k])]]] + 1L]]; repeat { x <- arr[np_integers(rng, 0, length(arr), 1) + 1L]; if (!(x %in% used)) break }; used <- c(used, x); M[s, k] <- x } }
  Z <- numeric(nn)
  for (p in pairs) Z <- Z + atanh(pmin(pmax(unlist(lapply(seq(1, nn, by = 250L), function(s) { e <- min(s + 249L, nn); rowSums(Rk[M[s:e, gi[[p[1]]]], , drop = FALSE] * Rk[M[s:e, gi[[p[2]]]], , drop = FALSE]) })), -CLIP), CLIP))
  Z <- Z / length(pairs); share <- vapply(genes, function(g) sum(vapply(pairs, function(p) g %in% p, TRUE)), 1)
  list(pairs = length(pairs), genes_in_pairs = length(genes), max_pairs_per_gene = as.integer(max(share)), z_obs = z_obs, z_null_graph = mean(Z), E_graph = z_obs - mean(Z), SE_graph = sd(Z))
}
fs <- function(x, nd) paste0(if (x < 0) "-" else "+", py_fmt_f(abs(x), nd))
a <- commandArgs(trailingOnly = TRUE); tags <- if (length(a)) a else TAGS
if (sys.nframe() == 0L) {
  rows <- list()
  for (tag in tags) {
    t0 <- Sys.time(); cells <- load_cells(tag, DATA); r <- graph_null(cells); s <- fromJSON(file.path(RES, sprintf("%s__primary__summary.json", tag)), simplifyVector = FALSE)$adjusted
    r$tag <- tag; r$E_prespecified <- if (is.null(s$E)) NA_real_ else s$E; r$SE_prespecified <- if (is.null(s$SE_E)) NA_real_ else s$SE_E; r$SE_ratio <- if (!is.na(r$SE_prespecified) && r$SE_prespecified != 0) r$SE_graph / r$SE_prespecified else NA_real_; rows[[length(rows) + 1]] <- r
    cat(sprintf("%-15s pairs %2d genes %2d (max %d pairs per gene) | E %+.3f -> %+.3f | SE %.3f -> %.3f (x%.2f) | %.0f s\n", tag, r$pairs, r$genes_in_pairs, r$max_pairs_per_gene, r$E_prespecified, r$E_graph, r$SE_prespecified, r$SE_graph, r$SE_ratio, as.numeric(difftime(Sys.time(), t0, units = "secs")))); rm(cells); gc(FALSE)
  }
  D <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))); D <- D[, c("pairs", "genes_in_pairs", "max_pairs_per_gene", "z_obs", "z_null_graph", "E_graph", "SE_graph", "tag", "E_prespecified", "SE_prespecified", "SE_ratio")]
  py_write_csv_repr(D, file.path(OUT, "gene_graph_by_dataset.csv"), eol = "\r\n"); H <- D[D$tag %in% HUMAN & D$pairs >= 3, ]; res <- list()
  L <- c("# WP4 add-on (post hoc): meta-analysis with a null that keeps the gene graph of the consensus pairs", "", "Same estimators and margin as meta_analysis.py; E_graph and SE_graph replace E and SE of the pre-specified analysis. The rule of 4.4 stays defined on the pre-specified SE; this is a sensitivity analysis.", "")
  for (spec in list(list("REML+HK", list(tau = "reml", ci = "hk", floor = FALSE)), list("REML+HK, factor floored at 1", list(tau = "reml", ci = "hk", floor = TRUE)), list("DL+HK", list(tau = "dl", ci = "hk")), list("DL, normal CI", list(tau = "dl", ci = "normal")), list("fixed effect", list(ci = "fixed")))) {
    r <- do.call(pool, c(list(H$E_graph, H$SE_graph), spec[[2]])); r$estimator <- spec[[1]]; r$verdict <- verdict(r$lo, r$hi); r$analysis <- "gene-graph null"; res[[length(res) + 1]] <- r
    L <- c(L, sprintf("* %-30s pooled E %s  95%% CI [%s, %s]  tau2 %s  I2 %.0f%%  -> %s", spec[[1]], fs(r$mu, 3), fs(r$lo, 3), fs(r$hi, 3), py_fmt_f(r$tau2, 4), 100 * r$I2, r$verdict))
  }
  L <- c(L, ""); for (t in H$tag) { R <- H[H$tag != t, ]; r <- pool(R$E_graph, R$SE_graph); L <- c(L, sprintf("  leave out %-15s pooled E %s  [%s, %s] -> %s", t, fs(r$mu, 3), fs(r$lo, 3), fs(r$hi, 3), verdict(r$lo, r$hi))) }
  Mt <- do.call(rbind, lapply(res, function(r) as.data.frame(r[c("mu", "se", "lo", "hi", "tau2", "I2", "Q", "k", "estimator", "verdict", "analysis")], stringsAsFactors = FALSE))); py_write_csv_repr(Mt, file.path(OUT, "gene_graph_meta.csv"), eol = "\r\n")
  write_text_bytes(paste(L, collapse = .NL), file.path(OUT, "gene_graph_report.md")); cat(paste(L, collapse = "\n"), "\n")
}
