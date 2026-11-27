## wider_search_analyze.R - R port of wp4_islets/wider_search_analyze.py: B4 and B5 of PRESPEC_WP4b_positive_control_and_wider_search.md for the kept series of wider_search/frozen_candidates.csv (fixed and fingerprinted before any value of the
## new series was read). Per analysable series: log2 fold change of the 17 genes (perturbed arm minus control arm as frozen; arrays on the log2 scale of the series matrix; RNA-seq log2(CPM + 1) of NCBI-generated counts), multiplied by the
## expected sign (-1 for MIRROR designs), compared with the INS-1 log2 fold change of the matching arm (wp5_reporting/effect_sizes.csv, RT-qPCR, 'as reported'; GSE60803 from perturbation/gse60803_fold_changes.csv): genes with both nonzero,
## exact binomial test of the equal signs, Spearman rho. Series whose data cannot be obtained by plain download (K4) are listed as not analysable. Pooling (B5): REML + Hartung-Knapp of the Fisher z of rho (variance 1.06 / (n - 3));
## outcome words: CI entirely above 0 "concordant", entirely below 0 "discordant", else "inconclusive"; fewer than 3 series: not applied.
## Input: WS_FROZEN (default the archived wider_search folder: frozen_candidates.csv, sample_records_<acc>.tsv), the raw data in WP4_DATA (GSE99503, GSE167073; wp4_http.R), WP4_PERT (folder of gse60803_fold_changes.csv, default results_R/wp4/perturbation).
## Output (WP4_WS_OUT, default results_R/wp4/wider_search): fold_changes_<acc>.csv, concordance_results.csv, concordance_meta.csv, concordance_report.md.   Usage: Rscript wider_search_analyze.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_stats.R")); source(file.path(HERE, "wp4_meta_lib.R")); options(py_json_nan = TRUE)
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); OUT <- Sys.getenv("WP4_WS_OUT", file.path(R_PORT, "results_R", "wp4", "wider_search")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
WS_FROZEN <- Sys.getenv("WS_FROZEN", file.path(WP4_DIR, "wider_search")); PERT <- Sys.getenv("WP4_PERT", file.path(R_PORT, "results_R", "wp4", "perturbation")); ES_CSV <- Sys.getenv("ES_CSV", normalizePath(file.path(WP4_DIR, "..", "wp5_reporting", "effect_sizes.csv")))
MOUSE <- c(Hes1 = "Hes1", Wnt5a = "Wnt5a", Wnt5b = "Wnt5b", Wnt9a = "Wnt9a", Lef1 = "Lef1", Tcf7 = "Tcf7", Tcf7l2 = "Tcf7l2", Itpr1 = "Itpr1", Ryr2 = "Ryr2", Ryr3 = "Ryr3", Cacna1c = "Cacna1c", Cacna1d = "Cacna1d", Kcnj11 = "Kcnj11", Ptbp1 = "Ptbp1", Glut2 = "Slc2a2", Ins1 = "Ins1", Glp1r = "Glp1r")
ins1 <- function() {
  E <- read.csv(ES_CSV, stringsAsFactors = FALSE, check.names = FALSE); E <- E[E$readout == "RT-qPCR" & E$variant == "as reported", ]; out <- list()
  for (a in list(c("DAPT", "DAPT vs CTRL"), c("DKK-1", "DKK-1 vs CTRL"))) { s <- E[E$contrast == a[2], ]; out[[a[1]]] <- list(fc = setNames(s$log2fc_geo[match(GENES, s$feature)], GENES), q = setNames(s$bh_q_17genes[match(GENES, s$feature)], GENES)) }
  out
}
titles <- function(acc) { t <- read.delim(file.path(WS_FROZEN, sprintf("sample_records_%s.tsv", acc)), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", colClasses = "character"); setNames(t$title, t$gsm) }
fc_gse60803 <- function() { d <- read.csv(file.path(PERT, "gse60803_fold_changes.csv"), stringsAsFactors = FALSE); setNames(d$log2fc_all_donors[match(GENES, d$gene)], GENES) }
fc_gse99503 <- function() {
  X <- read.delim(file.path(DATA, "GSE99503", "GSE99503_counts_GRCh38.p13_NCBI.tsv.gz"), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1); T <- titles("GSE99503"); ctrl <- names(T)[grepl("^Neg siRNA #[0-9]$", T)]; pert <- names(T)[grepl("^Neg siRNA \\+ 1-AKP #[0-9]$", T)]
  stopifnot(length(ctrl) == 3, length(pert) == 3, all(c(ctrl, pert) %in% colnames(X))); cs <- as.matrix(X[, c(ctrl, pert)]); L <- log2(sweep(cs, 2, colSums(cs), "/") * 1e6 + 1); fc <- rowMeans(L[, pert]) - rowMeans(L[, ctrl])
  list(fc = setNames(vapply(GENES, function(g) if (as.character(ENTREZ[[g]]) %in% names(fc)) unname(fc[as.character(ENTREZ[[g]])]) else NA_real_, 1), GENES), info = list(genes_found = sum(as.character(ENTREZ[GENES]) %in% names(fc)), ctrl = length(ctrl), pert = length(pert)))
}
fc_gse167073 <- function() {
  d <- file.path(DATA, "GSE167073"); T <- titles("GSE167073"); ctrl <- names(T)[grepl("Fzd9 wt", T, fixed = TRUE)]; pert <- names(T)[grepl("Fzd9 ko", T, fixed = TRUE)]; stopifnot(length(ctrl) == 4, length(pert) == 4); tgt <- names(MOUSE)[match(MOUSE, MOUSE)]
  con <- file(file.path(d, "GPL17400_data.txt"), "r", encoding = "UTF-8"); pm <- list(); ig <- NULL; repeat {
    ch <- readLines(con, n = 20000L, warn = FALSE); if (!length(ch)) break
    for (l in ch) {
      if (is.null(ig)) { if (startsWith(l, "ID\t")) { hdr <- strsplit(l, "\t", fixed = TRUE)[[1]]; ig <- match("gene_assignment", hdr) }; next }
      p <- strsplit(l, "\t", fixed = TRUE)[[1]]; if (length(p) < ig || p[ig] %in% c("---", "")) next
      for (a in strsplit(p[ig], " /// ", fixed = TRUE)[[1]]) { q <- strsplit(a, " // ", fixed = TRUE)[[1]]; if (length(q) > 1 && trimws(q[2]) %in% MOUSE) pm[[p[1]]] <- union(pm[[p[1]]], trimws(q[2])) }
    }
  }
  close(con); sm <- read_series_matrix(file.path(d, "GSE167073_series_matrix.txt.gz")); X <- as.matrix(sm$table[, -1]); rownames(X) <- sm$table[[1]]; storage.mode(X) <- "double"; cols <- colnames(X)
  ci <- match(ctrl, cols); pi_ <- match(pert, cols); out <- setNames(rep(NA_real_, length(GENES)), GENES); probes <- list()
  for (g in GENES) {
    cand <- names(pm)[vapply(pm, function(s) MOUSE[[g]] %in% s, TRUE)]; cand <- cand[cand %in% rownames(X)]; if (!length(cand)) next
    best <- cand[which.max(vapply(cand, function(x) mean(X[x, ci]), 1))]; v <- X[best, ]; stopifnot("matrix not on a log2 scale" = max(v) < 30); out[[g]] <- mean(v[pi_]) - mean(v[ci]); probes[[g]] <- list(best, length(cand))
  }
  list(fc = out, info = list(genes_found = sum(!is.na(out)), probes = probes))
}
concord <- function(fc, sign, vec) {
  x <- fc * sign; ok <- !is.na(x) & !is.na(vec) & x != 0 & vec != 0; n <- sum(ok); if (n < 3) return(list(n_genes = n))
  xs <- x[ok]; vs <- vec[ok]; k <- sum(sign(xs) == sign(vs)); sp <- py_spearmanr(xs, vs)
  list(n_genes = n, n_same_sign = k, binomial_p = py_binom_p(k, n), spearman_rho = sp$rho, spearman_p = sp$p, same_sign_genes = names(xs)[sign(xs) == sign(vs)])
}
fs <- function(x, nd) paste0(if (x < 0) "-" else "+", py_fmt_f(abs(x), nd))
run <- function() {
  F <- read.csv(file.path(WS_FROZEN, "frozen_candidates.csv"), stringsAsFactors = FALSE, check.names = FALSE, na.strings = ""); K <- F[F$decision == "keep", ]; V <- ins1(); res <- list()
  for (i in seq_len(nrow(K))) {
    k <- K[i, ]; acc <- k$accession; base <- list(accession = acc, arm = k$arm, expected_sign = as.integer(k$expected_sign), mirror = as.integer(k$mirror), weak_match_note = k$weak_match_note, title = F$title[F$accession == acc][1], k4 = k$k4)
    if (startsWith(as.character(k$k4), "NOT")) { res[[length(res) + 1]] <- c(base, list(status = paste0("not analysable (K4): ", k$data_route))); cat(acc, "not analysable\n"); next }
    got <- tryCatch({ if (acc == "GSE60803") list(fc = fc_gse60803(), info = list()) else if (acc == "GSE99503") fc_gse99503() else if (acc == "GSE167073") fc_gse167073() else stop("no analysis for ", acc) }, error = function(e) e)
    if (inherits(got, "error")) { res[[length(res) + 1]] <- c(base, list(status = sprintf("not analysable (data step failed: %s)", substr(conditionMessage(got), 1, 150)))); cat(acc, "FAILED", conditionMessage(got), "\n"); next }
    fc <- got$fc; info <- got$info; sg <- as.integer(k$expected_sign); vec <- V[[k$arm]]$fc
    py_write_csv_repr(setNames(data.frame(GENES, unname(fc), sg, unname(fc * sg), unname(vec), stringsAsFactors = FALSE), c("gene", "log2fc_perturbed_minus_control", "expected_sign", "log2fc_signed", paste0("ins1_log2fc_", k$arm))), file.path(OUT, sprintf("fold_changes_%s.csv", acc)), eol = "\r\n")
    r <- concord(fc, sg, vec); q <- V[[k$arm]]$q; strong <- names(q)[!is.na(q) & q < 0.05]; sgn <- (fc * sg)[strong]; sgn <- sgn[!is.na(sgn) & sgn != 0]
    r <- c(r, base); r$status <- "analysed"; r$genes_found <- if (!is.null(info$genes_found)) info$genes_found else sum(!is.na(fc)); r$strong_ins1_genes <- length(sgn); r$strong_same_sign <- sum(sign(sgn) == sign(vec[names(sgn)])); res[[length(res) + 1]] <- r
    cat(sprintf("%-10s %-6s n_genes %2d same sign %2d P %.3f rho %+.3f (P %.3f) | strong INS-1 genes same sign %d of %d\n", acc, k$arm, if (is.null(r$n_genes)) 0L else r$n_genes, if (is.null(r$n_same_sign)) -1L else r$n_same_sign, if (is.null(r$binomial_p)) NaN else r$binomial_p,
                if (is.null(r$spearman_rho)) NaN else r$spearman_rho, if (is.null(r$spearman_p)) NaN else r$spearman_p, r$strong_same_sign, r$strong_ins1_genes))
  }
  cols <- c("accession", "arm", "expected_sign", "mirror", "weak_match_note", "title", "k4", "status", "n_genes", "n_same_sign", "binomial_p", "spearman_rho", "spearman_p", "same_sign_genes", "genes_found", "strong_ins1_genes", "strong_same_sign")
  R <- do.call(rbind, lapply(res, function(r) { x <- lapply(cols, function(k) { v <- r[[k]]; if (is.null(v)) NA else if (k == "same_sign_genes") py_list_repr(v) else v }); names(x) <- cols; as.data.frame(x, stringsAsFactors = FALSE) }))
  py_write_csv_repr(R, file.path(OUT, "concordance_results.csv"), eol = "\r\n"); A <- R[R$status == "analysed", ]; A <- A[A$n_genes >= 12, ]; metas <- list()
  meta <- function(label, S) {
    if (nrow(S) < 3) return(list(set = label, k = nrow(S), applied = FALSE, pooled_rho = NA_real_, lo = NA_real_, hi = NA_real_, I2 = NA_real_, tau2 = NA_real_, outcome = "rule not applied (fewer than 3 analysable series)"))
    z <- atanh(pmin(pmax(S$spearman_rho, -0.999), 0.999)); se <- sqrt(1.06 / (S$n_genes - 3)); r <- pool(z, se); lo <- tanh(r$lo); hi <- tanh(r$hi)
    list(set = label, k = nrow(S), applied = TRUE, pooled_rho = tanh(r$mu), lo = lo, hi = hi, I2 = r$I2, tau2 = r$tau2, outcome = if (lo > 0) "concordant" else if (hi < 0) "discordant" else "inconclusive")
  }
  metas[[1]] <- meta("PRIMARY: all analysable series", A); metas[[2]] <- meta("Wnt series (DKK-1 arm)", A[A$arm == "DKK-1", ]); metas[[3]] <- meta("Notch series (DAPT arm)", A[A$arm == "DAPT", ])
  metas[[4]] <- meta("mirror designs excluded", A[A$mirror == 0, ]); metas[[5]] <- meta("GSE60803 excluded", A[A$accession != "GSE60803", ])
  M <- do.call(rbind, lapply(metas, as.data.frame, stringsAsFactors = FALSE)); py_write_csv_repr(M, file.path(OUT, "concordance_meta.csv"), eol = "\r\n")
  L <- c("# Wider search: concordance of kept perturbation series with the INS-1 arms (rules fixed in PRESPEC_WP4b... B4/B5 before any value was read)", "", "| series | arm | expected sign | status | genes | same sign | exact binomial P | Spearman rho (P) | strong INS-1 genes same sign |", "|---|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(R))) { r <- R[i, ]; sgn <- sprintf("%+d%s", r$expected_sign, if (r$mirror == 1) " (mirror)" else "")
    L <- c(L, if (r$status == "analysed") sprintf("| %s | %s | %s | analysed | %d | %d | %s | %s (%s) | %d of %d |", r$accession, r$arm, sgn, as.integer(r$n_genes), as.integer(r$n_same_sign), py_fmt_f(r$binomial_p, 3), fs(r$spearman_rho, 3), py_fmt_f(r$spearman_p, 3), as.integer(r$strong_same_sign), as.integer(r$strong_ins1_genes))
                   else sprintf("| %s | %s | %s | %s | - | - | - | - | - |", r$accession, r$arm, sgn, substr(r$status, 1, 110))) }
  L <- c(L, "", "## Pooled (REML + Hartung-Knapp on Fisher z)", "", "| set | k | pooled rho | 95% CI | I2 | outcome |", "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(M))) { m <- M[i, ]; L <- c(L, sprintf("| %s | %d | %s | %s | %s | %s |", m$set, as.integer(m$k), if (m$applied) fs(m$pooled_rho, 3) else "-", if (m$applied) sprintf("[%s, %s]", fs(m$lo, 3), fs(m$hi, 3)) else "-", if (m$applied) sprintf("%.0f%%", 100 * m$I2) else "-", m$outcome)) }
  write_text_bytes(paste0(paste(L, collapse = .NL), .NL), file.path(OUT, "concordance_report.md")); cat(paste(L, collapse = "\n"), "\n")
}
if (sys.nframe() == 0L) run()
