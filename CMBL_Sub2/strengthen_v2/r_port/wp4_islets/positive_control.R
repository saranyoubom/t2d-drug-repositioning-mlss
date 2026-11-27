## positive_control.R - R port of wp4_islets/positive_control.py: WP4 extension A1-A3 of PRESPEC_WP4b_positive_control_and_wider_search.md (fixed before the Python script ran): what could the islet cross-check have detected?
## A1  synthetic spike-in: one latent factor per replicate is added (donor-centred log scale, before the depth adjustment) to the distinct genes of the evaluable consensus pairs of a dataset with loading lambda_g = sqrt(rho / (1 - rho)) x sd_g,
##     nominal rho in {0.05, 0.1, 0.2, 0.3, 0.5}; mode 'detected' (primary: only in cells where the gene is detected, zeros stay zeros) and 'all' (upper bound); the unchanged pipeline (depth-adjusted partial Spearman, excess E over the 2000
##     matched random pair sets of the primary run) is applied; 200 replicates; the five human datasets are pooled per replicate with the primary estimator (REML + Hartung-Knapp, margin +/-0.10). Gate: without injection the code must
##     reproduce E and SE of results/<tag>__primary__summary.json (adjusted) to 1e-9.
## A2  real-gene control: 100 random sets of 11 pairs of cytosolic ribosomal protein genes; E_RP per set against matched random pairs (deciles of detection; pool without any ribosomal gene and without the 17 genes; 2000 draws per decile combination).
## A3  margin sensitivity of the primary pooled interval (descriptive).
## Output (WP4_PC_OUT, default results_R/wp4/positive_control): inj_<tag>.rds (the R counterpart of inj_<tag>.npz), rp_<tag>.json, positive_control_by_dataset.csv, positive_control_pooled.csv, rp_control.csv, margin_sensitivity.csv, positive_control_report.md.
## Usage: Rscript positive_control.R [tag ...]   run the datasets named (default all six) and then the report;   Rscript positive_control.R report   only pool and write the report from the saved files. PC_TEST=1: 5 replicates (code test only).
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); R_PORT <- normalizePath(file.path(HERE, ".."))
source(file.path(HERE, "wp4_core.R")); source(file.path(HERE, "wp4_parse.R")); source(file.path(HERE, "wp4_meta_lib.R")); options(py_json_nan = TRUE)
DATA <- Sys.getenv("WP4_DATA", file.path(R_PORT, "results_R", "wp4_data")); RES <- Sys.getenv("WP4_RES", normalizePath(file.path(R_PORT, "..", "wp4_islets", "results")))
META <- Sys.getenv("WP4_META", normalizePath(file.path(R_PORT, "..", "wp4_islets", "meta"))); OUT <- Sys.getenv("WP4_PC_OUT", file.path(R_PORT, "results_R", "wp4", "positive_control"))
TAGS <- c("gse81608", "gse84133_human", "gse86469", "gse81547", "emtab5061", "gse84133_mouse"); HUMAN <- TAGS[1:5]; RHOS <- c(0.05, 0.1, 0.2, 0.3, 0.5); MODES <- c("detected", "all"); NREP <- 200L; NSETS <- 100L; NPAIRS <- 11L
if (nzchar(Sys.getenv("PC_TEST"))) { NREP <- 5L; NSETS <- 5L; OUT <- paste0(OUT, "_test") }; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RP <- strsplit(paste("RPL3 RPL4 RPL5 RPL6 RPL7 RPL7A RPL8 RPL9 RPL10 RPL10A RPL11 RPL12 RPL13 RPL13A RPL14 RPL15 RPL17 RPL18 RPL18A RPL19 RPL21 RPL22 RPL23 RPL23A RPL24 RPL26 RPL27 RPL27A RPL28 RPL29 RPL30 RPL31 RPL32 RPL34 RPL35 RPL35A RPL36",
                    "RPL37 RPL37A RPL38 RPL39 RPLP0 RPLP1 RPLP2 RPS2 RPS3 RPS3A RPS4X RPS5 RPS6 RPS7 RPS8 RPS9 RPS10 RPS11 RPS12 RPS13 RPS14 RPS15 RPS15A RPS16 RPS17 RPS18 RPS19 RPS20 RPS21 RPS23 RPS24 RPS25 RPS26 RPS27 RPS27A RPS28 RPS29 RPSA FAU UBA52"), " ", fixed = TRUE)[[1]]
mouse_symbol <- function(s) paste0(substr(s, 1, 1), tolower(substring(s, 2)))
clipv <- function(x) pmin(pmax(x, -CLIP), CLIP)
fs <- function(x, nd) paste0(if (x < 0) "-" else "+", py_fmt_f(abs(x), nd))

prepare <- function(cells) {
  P <- prep_adjusted(cells); sym <- P$sym; idx <- P$idx; pairs <- P$pairs; Rk <- P$Rk; d <- draws_main(P$pool, P$bins, P$pidx, NN, SEED_NULL); V <- null_values(Rk, d$IA, d$IC); zm <- rowMeans(atanh(clipv(V)))
  genes <- GENES[vapply(GENES, function(g) any(vapply(pairs, function(p) g %in% p, TRUE)), TRUE)]; cols <- vapply(genes, function(g) idx[[sym(g)]], 1L); gi <- setNames(seq_along(genes), genes)
  G <- P$Xc[, cols, drop = FALSE]; M <- as.matrix(cells$B[, P$elig[cols], drop = FALSE] > 0); rho0 <- vapply(P$pidx, function(p) sum(Rk[p[1], ] * Rk[p[2], ]), 1); z_obs <- mean(atanh(clipv(rho0)))
  list(n = nrow(cells$B), donor = cells$donor, Z = P$Z, G = G, M = M, sd = apply(G, 2, np_std), pair_k = lapply(pairs, function(p) c(gi[[p[1]]], gi[[p[2]]])), zm = zm, zm_mean = mean(zm), SE = sd(zm), E_primary = z_obs - mean(zm),
       npairs = length(pairs), ngenes = length(genes), own = P$own, idx = idx, labels_e = P$labels_e, det_e = P$det_e, Rk = Rk)
}
one_rep <- function(P, rho, mode, rng) {
  z <- np_standard_normal(rng, P$n); z <- center_by_donor(matrix(z, ncol = 1), P$donor, "count")[, 1]; z <- z / np_std(z); add <- outer(z, sqrt(rho / (1 - rho)) * P$sd); if (mode == "detected") add <- add * P$M
  R <- adjusted_ranks(center_by_donor(P$G + add, P$donor, "count"), P$Z); zobs <- mean(atanh(clipv(vapply(P$pair_k, function(k) sum(R[k[1], ] * R[k[2], ]), 1))))
  c(zobs - P$zm_mean, (1 + sum(P$zm >= zobs)) / (NN + 1))
}
null_cell_mean <- function(Rk, pool, da, dc, seed, n = NN) {
  rng <- np_rng(seed); pa <- pool[[da + 1L]]; pc <- pool[[dc + 1L]]; ia <- pa[np_integers(rng, 0, length(pa), n) + 1L]; ic <- pc[np_integers(rng, 0, length(pc), n) + 1L]; bad <- ia == ic; k <- 0L
  while (any(bad) && k < 1000L) { ic[bad] <- pc[np_integers(rng, 0, length(pc), sum(bad)) + 1L]; bad <- ia == ic; k <- k + 1L }
  v <- unlist(lapply(seq(1, n, by = 500L), function(s) { e <- min(s + 499L, n); rowSums(Rk[ia[s:e], , drop = FALSE] * Rk[ic[s:e], , drop = FALSE]) })); mean(atanh(clipv(v)))
}
rp_control <- function(cells, P, di) {
  ids <- read.csv(file.path(WP4_DIR, "rp_ids_human.csv"), stringsAsFactors = FALSE, colClasses = "character")
  lab <- if (cells$id_type == "entrez") setNames(ids$entrez, ids$symbol)[RP] else if (cells$id_type == "ensembl") setNames(ids$ensembl, ids$symbol)[RP] else if (cells$species == "mouse") setNames(vapply(RP, mouse_symbol, ""), RP) else setNames(RP, RP)
  present <- RP[lab[RP] %in% cells$genes]; eligible <- unname(lab[present][lab[present] %in% names(P$idx)]); mp <- make_pool(P$det_e, P$labels_e, c(P$own, unname(lab[RP]))); bins <- mp$bins; pool <- mp$pool
  rng <- np_rng(2030L + di); table <- list(); Es <- numeric(0); rhos <- numeric(0)
  for (s in seq_len(NSETS)) {
    zs <- numeric(0); ns <- numeric(0)
    for (q in seq_len(NPAIRS)) {
      ac <- np_choice_noreplace_idx(rng, length(eligible), 2) + 1L; ia <- P$idx[[eligible[ac[1]]]]; ic <- P$idx[[eligible[ac[2]]]]; r <- sum(P$Rk[ia, ] * P$Rk[ic, ]); zs <- c(zs, atanh(clipv(r))); rhos <- c(rhos, r)
      key <- paste(bins[ia], bins[ic]); if (is.null(table[[key]])) table[[key]] <- null_cell_mean(P$Rk, pool, bins[ia], bins[ic], c(2030, di, bins[ia], bins[ic])); ns <- c(ns, table[[key]])
    }
    Es <- c(Es, mean(zs) - mean(ns))
  }
  list(rp_genes_in_data = length(present), rp_genes_eligible = length(eligible), E_RP_mean = mean(Es), E_RP_sd = sd(Es), share_E_RP_ge_0_10 = mean(Es >= 0.10), mean_rho_RP_pairs = mean(rhos))
}
run_dataset <- function(tag) {
  di <- match(tag, TAGS) - 1L; t0 <- Sys.time(); cells <- load_cells(tag, DATA); P <- prepare(cells); s <- fromJSON(file.path(RES, sprintf("%s__primary__summary.json", tag)), simplifyVector = FALSE)$adjusted
  if (!(abs(P$E_primary - s$E) < 1e-9 && abs(P$SE - s$SE_E) < 1e-9)) stop(sprintf("GATE FAILED %s: E %.12f vs %.12f, SE %.12f vs %.12f", tag, P$E_primary, s$E, P$SE, s$SE_E))
  cat(sprintf("%-15s gate ok (E %+.4f, SE %.4f; %d pairs, %d genes) | prepared in %.0f s\n", tag, P$E_primary, P$SE, P$npairs, P$ngenes, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  Erep <- array(0, c(length(RHOS), length(MODES), NREP)); Pemp <- array(0, c(length(RHOS), length(MODES), NREP))
  for (ri in seq_along(RHOS)) for (mi in seq_along(MODES)) for (rep in seq_len(NREP)) { o <- one_rep(P, RHOS[ri], MODES[mi], np_rng(c(2029, di, ri - 1, mi - 1, rep - 1))); Erep[ri, mi, rep] <- o[1]; Pemp[ri, mi, rep] <- o[2] }
  saveRDS(list(Erep = Erep, Pemp = Pemp, E_primary = P$E_primary, SE = P$SE, rhos = RHOS, nrep = NREP), file.path(OUT, sprintf("inj_%s.rds", tag)))
  cat(sprintf("%-15s injection done (%.0f s): recovered E at nominal rho 0.3: detected %+.3f, all %+.3f\n", tag, as.numeric(difftime(Sys.time(), t0, units = "secs")), mean(Erep[match(0.3, RHOS), 1, ]), mean(Erep[match(0.3, RHOS), 2, ])))
  rp <- rp_control(cells, P, di); write_json_typed(list(rp_genes_in_data = as.integer(rp$rp_genes_in_data), rp_genes_eligible = as.integer(rp$rp_genes_eligible), E_RP_mean = rp$E_RP_mean, E_RP_sd = rp$E_RP_sd, share_E_RP_ge_0_10 = rp$share_E_RP_ge_0_10, mean_rho_RP_pairs = rp$mean_rho_RP_pairs), file.path(OUT, sprintf("rp_%s.json", tag)))
  cat(sprintf("%-15s ribosomal control: %d genes eligible | E_RP mean %+.3f (SD %.3f), share >= 0.10: %.2f | %.0f s\n", tag, rp$rp_genes_eligible, rp$E_RP_mean, rp$E_RP_sd, rp$share_E_RP_ge_0_10, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
report <- function() {
  D <- list(); for (t in TAGS) if (file.exists(file.path(OUT, sprintf("inj_%s.rds", t)))) D[[t]] <- readRDS(file.path(OUT, sprintf("inj_%s.rds", t))); rows <- list()
  for (t in names(D)) for (mi in seq_along(MODES)) for (ri in seq_along(RHOS)) {
    d <- D[[t]]; E <- d$Erep[ri, mi, ]; zn <- atanh(RHOS[ri])
    rows[[length(rows) + 1]] <- list(tag = t, mode = MODES[mi], rho = RHOS[ri], nominal_z = zn, E_primary = d$E_primary, SE_E = d$SE, E_recovered_mean = mean(E), E_recovered_sd = sd(E), attenuation = (mean(E) - d$E_primary) / zn, share_empirical_P_lt_0_05 = mean(d$Pemp[ri, mi, ] < 0.05))
  }
  A <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)); py_write_csv_repr(A, file.path(OUT, "positive_control_by_dataset.csv"), eol = "\r\n"); Pp <- NULL
  if (all(HUMAN %in% names(D))) {
    se <- vapply(HUMAN, function(t) D[[t]]$SE, 1); prow <- list(); NO <- "no meaningful excess co-expression"; MORE <- "consensus pairs are co-expressed more than random pairs"
    for (mi in seq_along(MODES)) for (ri in seq_along(RHOS)) {
      V <- character(NREP); mu <- numeric(NREP); lo <- numeric(NREP); hi <- numeric(NREP)
      for (rep in seq_len(NREP)) { r <- pool(vapply(HUMAN, function(t) D[[t]]$Erep[ri, mi, rep], 1), se); V[rep] <- verdict(r$lo, r$hi); mu[rep] <- r$mu; lo[rep] <- r$lo; hi[rep] <- r$hi }
      prow[[length(prow) + 1]] <- list(mode = MODES[mi], rho = RHOS[ri], replicates = NREP, share_no_meaningful_excess = mean(V == NO), share_inconclusive = mean(V == "inconclusive"), share_co_expressed_more = mean(V == MORE),
                                       pooled_E_mean = mean(mu), pooled_E_sd = sd(mu), pooled_lo_mean = mean(lo), pooled_hi_mean = mean(hi))
    }
    Pp <- do.call(rbind, lapply(prow, as.data.frame, stringsAsFactors = FALSE)); py_write_csv_repr(Pp, file.path(OUT, "positive_control_pooled.csv"), eol = "\r\n")
  } else py_write_csv_repr(data.frame(), file.path(OUT, "positive_control_pooled.csv"), eol = "\r\n")
  rpl <- list(); for (t in TAGS) if (file.exists(file.path(OUT, sprintf("rp_%s.json", t)))) rpl[[length(rpl) + 1]] <- c(list(tag = t), fromJSON(file.path(OUT, sprintf("rp_%s.json", t)), simplifyVector = FALSE))
  R <- do.call(rbind, lapply(rpl, as.data.frame, stringsAsFactors = FALSE)); R$SE_E_consensus <- vapply(R$tag, function(t) if (t %in% names(D)) D[[t]]$SE else NA_real_, 1); R$E_RP_in_SE_units <- R$E_RP_mean / R$SE_E_consensus
  py_write_csv_repr(R, file.path(OUT, "rp_control.csv"), eol = "\r\n")
  ms <- fromJSON(file.path(META, "meta_summary.json"), simplifyVector = FALSE); vv <- Filter(function(v) v$variant == "primary" && v$mode == "adjusted", ms$variants)[[1]]; e <- Filter(function(x) startsWith(x$estimator, "REML+HK (primary)"), vv$estimators)[[1]]
  M <- data.frame(margin = c(0.05, 0.075, 0.10, 0.15, 0.20), verdict = vapply(c(0.05, 0.075, 0.10, 0.15, 0.20), function(m) verdict(e$lo, e$hi, m), ""), stringsAsFactors = FALSE); py_write_csv_repr(M, file.path(OUT, "margin_sensitivity.csv"), eol = "\r\n"); mmin <- max(abs(e$lo), abs(e$hi))
  L <- c("# WP4 extension A1-A3: positive control of the islet cross-check (reading rules fixed in PRESPEC_WP4b_positive_control_and_wider_search.md before this ran)", "")
  if (!is.null(Pp)) {
    L <- c(L, sprintf("## A1. Pooled verdict over the five human datasets, per nominal rho (200 replicates; no injection = the primary result, pooled E %s [%s, %s])", fs(e$mu, 3), fs(e$lo, 3), fs(e$hi, 3)), "",
           "| mode | nominal rho | no meaningful excess | inconclusive | co-expressed more | pooled E mean (SD) | pooled CI mean |", "|---|---|---|---|---|---|---|")
    for (i in seq_len(nrow(Pp))) { r <- Pp[i, ]; L <- c(L, sprintf("| %s | %s | %s | %s | %s | %s (%s) | [%s, %s] |", r$mode, py_fmt_f(r$rho, 2), py_fmt_f(r$share_no_meaningful_excess, 2), py_fmt_f(r$share_inconclusive, 2), py_fmt_f(r$share_co_expressed_more, 2), fs(r$pooled_E_mean, 3), py_fmt_f(r$pooled_E_sd, 3), fs(r$pooled_lo_mean, 3), fs(r$pooled_hi_mean, 3))) }
    d <- Pp[Pp$mode == "detected", ]; fr <- setNames(d$share_no_meaningful_excess, d$rho); r1 <- fr[["0.3"]] <= 0.20; ok <- as.numeric(names(fr))[fr <= 0.20]
    L <- c(L, "", sprintf("**Reading rule R1 (mode \"detected\", nominal rho 0.3, false-reassurance rate <= 20%%): %s** (false-reassurance rate at 0.3: %s).", if (r1) "MET: the cross-check is SENSITIVE" else "NOT MET", py_fmt_f(fr[["0.3"]], 2)),
           if (r1) "R2 not needed." else sprintf("R2: the smallest nominal rho of the grid with a false-reassurance rate <= 20%% is %s.", if (length(ok)) py_fmt_f(min(ok), 2) else "not reached up to 0.5"), "")
  }
  L <- c(L, "## A1. Per dataset (recovered excess and attenuation; mode \"detected\" unless stated)", "", "| dataset | mode | nominal rho | E primary | E recovered mean (SD) | attenuation | share with empirical P < 0.05 |", "|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(A))) { r <- A[i, ]; L <- c(L, sprintf("| %s | %s | %s | %s | %s (%s) | %s | %s |", r$tag, r$mode, py_fmt_f(r$rho, 2), fs(r$E_primary, 3), fs(r$E_recovered_mean, 3), py_fmt_f(r$E_recovered_sd, 3), py_fmt_f(r$attenuation, 2), py_fmt_f(r$share_empirical_P_lt_0_05, 2))) }
  if (nrow(R)) {
    n_ok <- sum(R$tag %in% HUMAN & R$E_RP_mean > 0.10)
    L <- c(L, "", "## A2. Ribosomal protein pairs (100 random sets of 11 pairs per dataset; excess E_RP over matched random pairs)", "", "| dataset | RP genes eligible | mean E_RP (SD) | E_RP in SE units of the consensus E | share of sets with E_RP >= 0.10 | mean rho of the RP pairs |", "|---|---|---|---|---|---|")
    for (i in seq_len(nrow(R))) { r <- R[i, ]; L <- c(L, sprintf("| %s | %d | %s (%s) | %s | %s | %s |", r$tag, r$rp_genes_eligible, fs(r$E_RP_mean, 3), py_fmt_f(r$E_RP_sd, 3), py_fmt_f(r$E_RP_in_SE_units, 1), py_fmt_f(r$share_E_RP_ge_0_10, 2), fs(r$mean_rho_RP_pairs, 3))) }
    L <- c(L, "", sprintf("**Reading rule (mean E_RP > +0.10 in at least 4 of the 5 human datasets): %s** (%d of 5).", if (n_ok >= 4) "MET" else "NOT MET", n_ok))
  }
  L <- c(L, "", sprintf("## A3. Margin sensitivity of the primary pooled interval [%s, %s] (the pre-set margin stays +/-0.10)", fs(e$lo, 3), fs(e$hi, 3)), "", "| margin | verdict of the primary interval |", "|---|---|", vapply(seq_len(nrow(M)), function(i) sprintf("| %s | %s |", py_fmt_f(M$margin[i], 3), M$verdict[i]), ""))
  L <- c(L, "", sprintf("The smallest margin at which the verdict is \"no meaningful excess\" is %s: the data exclude an excess in Fisher z above that value.", py_fmt_f(mmin, 3)))
  write_text_bytes(paste0(paste(L, collapse = .NL), .NL), file.path(OUT, "positive_control_report.md")); cat(paste(L, collapse = "\n"), "\n")
}
a <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) { if (!identical(a, "report")) for (t in (if (length(a)) a else TAGS)) run_dataset(t); report() }
