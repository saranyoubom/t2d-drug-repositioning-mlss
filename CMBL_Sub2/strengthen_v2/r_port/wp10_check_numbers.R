## wp10_check_numbers.R - R port of wp10_versions/check_wp10_numbers.py: independent rebuild of the WP10 numbers from the raw inputs with plain readers only (no import of the analysis scripts): the STRING classification (v12.0 export of the reported run,
## v12.5 network from the API), the DrugCentral counts and potencies (old file, new file), the MLSS combination counts, and the rank statements taken from mlss_comparison.csv; then compares them with the result files of the R run
## (results_R/wp10, or R_OUT_DIR) and with the numbers stated in the archived README.md of wp10_versions. Prints OK or FAIL per check, '%d checks, %d failures', and writes check_wp10_last_run.txt (into WP10_CHECK_OUT, default results_R/wp10) last.
## Usage: Rscript wp10_check_numbers.R
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(HERE, "pycompat.R")); suppressMessages(library(jsonlite))
V2 <- normalizePath(file.path(HERE, "..")); PROJECT <- normalizePath(file.path(V2, "..", "..")); R17 <- Sys.getenv("R17_DIR", normalizePath(file.path(V2, "..", "R_17genes"))); ARC <- file.path(V2, "wp10_versions")
RES <- Sys.getenv("R_OUT_DIR", file.path(HERE, "results_R", "wp10")); CHK <- Sys.getenv("WP10_CHECK_OUT", RES); dir.create(CHK, showWarnings = FALSE, recursive = TRUE)
count <- 0L; fails <- character(0); ok <- function(c, m) { count <<- count + 1L; cat(if (isTRUE(c)) "OK   " else "FAIL ", m, "\n", sep = ""); if (!isTRUE(c)) fails <<- c(fails, m) }
rb <- readBin(file.path(ARC, "README.md"), "raw", file.size(file.path(ARC, "README.md"))); rd <- rawToChar(rb); Encoding(rd) <- "UTF-8"; rd <- gsub("\r\n", "\n", rd, fixed = TRUE)
has <- function(t, m) ok(grepl(enc2utf8(t), rd, fixed = TRUE, useBytes = TRUE), sprintf("README states %s: \"%s\"", m, t))
pn <- function(a, b) vapply(seq_along(a), function(i) paste(py_sort(c(a[i], b[i])), collapse = "-"), "")
rdcsv <- function(f, sep = ",", enc = "UTF-8-BOM") read.csv(f, sep = sep, quote = "\"", colClasses = "character", na.strings = character(0), check.names = FALSE, fileEncoding = enc, comment.char = "", stringsAsFactors = FALSE)
## ---------------- STRING
## 27 Nov 2026: the retrieved network is read from WP10_RETRIEVED_DIR when set (run_all_v2.R: the frozen retrieval of the repository), else from RES as before
old <- rdcsv(file.path(R17, "string_interactions.csv")); nl <- strsplit(trimws(paste(readLines(file.path(Sys.getenv("WP10_RETRIEVED_DIR", RES), "string_v12_5_network.tsv"), encoding = "UTF-8", warn = FALSE), collapse = "\n")), "\n", fixed = TRUE)[[1]]; hdr <- strsplit(nl[1], "\t", fixed = TRUE)[[1]]
new <- as.data.frame(do.call(rbind, lapply(nl[-1], function(l) { p <- strsplit(l, "\t", fixed = TRUE)[[1]]; length(p) <- length(hdr); p })), stringsAsFactors = FALSE); names(new) <- hdr
cp <- rdcsv(file.path(R17, "Fig_4d_consensus_pairs.csv"))$Consensus_Pairs; cons <- unique(vapply(strsplit(cp, "_", fixed = TRUE), function(s) pn(s[1], s[2]), ""))
f4 <- read.csv(file.path(R17, "Fig_4e_Hes1_significant_correlations_only.csv"), header = FALSE, colClasses = "character", na.strings = character(0), fileEncoding = "UTF-8-BOM", check.names = FALSE); head <- unlist(f4[1, -1], use.names = FALSE); pairs <- list()
for (i in 2:nrow(f4)) for (j in seq_along(head)) { v <- f4[i, j + 1]; if (!(v %in% c("", "NA")) && f4[i, 1] != head[j] && abs(as.numeric(v)) >= 0.7) pairs[[pn(f4[i, 1], head[j])]] <- as.numeric(v) }
ok(length(pairs) == 120 && length(cons) == 21, "120 pairs with |r| >= 0.7 and 21 consensus pairs")
ALIAS <- c(Slc2a2 = "Glut2"); al <- function(x) ifelse(x %in% names(ALIAS), ALIAS[x], x)
sset <- function(a, b, alias) { if (alias) { a <- al(a); b <- al(b) }; unique(pn(a, b)) }
cls <- function(ss) { c <- c("Known Validated" = 0L, "Novel Discovery" = 0L, "STRING Only" = 0L, "Not Validated" = 0L); for (p in names(pairs)) { i <- p %in% ss; v <- p %in% cons; k <- if (i && v) "Known Validated" else if (v) "Novel Discovery" else if (i) "STRING Only" else "Not Validated"; c[[k]] <- c[[k]] + 1L }; c }
oa <- old$node1; ob <- old$node2; na <- new$preferredName_A; nb <- new$preferredName_B; res <- list()
for (nm in c("lit", "alias")) { alias <- nm == "alias"; res[[paste(nm, "v120")]] <- cls(sset(oa, ob, alias)); res[[paste(nm, "v125")]] <- cls(sset(na, nb, alias)) }
SM <- fromJSON(file.path(RES, "string_summary.json"), simplifyVector = FALSE)
for (z in list(c("lit", "as_reported"), c("alias", "alias"))) for (v in list(c("v120", "v12_0"), c("v125", "v12_5"))) { c <- res[[paste(z[1], v[1])]]; s <- SM[[z[2]]][[v[2]]]
  ok(all(vapply(names(c), function(k) c[[k]] == s[[k]], TRUE)) && c[["Known Validated"]] + c[["STRING Only"]] == s$in_STRING, sprintf("STRING %s %s counts equal the result file: {%s}", z[2], v[2], paste(sprintf("'%s': %d", names(c), c), collapse = ", "))) }
ok(length(sset(oa, ob, FALSE)) == 59 && length(sset(na, nb, FALSE)) == 59, "59 distinct pairs in the v12.0 export and 59 in the v12.5 network")
g17 <- c("Wnt5a", "Wnt5b", "Wnt9a", "Lef1", "Tcf7", "Tcf7l2", "Hes1", "Itpr1", "Ryr2", "Ryr3", "Cacna1c", "Cacna1d", "Kcnj11", "Ptbp1", "Glut2", "Ins1", "Glp1r")
e17 <- function(a, b) { a2 <- al(a); b2 <- al(b); k <- a2 %in% g17 & b2 %in% g17; unique(pn(a2[k], b2[k])) }; a0 <- e17(oa, ob); a1 <- e17(na, nb)
ok(length(a0) == 38 && length(a1) == 39 && identical(py_sort(setdiff(a1, a0)), c("Kcnj11-Ryr2", "Wnt5a-Wnt9a", "Wnt5b-Wnt9a")) && identical(py_sort(setdiff(a0, a1)), c("Cacna1c-Itpr1", "Cacna1d-Itpr1")), "edges among the 17 genes: 38 in 12.0, 39 in 12.5, 3 gained, 2 lost")
glut <- Filter(function(p) "Glut2" %in% strsplit(p, "-", fixed = TRUE)[[1]], names(pairs)); so_t <- sset(oa, ob, TRUE); so_f <- sset(oa, ob, FALSE)
ok(length(glut) == 15 && identical(py_sort(glut[glut %in% so_t]), c("Glp1r-Glut2", "Glut2-Ins1", "Glut2-Kcnj11", "Glut2-Tcf7l2")) && !any(glut %in% so_f), "Glut2: 15 pairs, 4 of them in STRING only under the alias, none under literal names")
has("29 (24.2%) | 70 (58.3%) | 34 (28.3%)", "reported STRING classification"); has("30 (25.0%) | 69 (57.5%) | 35 (29.2%)", "v12.5 classification"); has("33 (27.5%) | 66 (55.0%) | 38 (31.7%)", "v12.0 alias classification"); has("34 (28.3%) | 65 (54.2%) | 39 (32.5%)", "v12.5 alias classification")
has("38 in 12.0 and 39 in 12.5", "edges among 17 genes"); has("Glp1r-Glut2, Glut2-Ins1, Glut2-Kcnj11, Glut2-Tcf7l2; 15 of the 120 pairs involve Glut2", "Glut2 pairs")
## ---------------- DrugCentral
odc <- rdcsv(file.path(R17, "drug.target.interaction.csv"), ",", "UTF-8-BOM"); NDC <- file.path(PROJECT, "Database", "drug.target.interaction.tsv.gz"); nsha <- sha256_file(NDC)
ndc <- read.delim(gzfile(NDC, encoding = "UTF-8"), sep = "\t", quote = "\"", colClasses = "character", na.strings = character(0), check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
ok(startsWith(nsha, "19086849"), "new DrugCentral file sha256 starts with 19086849"); ok(nrow(odc) == 19378 && nrow(ndc) == 22364, "rows: 19378 before and 22364 after")
ok(length(unique(odc$DRUG_NAME)) == 2587 && length(unique(ndc$DRUG_NAME)) == 2897, "distinct drugs: 2587 before and 2897 after")
cen <- rdcsv(file.path(R17, "Fig_6b_Gene_Centrality_Summary.csv")); hubs <- toupper(cen$Gene[as.numeric(cen$Degree) >= 5]); ok(identical(py_sort(hubs), c("CACNA1D", "GLP1R", "PTBP1", "TCF7L2")), paste0("hubs (degree >= 5): ", paste(py_sort(hubs), collapse = ", ")))
targets <- function(rows) { k <- toupper(rows$GENE) %in% hubs & rows$ORGANISM == "Homo sapiens"; d <- rows$DRUG_NAME[k]; g <- toupper(rows$GENE[k]); t <- list(); for (i in seq_along(d)) if (!(g[i] %in% t[[d[i]]])) t[[d[i]]] <- c(t[[d[i]]], g[i]); t }
pot <- function(rows) { g <- toupper(rows$GENE); k <- g %in% hubs & rows$ACT_TYPE %in% c("IC50", "Ki", "EC50") & rows$ORGANISM == "Homo sapiens" & rows$RELATION == "=" & !(rows$ACT_VALUE %in% c("", "NA")); key <- paste(rows$DRUG_NAME[k], g[k], sep = "\r"); v <- as.numeric(rows$ACT_VALUE[k]); p <- list(); for (i in seq_along(key)) p[[key[i]]] <- max(if (is.null(p[[key[i]]])) -1e9 else p[[key[i]]], v[i]); p }
t0 <- targets(odc); t1 <- targets(ndc); DS <- fromJSON(file.path(RES, "drugcentral_summary.json"), simplifyVector = FALSE)
ok(length(t0) == 18 && length(t1) == 19 && identical(py_sort(setdiff(names(t1), names(t0))), "tirzepatide") && !length(setdiff(names(t0), names(t1))), "MLSS drugs 18 before, 19 after, only tirzepatide added, none removed")
cnt <- function(t, g) sum(vapply(t, function(gs) g %in% gs, TRUE))
ok(identical(c(cnt(t0, "GLP1R"), cnt(t1, "GLP1R"), cnt(t0, "CACNA1D"), cnt(t1, "CACNA1D"), cnt(t0, "TCF7L2") + cnt(t1, "TCF7L2") + cnt(t0, "PTBP1") + cnt(t1, "PTBP1")), c(7L, 8L, 11L, 11L, 0L)), "GLP1R drugs 7 then 8, CACNA1D drugs 11 then 11, none for TCF7L2 and PTBP1")
ok(all(vapply(c(t0, t1), length, 1L) == 1L), "every drug has exactly one hub target")
n0 <- length(t0); n1 <- length(t1); cross <- function(t) { cb <- combn(length(t), 2); sum(apply(cb, 2, function(ix) !identical(t[[ix[1]]], t[[ix[2]]]))) }
ok(identical(c((n0 * (n0 - 1L)) %/% 2L, (n1 * (n1 - 1L)) %/% 2L, as.integer(cross(t0)), as.integer(cross(t1))), c(153L, 171L, 77L, 88L)), "combinations 153 then 171; cross-target 77 then 88")
p0 <- pot(odc); p1 <- pot(ndc); keys <- union(names(p0), names(p1)); in_t <- vapply(strsplit(keys, "\r", fixed = TRUE), function(k) k[1] %in% union(names(t0), names(t1)), TRUE)
gv <- function(p, k) if (is.null(p[[k]])) -99 else p[[k]]; dk <- py_sort(keys[in_t & vapply(keys, function(k) abs(gv(p0, k) - gv(p1, k)) > 1e-9, TRUE)])
ok(identical(dk, "tirzepatide\rGLP1R") && abs(p1[["tirzepatide\rGLP1R"]] - 8.37) < 1e-9, "the only potency that differs between the files is tirzepatide at GLP1R = 8.37")
for (kv in list(list("mlss_drugs_old", 18), list("mlss_drugs_new", 19), list("combinations_old", 153), list("combinations_new", 171), list("cross_target_old", 77), list("cross_target_new", 88), list("combinations_shared", 153), list("score_changed_shared", 0), list("top15_old_still_in_top15_new", 14),
                list("combinations_with_new_drugs", 18), list("pic50_changes", 1), list("rows_old", 19378), list("rows_new", 22364))) ok(DS[[kv[[1]]]] == kv[[2]], sprintf("drugcentral_summary.json %s = %s", kv[[1]], kv[[2]]))
ok(DS$spearman_scores_shared == 1 && DS$spearman_ranks_shared == 1 && !isTRUE(DS$tie_straddles_rank15_old) && !isTRUE(DS$tie_straddles_rank15_new), "Spearman 1.0 for scores and ranks of the shared combinations; no tie across rank 15")
cm <- rdcsv(file.path(RES, "mlss_comparison.csv"), ",", "UTF-8"); rownames(cm) <- cm$combination; ro <- suppressWarnings(as.integer(cm$rank_old)); top <- cm[!is.na(ro) & ro <= 15, ]; top <- top[order(as.integer(top$rank_old)), ]
ok(nrow(top) == 15 && sum(as.integer(top$rank_new) <= 15) == 14, "14 of the 15 reported top combinations are in the new top 15")
ok(all(as.integer(top$rank_new)[as.integer(top$rank_old) <= 12] == as.integer(top$rank_old)[as.integer(top$rank_old) <= 12]) && identical(as.integer(top$rank_new)[13:15], c(14L, 15L, 16L)), "old ranks 1 to 12 unchanged; old ranks 13, 14, 15 become 14, 15, 16")
it <- cm["isradipine + tirzepatide", ]; ok(it$rank_new == "13" && it$status == "new only" && abs(as.numeric(it$MLSS_new) - 9.9158) < 5e-5, "isradipine + tirzepatide enters at rank 13 with MLSS 9.9158")
tz <- sort(as.integer(cm$rank_new[grepl("tirzepatide", cm$combination, fixed = TRUE)])); ok(length(tz) == 18 && identical(tz[1:11], c(13L, 26L, 27L, 30L, 33L, 69L, 70L, 71L, 72L, 73L, 74L)) && min(tz[12:length(tz)]) == 99, "tirzepatide combinations: 18; cross-target ranks 13, 26, 27, 30, 33, 69 to 74; same-target from rank 99")
sm <- cm[cm$status == "both" & cm$MLSS_old != "" & abs(as.numeric(cm$MLSS_old) - as.numeric(cm$MLSS_new)) > 1e-9, ]; ok(nrow(sm) == 0, "no shared combination changed its score")
for (tw in list(c("Spearman rho 1.0", "rank correlation"), c("isradipine + tirzepatide enters at rank 13 (MLSS 9.9158)", "new top-15 entrant"), c("**14 of the 15 reported combinations stay in the top 15**", "top-15 overlap"), c("153 become 171; cross-target pairs 77 become 88", "combination counts"),
                c("18 become 19", "drug count"), c("(8.37)", "tirzepatide potency"), c("22364 rows against 19378", "rows"), c("2897 distinct drugs against 2587", "distinct drugs"),
                c("11 cross-target at ranks 13, 26, 27, 30, 33, 69, 70, 71, 72, 73 and 74; 7 same-target at rank 99 or lower", "tirzepatide combination ranks"), c("0 changed", "changed scores"))) has(tw[1], tw[2])
fpl <- readLines(file.path(ARC, "FINGERPRINTS.txt"), warn = FALSE, encoding = "UTF-8"); ok(length(fpl) == 10, "FINGERPRINTS.txt holds 5 fingerprints")
ok(any(grepl(sha256_file(file.path(ARC, "PRESPEC_WP10_database_versions.md")), fpl, fixed = TRUE)), "the pre-specification equals its fingerprint (unchanged since it was fingerprinted)")
cat(sprintf("\n%d checks, %d failures\n", count, length(fails)))
write_text_bytes(gsub("\n", .NL, sprintf("%d checks, %d failures%s\n", count, length(fails), if (length(fails)) paste0(": ", paste(fails, collapse = "; ")) else ""), fixed = TRUE), file.path(CHK, "check_wp10_last_run.txt"))
