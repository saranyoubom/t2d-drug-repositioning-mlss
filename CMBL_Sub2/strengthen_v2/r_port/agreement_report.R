## agreement_report.R - the agreement report of R_PORT_PLAN.md step 5: runs every compare script of the port (each compares the R results in results_R/ with the saved Python results of the repository) and writes
## results_R/agreement_report.md (per package: items agreeing, files byte-identical, items differing, values compared, largest relative difference, check of the line counts; every line that is not an agreement, with the
## reason known for it), results_R/agreement_summary.csv and results_R/agreement_lines.txt (every line that is not an agreement, as 'package: line', the form in which the EXPLAIN keys below are written).
## Gate rules, tightened 26 Nov 2026 after the independent audit (audit_261126/AUDIT_REPORT.md, C1):
##  1. a compare script that exits with a nonzero status gives an ERROR line (whatever it would have printed after the failure is lost);
##  2. only DIFFER lines can carry a known reason; MISSING, ERROR, FAIL, COUNT and mismatch lines are always UNEXPLAINED;
##  3. a known reason belongs to one exact item: the package, the whole DIFFER line (counts, largest difference) and the detail lines the compare script prints under it (which columns, rows or lines differ, with
##     their first values), so any change of a differing item comes back as UNEXPLAINED and has to be looked at again; the reasons are kept in agreement_known_reasons.tsv (columns key and reason);
##  4. each package must print the expected numbers of AGREE, SAME and DIFFER lines (EXPECT, from the reference run); any other number gives a COUNT line: an item that disappeared, a crash half-way, a new item.
## Usage: Rscript agreement_report.R [package ...]   (default: all; names as in the first column of RUNS; a run of some packages writes a report of those packages only)
R_PORT <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))); source(file.path(R_PORT, "pycompat.R"))
RS <- Sys.getenv("RSCRIPT", file.path(R.home("bin"), "Rscript.exe")); if (!file.exists(RS)) RS <- file.path(R.home("bin"), "Rscript")
RUNS <- list(                                                                                                                         # package, compare script, arguments[, environment variables for the script]
  list("WP1", "compare_wp1.R", character(0)), list("WP3", "compare_wp3.R", character(0)), list("WP5a", "compare_wp5.R", character(0)), list("WP9", "compare_wp9.R", character(0)), list("WP10", "compare_wp10.R", character(0)), list("WP12", "compare_wp12.R", character(0)),
  list("WP13", "compare_wp13.R", character(0)), list("WP14", "compare_wp14.R", character(0)), list("reanalysis_17genes", "compare_reanalysis.R", character(0)),
  list("route 3 (full run, 58 series)", "compare_route3.R", character(0), c(RB_OUT = file.path(R_PORT, "results_R", "benchmark_route3_full"))), list("route 3 selection", "compare_select.R", character(0)),
  list("WP2 post-processing", "compare_wp2.R", c("post", "design")), list("WP2 gate (extract)", "compare_wp2.R", "gate"), list("WP2 held-out arrays", "compare_wp2.R", "heldout"), list("WP2 RNA-seq", "compare_wp2.R", "rnaseq"), list("WP2 calibration", "compare_wp2.R", "calib"),
  list("WP2 analyse_matrix on 74 matrices", "compare_wp2.R", "core"), list("post hoc checks", "compare_posthoc.R", character(0)),
  list("WP4 datasets", "compare_wp4.R", "results"), list("WP4 meta-analysis", "compare_wp4.R", "meta"), list("WP4 gene-graph null", "compare_wp4.R", "meta_graph"), list("WP4 metacells", "compare_wp4.R", "metacell"), list("WP4 positive control", "compare_wp4.R", "positive"),
  list("WP4 GSE60803", "compare_wp4.R", "perturbation"), list("WP4 wider search", "compare_wp4.R", "wider"), list("WP4 live queries", "compare_wp4.R", "scout"), list("Additional file 6 (built from the archived results)", "compare_af6.R", character(0)))
## expected numbers of AGREE/SAME/DIFFER lines per package (reference run of 26 Nov 2026, 07:40, after the audit fixes; change only after looking at why a count changed)
EXPECT <- c("WP1" = "12/0/2", "WP3" = "4/0/0", "WP5a" = "2/0/1", "WP9" = "9/9/0", "WP10" = "8/0/0", "WP12" = "5/5/0", "WP13" = "2/0/0", "WP14" = "1/0/0", "reanalysis_17genes" = "14/0/0",
            "route 3 (full run, 58 series)" = "58/0/3", "route 3 selection" = "1/0/1", "WP2 post-processing" = "12/0/0", "WP2 gate (extract)" = "3/0/0", "WP2 held-out arrays" = "0/0/2", "WP2 RNA-seq" = "1/0/1",
            "WP2 calibration" = "1/0/5", "WP2 analyse_matrix on 74 matrices" = "0/0/1", "post hoc checks" = "8/0/6", "WP4 datasets" = "69/0/1", "WP4 meta-analysis" = "8/1/0", "WP4 gene-graph null" = "2/1/0",
            "WP4 metacells" = "4/0/0", "WP4 positive control" = "14/1/2", "WP4 GSE60803" = "2/0/0", "WP4 wider search" = "5/1/0", "WP4 live queries" = "3/8/0", "Additional file 6 (built from the archived results)" = "1/1/0")
## known reasons: agreement_known_reasons.tsv, tab-separated columns 'key' (an exact item as written in results_R/agreement_lines.txt) and 'reason' (checked for that item; see R_PORT_PLAN.md, the memory notes, audit_261126/AUDIT_REPORT.md)
KR_FILE <- file.path(R_PORT, "agreement_known_reasons.tsv")
KR <- if (file.exists(KR_FILE)) read.delim(KR_FILE, quote = "", comment.char = "", stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character") else data.frame(key = character(0), reason = character(0))
lines <- list(); summ <- list(); dump <- character(0); sel <- commandArgs(trailingOnly = TRUE)
for (run in RUNS) {
  pk <- run[[1]]; if (length(sel) && !(pk %in% sel)) next; cat("== ", pk, "\n", sep = ""); t0 <- Sys.time()
  env <- if (length(run) >= 4) run[[4]] else character(0); old <- if (length(env)) Sys.getenv(names(env), unset = NA, names = TRUE) else character(0)
  if (length(env)) do.call(Sys.setenv, as.list(env))
  out <- tryCatch(suppressWarnings(system2(RS, c("--vanilla", file.path(R_PORT, run[[2]]), run[[3]]), stdout = TRUE, stderr = FALSE)), error = function(e) structure(paste("ERROR", conditionMessage(e)), status = -1L))
  for (k in names(env)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, setNames(list(old[[k]]), k))
  st <- attr(out, "status"); out <- as.character(out)
  keep <- grepl("^(AGREE|SAME|DIFFER|BYTES|MISSING|PASS|FAIL|ERROR)\\b", out) | grepl("^\\s+mismatch\\b", out); idx <- which(keep); L <- trimws(out[idx])
  D <- vapply(idx, function(i) { j <- i + 1L; det <- character(0)                                                     # the indented detail lines printed under an item: part of its key
    while (j <= length(out) && grepl("^\\s+\\S", out[j]) && !grepl("^\\s+mismatch\\b", out[j])) { det <- c(det, trimws(out[j])); j <- j + 1L }; paste(det, collapse = " // ") }, "")
  if (!is.null(st) && !identical(as.integer(st), 0L)) { L <- c(L, sprintf("ERROR %s exited with status %s: whatever it would have printed after the failure is missing", paste(c(run[[2]], run[[3]]), collapse = " "), st)); D <- c(D, "") }
  agr <- L[startsWith(L, "AGREE")]; mv <- regmatches(agr, regexec("([0-9]+) values compared", agr)); nv <- sum(vapply(mv, function(z) if (length(z) > 1) as.numeric(z[2]) else 0, 1))
  rel <- suppressWarnings(as.numeric(sub("^.*largest relative difference ([0-9.eE+-]+).*$", "\\1", L[grepl("largest relative difference", L)]))); mx <- if (any(!is.na(rel))) max(rel, na.rm = TRUE) else NA_real_
  nA <- sum(startsWith(L, "AGREE")); nS <- sum(startsWith(L, "SAME")); nD <- sum(startsWith(L, "DIFFER")); nM <- sum(startsWith(L, "MISSING")); nE <- sum(startsWith(L, "ERROR")) + sum(startsWith(L, "FAIL")) + sum(startsWith(L, "mismatch"))
  obs <- sprintf("%d/%d/%d", nA, nS, nD); exp_n <- if (pk %in% names(EXPECT)) EXPECT[[pk]] else NA_character_
  cnt <- if (is.na(exp_n)) sprintf("no expected numbers (observed %s)", obs) else if (obs == exp_n) "as expected" else sprintf("observed %s, expected %s", obs, exp_n)
  if (cnt != "as expected") { L <- c(L, sprintf("COUNT AGREE/SAME/DIFFER lines: %s", cnt)); D <- c(D, "") }
  summ[[length(summ) + 1L]] <- data.frame(package = pk, compare = paste(c(run[[2]], run[[3]]), collapse = " "), items_agree = nA, files_byte_identical = nS, items_differ = nD, missing = nM, errors = nE, line_counts = cnt,
                                          values_compared = nv, largest_relative_difference = mx, seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), stringsAsFactors = FALSE)
  for (k in which(grepl("^(DIFFER|MISSING|ERROR|FAIL|COUNT|mismatch)\\b", L))) {
    key <- paste0(pk, ": ", L[k], if (nzchar(D[k])) paste0(" // ", D[k]) else ""); why <- if (startsWith(L[k], "DIFFER")) KR$reason[match(key, KR$key)] else NA_character_
    lines[[length(lines) + 1L]] <- c(pk, L[k], if (is.na(why)) "UNEXPLAINED" else why); dump <- c(dump, key)
  }
  cat(sprintf("   agree %d, byte-identical %d, differ %d, missing %d, error %d, line counts %s; %d values compared, largest relative difference %.2g\n", nA, nS, nD, nM, nE, cnt, nv, mx))
}
S <- do.call(rbind, summ); py_write_csv_repr(S, file.path(R_PORT, "results_R", "agreement_summary.csv"), eol = "\r\n")
write_text_bytes(paste0(paste(dump, collapse = "\n"), "\n"), file.path(R_PORT, "results_R", "agreement_lines.txt"))
unexpl <- sum(vapply(lines, function(l) l[3] == "UNEXPLAINED", TRUE)); ncnt <- sum(S$line_counts != "as expected")
md <- c("# Agreement of the R port with the saved Python results", "",
        sprintf("Made by `agreement_report.R` on %s. Each row runs one compare script against the archived Python results of `strengthen_v2/`.", format(Sys.time(), "%d %B %Y, %H:%M")), "",
        "- An item (a file, or one record where a script compares record by record, as the route 3 script does per series) agrees when every number is equal or |a - b| <= 1e-12 + tol * max(|a|, |b|) (tol 2e-6 for CSV files and 2e-9 for JSON files unless the script names another), and every string and flag is equal.",
        "- Values compared: the values of the agreeing items only. Largest relative difference: |a - b| / max(|a|, |b|) over the values that differ by more than 1e-12, in agreeing and differing items.",
        "- Gate rules (26 Nov 2026): a compare script that exits with an error gives an ERROR line; only a DIFFER line can carry a known reason, and only the exact item that was checked (the line and the detail lines under it, `agreement_known_reasons.tsv`); each package must print its expected numbers of AGREE / SAME / DIFFER lines, otherwise a COUNT line is added. MISSING, ERROR, FAIL, COUNT and mismatch lines never carry a known reason.", "",
        "| Package | Items agreeing | Files byte-identical | Items differing | Missing | Errors | Line counts | Values compared | Largest relative difference |", "|---|---:|---:|---:|---:|---:|---|---:|---:|")
for (i in seq_len(nrow(S))) md <- c(md, sprintf("| %s | %d | %d | %d | %d | %d | %s | %s | %s |", S$package[i], S$items_agree[i], S$files_byte_identical[i], S$items_differ[i], S$missing[i], S$errors[i], S$line_counts[i],
                                                format(S$values_compared[i], big.mark = ","), if (is.na(S$largest_relative_difference[i])) "" else sprintf("%.1e", S$largest_relative_difference[i])))
md <- c(md, sprintf("| **All** | %d | %d | %d | %d | %d | %d not as expected | %s | %s |", sum(S$items_agree), sum(S$files_byte_identical), sum(S$items_differ), sum(S$missing), sum(S$errors), ncnt,
                    format(sum(S$values_compared), big.mark = ","), sprintf("%.1e", max(S$largest_relative_difference, na.rm = TRUE))), "", "## Lines that are not agreements, with the reason", "")
if (length(lines)) { md <- c(md, "| Package | Line of the compare script | Reason |", "|---|---|---|"); for (l in lines) md <- c(md, sprintf("| %s | %s | %s |", l[1], gsub("\\|", "/", substr(l[2], 1, 230)), l[3])) } else md <- c(md, "None.")
md <- c(md, "", sprintf("Lines without a known reason: %d. Packages whose line counts are not as expected: %d.", unexpl, ncnt))
write_text_bytes(paste0(paste(md, collapse = "\n"), "\n"), file.path(R_PORT, "results_R", "agreement_report.md")); cat("agreement report written; lines without a known reason:", unexpl, "; packages with unexpected line counts:", ncnt, "\n")
