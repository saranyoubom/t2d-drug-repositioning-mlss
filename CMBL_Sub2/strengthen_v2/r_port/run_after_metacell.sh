#!/bin/bash
## run_after_metacell.sh - queued checks of 25 Nov 2026: waits for the running metacell chain (run_wp4_chain.sh), then runs one after the other: metacell comparison, README check on the R results, a regression rerun of the light analyses after the
## number-formatter fix (WP5a, WP9, WP10, WP11, WP12, WP13, WP14), the AF6 workbook built from the R results (compared with the submitted one), and the agreement report. Log: results_R/after_metacell.log
cd "$(dirname "$0")" || exit 1
RS="${RSCRIPT:-/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe}"; LOG=results_R/after_metacell.log; : > "$LOG"
until awk '/START metacell$/{s=NR; d=0} /chain finished/{if (s) d=1} END{exit !d}' results_R/wp4_chain.log; do sleep 60; done
run() { echo "=== $(date '+%H:%M:%S') $*" >> "$LOG"; "$RS" "$@" 2>&1 | grep -v "Warning message\|built under" >> "$LOG"; echo "--- exit ${PIPESTATUS[0]}" >> "$LOG"; }
run compare_wp4.R metacell
CHECK_SOURCE=R run wp4_islets/check_readme_numbers.R
## regression of the formatter fix
ES_OUT="$PWD/results_R/wp5" ES_LABEL=effect_sizes.py run effect_sizes.R
run wp9_safety.R analyze; run wp9_safety.R fallback-analyze; run wp9_check_numbers.R
run wp10_string_sensitivity.R analyze; run wp10_drugcentral_sensitivity.R; run wp10_string_web_export_check.R; run wp10_check_numbers.R
run wp12_annotations.R analyze; run wp13_gapdh_versions.R; run wp14_effsens_analyzed_matrix.R
for s in reconstruct_threshold check_threshold150 check_refA_vs_plate1 check_aug26_nov09_order export_tables fingerprint_inputs check_wp11_numbers; do run wp11_plate_provenance/$s.R; done
echo "=== byte comparisons of the check outputs with the archive" >> "$LOG"
A=..
for f in threshold_fit.txt threshold150_matches.csv threshold150_summary.txt refA_vs_plate1.txt aug26_nov09_order_check.txt aug26_nov09_order_check.json aug26_workbook_vs_instrument.csv threshold_validation.csv refA_block_residuals.csv refA_dcq_comparison.csv; do cmp -s results_R/wp11/$f $A/wp11_plate_provenance/results/$f && echo "SAME  wp11 $f" >> "$LOG" || echo "DIFF  wp11 $f" >> "$LOG"; done
for p in "results_R/wp11/check_wp11_last_run.txt wp11_plate_provenance/check_wp11_last_run.txt" "results_R/wp9/check_wp9_last_run.txt wp9_safety/check_wp9_last_run.txt" "results_R/wp10/check_wp10_last_run.txt wp10_versions/check_wp10_last_run.txt" "results_R/wp10/string_web_export_check.json wp10_versions/results/string_web_export_check.json"; do
  set -- $p; cmp -s "$1" "$A/$2" && echo "SAME  $1" >> "$LOG" || echo "DIFF  $1" >> "$LOG"; done
## Additional file 6 from the R results
AF6_SOURCE=R AF6_OUT="$PWD/results_R/af6/Additional file 6 from R results.xlsx" run af6/build_af6.R
AF6_OUT="$PWD/results_R/af6/Additional file 6 from R results.xlsx" run compare_af6.R 4
run agreement_report.R
echo "=== $(date '+%H:%M:%S') all queued checks finished" >> "$LOG"
