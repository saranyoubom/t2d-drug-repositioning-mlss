#!/bin/bash
# run_o6_calib.sh - 27 Nov 2026, impact check O6, steps 3 to 7 (run after gate_part1.R extract with ANOVA_EXACT=1 and prep_exact_calib.R prep; analysis only, nothing in the package changes).
# All outputs go to results_R/o6_anova_exact (WP2_OUT=.../wp2_exact for the development series); results_R/wp2 and the archive are only read.
#   3 gate_part1.R consistency   analyse_matrix on the sets whose matrix changed (the other lines were copied by prep), then prep_exact_calib.R order
#   4 calibrate_t2.R run         the series without a copied calibration file (2 workers), then summarize (WP2_DEV_CORE = the exact-rule core results)
#   5 apply_dev.R                development application on the exact-rule core results
#   6 run_heldout.R arrays full  GSE68424 alone, numpy order (gate: must equal the line of the R run) and exact rule
#   7 summarize_validation.R     development (exact) + held-out arrays of the R run with the GSE68424 line replaced + RNA-seq of the R run (unchanged, no constant genes)
cd "$(dirname "$0")/.."
R="/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe"
O="$PWD/results_R/o6_anova_exact"; W="$O/wp2_exact"; RUN="$PWD/results_R/wp2"
set -o pipefail
echo "$(date +%T) consistency"; ANOVA_EXACT=1 WP2_OUT="$W" "$R" wp2_benchmark2/gate_part1.R consistency > "$O/calib_consistency.log" 2>&1 || { echo "consistency failed"; exit 1; }
"$R" o6_anova_exact/prep_exact_calib.R order >> "$O/calib_consistency.log" 2>&1 || { echo "order failed"; exit 1; }
echo "$(date +%T) calibration"; for w in 1 2; do ANOVA_EXACT=1 WP2_OUT="$W" "$R" wp2_benchmark2/calibrate_t2.R run $w 2 > "$O/calib_run_w$w.log" 2>&1 & done; wait
ANOVA_EXACT=1 WP2_OUT="$W" WP2_DEV_CORE="$W/compact/dev_core_results.jsonl" "$R" wp2_benchmark2/calibrate_t2.R summarize > "$O/calib_summary.log" 2>&1 || { echo "calibration summary failed"; exit 1; }
echo "$(date +%T) apply_dev"; WP2_OUT="$W" WP2_DEV_CORE="$W/compact/dev_core_results.jsonl" "$R" wp2_benchmark2/apply_dev.R > "$O/apply_dev.log" 2>&1 || { echo "apply_dev failed"; exit 1; }
echo "$(date +%T) heldout GSE68424"
WP2_ONLY=GSE68424 WP2_OUT="$O/heldout68424_numpy" "$R" wp2_benchmark2/run_heldout.R arrays full > "$O/heldout68424_numpy.log" 2>&1 || { echo "heldout numpy failed"; exit 1; }
ANOVA_EXACT=1 WP2_ONLY=GSE68424 WP2_OUT="$O/heldout68424_exact" "$R" wp2_benchmark2/run_heldout.R arrays full > "$O/heldout68424_exact.log" 2>&1 || { echo "heldout exact failed"; exit 1; }
"$R" o6_anova_exact/merge_heldout.R > "$O/merge_heldout.log" 2>&1 || { echo "merge failed"; exit 1; }
echo "$(date +%T) summarize_validation"; WP2_OUT="$W" WP2_DEV_CSV="$W/dev_application_per_series.csv" WP2_HELDOUT="$O/heldout_arrays_results_exact.jsonl" WP2_RNASEQ="$RUN/heldout_rnaseq_results.jsonl" "$R" wp2_benchmark2/summarize_validation.R > "$O/summarize_validation.log" 2>&1 || { echo "summarize_validation failed"; exit 1; }
#   8 baseline for the comparison: the same two summaries on the R run with numpy's order (results_R/wp2 inputs; the archived Python files differ from them only in the BLAS-noise numbers of
#     the constant-gene matrices, R_PORT_PLAN.md 25 Nov); the calibration baseline is results_R/wp2/calibration_verdict.txt
B="$O/wp2_numpy_R"; mkdir -p "$B"
echo "$(date +%T) baseline"; WP2_OUT="$B" WP2_DEV_CORE="$RUN/compact/dev_core_results.jsonl" "$R" wp2_benchmark2/apply_dev.R > "$O/baseline_apply_dev.log" 2>&1 || { echo "baseline apply_dev failed"; exit 1; }
WP2_OUT="$B" WP2_DEV_CSV="$B/dev_application_per_series.csv" WP2_HELDOUT="$RUN/heldout_arrays_results.jsonl" WP2_RNASEQ="$RUN/heldout_rnaseq_results.jsonl" "$R" wp2_benchmark2/summarize_validation.R > "$O/baseline_summarize_validation.log" 2>&1 || { echo "baseline summarize_validation failed"; exit 1; }
echo "$(date +%T) done"
