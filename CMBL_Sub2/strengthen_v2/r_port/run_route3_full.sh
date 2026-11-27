#!/bin/bash
## run_route3_full.sh - 25 Nov 2026: waits for the running WP2 gate compare (gate_part1.R compare, log results_R/wp2/gate_compare.log), writes the gate report, then runs the
## full route 3 benchmark in R on the cached GEO files (results_R/wp2_data/bench_data; frozen list, 20 analysable per platform), summarizes it and compares it with the archive.
## GSE31851: if its own matrix download fails, RB_ARCHIVED + RB_ARCHIVED_ONLY make the runner take the archived line (listed in archived_records_used.txt, printed as MISSING by compare_route3.R). Log: results_R/route3_full.log
## 26 Nov 2026: RB_ARCHIVED_ONLY added (the runner now needs it); the run of 25 Nov did not use the fallback (GSE31851 downloaded and was computed by R).
cd "$(dirname "$0")" || exit 1
RS="${RSCRIPT:-/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe}"; LOG=results_R/route3_full.log; : > "$LOG"
until grep -q "^exit " results_R/wp2/gate_compare.log; do sleep 60; done
run() { echo "=== $(date '+%H:%M:%S') $*" >> "$LOG"; "$RS" "$@" 2>&1 | grep -v "Warning message\|built under" >> "$LOG"; echo "--- exit ${PIPESTATUS[0]}" >> "$LOG"; }
run wp2_benchmark2/gate_part1.R report
F="$PWD/results_R/benchmark_route3_full"; ARC="$(cd ../../benchmark_route3 && pwd -W)"
RB_OUT="$F" RB_ARCHIVED="$ARC/results_per_dataset.jsonl" RB_ARCHIVED_ONLY=GSE31851 run benchmark_route3/run_benchmark.R "$PWD/results_R/wp2_data/bench_data" 20
RB_RESULTS="$F/results_per_dataset.jsonl" RB_OUT="$F" run benchmark_route3/summarize_benchmark.R
RB_OUT="$F" run compare_route3.R
echo "=== byte comparisons with the archive" >> "$LOG"
for f in results_per_dataset.jsonl processing_log.csv per_dataset_results.csv summary_endpoints.json; do cmp -s "$F/$f" "../../benchmark_route3/$f" && echo "SAME  $f" >> "$LOG" || echo "DIFF  $f" >> "$LOG"; done
echo "=== $(date '+%H:%M:%S') route 3 full run finished" >> "$LOG"
