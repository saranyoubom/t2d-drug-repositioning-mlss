#!/bin/bash
## run_o6_route3.sh - 26 Nov 2026, impact check O6 of the independent audit (audit_261126/AUDIT_REPORT.md): which route 3 series are analysable (>= 34 DE genes) when the ANOVA uses exact arithmetic
## for genes that are constant in all 12 samples (ANOVA_EXACT=1: p = 1 instead of numpy's p = 0). Dry runs only (DE counts and the frozen selection loop with its duplicate rule and 20 per platform),
## on the cached GEO files; series never screened in the archive are downloaded into the same cache. Two runs:
##   numpy   archived behaviour (gate: must equal the archived dry_run_per_dataset.jsonl in status and n_de for every series)
##   exact   ANOVA_EXACT=1
## Outputs: results_R/o6_anova_exact/route3_dry_{numpy,exact}/dry_run_per_dataset.jsonl, log results_R/o6_anova_exact/route3_dry.log. Nothing archived is written.
cd "$(dirname "$0")/.." || exit 1
RS="${RSCRIPT:-/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe}"; O=results_R/o6_anova_exact; mkdir -p "$O"; LOG=$O/route3_dry.log; : > "$LOG"
run() { echo "=== $(date '+%H:%M:%S') $*" >> "$LOG"; "$RS" "$@" 2>&1 | grep -v "Warning message\|built under" >> "$LOG"; echo "--- exit ${PIPESTATUS[0]}" >> "$LOG"; }
RB_OUT="$PWD/$O/route3_dry_numpy" run benchmark_route3/run_benchmark.R "$PWD/results_R/wp2_data/bench_data" 20 dry
ANOVA_EXACT=1 RB_OUT="$PWD/$O/route3_dry_exact" run benchmark_route3/run_benchmark.R "$PWD/results_R/wp2_data/bench_data" 20 dry
echo "=== $(date '+%H:%M:%S') finished" >> "$LOG"
