#!/bin/bash
## run_o6_heldout.sh - 26 Nov 2026, impact check O6 (see run_o6_route3.sh): dry runs of the WP2 held-out arrays on the ARCHIVED frozen held-out list, numpy order (gate against the archived
## heldout_arrays_dry.jsonl) and ANOVA_EXACT=1; then the RNA-seq constant-gene count and compare_o6.R. Series never screened in the archive are downloaded into the usual cache.
## The held-out list itself excludes series that share samples with the development series; if the exact rule changes the development set, that exclusion is checked separately.
cd "$(dirname "$0")/.." || exit 1
RS="${RSCRIPT:-/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe}"; O=results_R/o6_anova_exact; mkdir -p "$O"; LOG=$O/heldout_dry.log; : > "$LOG"
run() { echo "=== $(date '+%H:%M:%S') $*" >> "$LOG"; "$RS" "$@" 2>&1 | grep -v "Warning message\|built under" >> "$LOG"; echo "--- exit ${PIPESTATUS[0]}" >> "$LOG"; }
WP2_OUT="$PWD/$O/heldout_numpy" run wp2_benchmark2/run_heldout.R arrays dry
ANOVA_EXACT=1 WP2_OUT="$PWD/$O/heldout_exact" run wp2_benchmark2/run_heldout.R arrays dry
run o6_anova_exact/rnaseq_constant_check.R
run o6_anova_exact/compare_o6.R
echo "=== $(date '+%H:%M:%S') finished" >> "$LOG"
