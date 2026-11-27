#!/bin/bash
# run_post_hoc.sh - the three post hoc runs of the benchmark (post_hoc_191126.R), one after the other (one heavy job at a time on the 8 GB laptop):
#   dev         development series, the 74 matrices of the R extract (results_R/wp2/compact/dev_sets.json, GSE31851 from the archived matrices)
#   ins1        the INS-1 matrix of the manuscript, 10,000 permutations
#   validation  the 59 validation sets rebuilt from the cached GEO files (results_R/wp2_data)
# Usage: bash run_post_hoc.sh   (logs in results_R/post_hoc/run_*.log; results in results_R/post_hoc/results_<mode>; compare with compare_posthoc.R)
cd "$(dirname "$0")"
R="/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe"
LOG=results_R/post_hoc; mkdir -p $LOG
echo "$(date +%T) dev"; WP2_DEV_SETS="$PWD/results_R/wp2/compact/dev_sets.json" "$R" post_hoc_191126.R dev > $LOG/run_dev.log 2>&1
echo "$(date +%T) ins1"; "$R" post_hoc_191126.R ins1 > $LOG/run_ins1.log 2>&1
echo "$(date +%T) validation"; "$R" post_hoc_191126.R validation > $LOG/run_validation.log 2>&1
echo "$(date +%T) done"
