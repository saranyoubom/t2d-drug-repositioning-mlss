#!/bin/bash
# run_wp2_validation.sh - the data-dependent verification runs of the WP2 port, one after the other (one heavy job at a time on the 8 GB laptop):
#   1 gate_part1.R extract         re-read the 58 part-1 series from the cached GEO files, compare status, DE counts, gene lists with the archive
#   2 gate_part1.R consistency     (already run from the archived matrices; skipped here)
#   3 run_heldout.R arrays dry|full   held-out arrays
#   4 run_rnaseq.R dry|full           held-out RNA-seq (human; mouse has no NCBI-generated counts)
# Usage: bash run_wp2_validation.sh   (waits until wp2_prefetch.R has finished; logs in results_R/wp2/validation_*.log)
cd "$(dirname "$0")"
R="/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe"
LOG=results_R/wp2
while ! grep -q "prefetch finished" results_R/wp2_data/prefetch_log_retry.txt 2>/dev/null; do sleep 20; done
echo "$(date +%T) extract"; WP2_DEV_SETS="$PWD/oracle/dev_sets_archive.json" "$R" wp2_benchmark2/gate_part1.R extract > $LOG/validation_extract.log 2>&1
echo "$(date +%T) heldout dry"; "$R" wp2_benchmark2/run_heldout.R arrays dry > $LOG/validation_heldout_dry.log 2>&1
echo "$(date +%T) heldout full"; "$R" wp2_benchmark2/run_heldout.R arrays full > $LOG/validation_heldout_full.log 2>&1
echo "$(date +%T) rnaseq dry"; "$R" wp2_benchmark2/run_rnaseq.R dry > $LOG/validation_rnaseq_dry.log 2>&1
echo "$(date +%T) rnaseq full"; "$R" wp2_benchmark2/run_rnaseq.R full > $LOG/validation_rnaseq_full.log 2>&1
echo "$(date +%T) done"
