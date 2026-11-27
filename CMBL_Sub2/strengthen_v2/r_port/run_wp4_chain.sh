#!/bin/bash
## run_wp4_chain.sh - the WP4 analyses that follow run_datasets.R, one after the other (one heavy job at a time: this laptop has little free memory), all on the R results of run_datasets.R:
## meta-analysis, gene-graph null, metacell sensitivity (and its gene tables), positive control (spike-in, ribosomal control, margin), wider search (concordance), summary tables.
## Logs: results_R/wp4_chain.log (start / end / exit status of each step) and results_R/wp4_chain_<step>.log. Compare afterwards with compare_wp4.R meta | meta_graph | metacell | positive | wider.
## Usage: bash run_wp4_chain.sh [step ...]     steps: meta gene_graph metacell metacell_genes positive wider tables (default: all, in this order)
cd "$(dirname "$0")" || exit 1
R="${RSCRIPT:-/c/Program Files/R/R-4.3.1/bin/x64/Rscript.exe}"
export WP4_RES="$PWD/results_R/wp4/results" WP4_META="$PWD/results_R/wp4/meta" WP4_TAB_CAND="$PWD/results_R/wp4/perturbation_decisions"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a results_R/wp4_chain.log; }
run() { name="$1"; shift; log "START $name"; "$R" "$@" > "results_R/wp4_chain_${name}.log" 2>&1; log "END $name exit $?"; }
steps=("$@"); [ ${#steps[@]} -eq 0 ] && steps=(meta gene_graph metacell metacell_genes positive wider tables)
for s in "${steps[@]}"; do
  case "$s" in
    meta)           run meta wp4_meta_analysis.R ;;
    gene_graph)     run gene_graph wp4_islets/gene_graph_null.R ;;
    metacell)       run metacell wp4_islets/metacell.R ;;
    metacell_genes) run metacell_genes wp4_islets/metacell.R genes ;;
    positive)       run positive wp4_islets/positive_control.R ;;
    wider)          run wider wp4_islets/wider_search_analyze.R ;;
    tables)         run tables wp4_islets/make_tables.R ;;
    *)              log "unknown step $s" ;;
  esac
done
log "chain finished"
