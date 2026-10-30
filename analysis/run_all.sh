#!/usr/bin/env bash
# End-to-end GO enrichment pipeline for the 1-month DDX41 cKO retina data.
# Usage: bash analysis/run_all.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${ROOT}/analysis/ref"
OUT="${ROOT}/results"
mkdir -p "$REF" "$OUT"

# 1. reference data (not committed; ~45 MB)
[ -f "${REF}/go-basic.obo" ] || curl -sS -o "${REF}/go-basic.obo" \
  https://current.geneontology.org/ontology/go-basic.obo
[ -f "${REF}/mgi.gaf.gz" ]   || curl -sS -o "${REF}/mgi.gaf.gz" \
  https://current.geneontology.org/annotations/mgi.gaf.gz

# 2. hypergeometric over-representation, up and down separately
python3 analysis/go_enrichment.py \
  "${ROOT}/1m-deseq.tabular" "${REF}/go-basic.obo" "${REF}/mgi.gaf.gz" "$OUT"

# 3. expression-matched permutation test (controls the expression confound)
python3 analysis/expression_matched_test.py \
  "${OUT}/go_enrichment_all.tsv" "${ROOT}/1m-deseq.tabular" \
  "${OUT}/gene2go_matrix.npz" "${OUT}/go_enrichment_exprmatched.tsv" 5000

# 4. pick the terms to plot, then draw
python3 analysis/select_terms.py \
  "${OUT}/go_enrichment_exprmatched.tsv" "${OUT}/go_terms_plot.tsv" 8
python3 analysis/plot_go_bubble.py \
  "${OUT}/go_terms_plot.tsv" "${OUT}/go_bubble_plot"

echo "done -> ${OUT}/go_bubble_plot.png"
