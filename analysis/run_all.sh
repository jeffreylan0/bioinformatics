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

# 4. the retina/vision subgraph of GO: curated roots plus their descendants
python3 analysis/retina_go_set.py \
  "${REF}/go-basic.obo" "${OUT}/retina_go_terms.tsv"

# 5. pick the terms to plot (top 8 BP, top 5 CC among retina terms), then draw
python3 analysis/select_terms.py \
  "${OUT}/go_enrichment_exprmatched.tsv" "${OUT}/go_terms_plot.tsv" \
  "BP:8,CC:5" "${OUT}/retina_go_terms.tsv"
# 6. the GO figure
python3 analysis/plot_go_bubble.py \
  "${OUT}/go_terms_plot.tsv" "${OUT}/go_bubble_plot"

# 7. why the photoreceptor terms rise: cell-population composition, then the
#    gene-level detail behind it
python3 analysis/plot_composition.py \
  "${ROOT}/1m-counts.tabular" "${ROOT}/1m-deseq.tabular" "${OUT}/cell_composition"
python3 analysis/plot_marker_genes.py \
  "${ROOT}/1m-deseq.tabular" "${OUT}/marker_genes"

echo "done -> ${OUT}/{go_bubble_plot,cell_composition,marker_genes}.svg + .pdf"
