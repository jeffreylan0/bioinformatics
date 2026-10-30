# GO enrichment — 1-month CHX10-Cre; Ddx41^fl/fl retina

Gene-ontology over-representation analysis of `1m-deseq.tabular`, with the result
figure in [`../results/go_bubble_plot.png`](../results/go_bubble_plot.png).

Run everything with `bash analysis/run_all.sh`.

## Input

`1m-deseq.tabular` is a headerless DESeq2 table: `gene, baseMean, log2FoldChange,
lfcSE, stat, pvalue, padj`, 26,264 rows. `1m-counts.tabular` holds DESeq2-normalized
counts for 3 control and 3 cKO retinas.

The sign convention was **checked against the counts matrix, not assumed**: DESeq2
`log2FoldChange` reproduces `log2(mean cKO / mean WT)` essentially exactly (r = 0.95
over all genes; identical to 3 decimals for individual genes such as `Abcd4` +3.00 and
`Pde6c` -1.20). **Positive log2FoldChange = higher in the cKO.**

`Ddx41` itself is unchanged at the whole-gene level (log2FC 0.15, padj 0.40) — expected
for a floxed-exon conditional allele in a bulk tissue where Cre is active in only part
of the retina.

## Method

| Step | Choice |
|---|---|
| Significance | DESeq2 `padj < 0.05`, split by sign of log2FoldChange |
| Annotations | GO mouse GAF from MGI, `taxon:10090`, NOT-qualified and `ND` evidence dropped; all other evidence codes (including IEA) kept |
| DAG propagation | each gene inherits all ancestors of its direct terms over `is_a` + `part_of` |
| Universe | the 17,444 genes DESeq2 actually tested (non-NA padj) **and** that carry a GO annotation — not all genes in the genome |
| Term size | 10–500 genes in the universe; ≥ 5 DE genes to be reported |
| Test | one-sided hypergeometric, Benjamini–Hochberg within each ontology and direction |
| Confound control | expression-matched permutation test (below), 5,000 permutations |
| Redundancy | a term is dropped if it shares > 50 % of its DE genes with an already-selected stronger term |
| Figure | top 8 surviving terms per ontology per direction |

Reference data: GO `releases/2026-07-26`, MGI GAF generated 2026-05-21.

### Why the expression-matched test is here

The DE calls in this dataset are strongly stratified by expression level:

| | median baseMean | median \|log2FC\| |
|---|---|---|
| Down in cKO (n = 3,138) | 35 | 2.64 |
| Up in cKO (n = 1,149) | 2,638 | 0.35 |

Genes detected at ≥ 10 normalized counts fall from ~19,900 (WT) to ~17,400 (cKO). So a
plain hypergeometric test can report a term simply because its genes are lowly (or
highly) expressed. `expression_matched_test.py` rebuilds the null by drawing background
gene sets that match the DE set's baseMean-decile distribution, and reports an empirical
p-value and FDR. The plotted terms must pass **both** tests at FDR < 0.05.

This matters: 310/390 "down" and 179/235 "up" terms survive, but the ones that fail are
not random. Every down-regulated *immune* term (immune effector process, adaptive immune
response, …) fails expression matching, while the up-regulated type I interferon,
complement and inflammasome terms survive — i.e. the real immune signal in this dataset
points **up**, and the plain test would have shown it in both directions.

## Files

| File | Contents |
|---|---|
| `go_enrichment.py` | OBO/GAF parsing, DAG propagation, hypergeometric test |
| `expression_matched_test.py` | expression-matched permutation null |
| `select_terms.py` | double-significance filter, redundancy pruning, top-N |
| `plot_go_bubble.py` | the bubble plot (PNG + PDF) |
| `run_all.sh` | runs the four steps in order, fetching reference data first |
| `../results/go_enrichment_exprmatched.tsv` | every tested term, both directions, with member genes and both sets of p/FDR values |
| `../results/go_terms_plot.tsv` | the terms drawn in the figure |

## Reading the figure

Each bubble is one GO term. x = fold enrichment (observed / expected fraction of DE
genes), bubble area ∝ number of DE genes in the term, colour = −log10 FDR from the
hypergeometric test. Rows are ontologies, columns are the two directions.

## Caveats

- This is bulk RNA-seq of a degenerating tissue with n = 3 per genotype. Over-representation
  here reflects **both** transcriptional regulation and shifts in cell composition, and the
  two cannot be separated from this data. The "up in cKO" photoreceptor-cilium and
  outer-segment signal is carried by rod genes (`Rp1`, `Pde6a`, `Cnga1`, `Sag`, `Gnat1`,
  `Ush2a`, `Bbs1/2/4/7/12`), whose *share* of the library rises while cone genes (`Pde6c`
  -1.20, `Arr3` -1.03, `Opn1sw` -0.83) fall — consistent with cones being lost
  disproportionately, not with rod genes being transcriptionally induced.
- The reduced gene detection in the cKO (~2,500 fewer genes at ≥ 10 counts) is itself
  either a phenotype or a library-complexity artifact; normalized counts alone cannot
  distinguish them. Raw counts, library sizes and mapping rates would settle it.
- IEA (electronically inferred) annotations are included, which is standard for mouse but
  makes broad terms broader.
