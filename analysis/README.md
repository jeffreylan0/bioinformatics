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
| Scope | the figure is restricted to retina/vision GO terms (see below); the enrichment itself is run genome-wide |
| Display fallback | a direction with no term at FDR < 0.05 shows nominal candidates (p < 0.10, fold enrichment > 1) as hollow bubbles |
| Redundancy | a term is dropped if its DE genes are > 50 % Jaccard-similar to an already-selected stronger term |
| Figure | top 8 biological-process and top 5 cellular-component terms per direction; molecular function is computed but not plotted |

Reference data: GO `releases/2026-07-26`, MGI GAF generated 2026-05-21.

### Defining "retina-related"

`retina_go_set.py` takes 22 curated entry points into the eye/vision part of the
ontology — eye and camera-type eye development, retina development and
homeostasis, photoreceptor differentiation/development/maintenance, visual
perception, phototransduction, opsin signalling, retinoid/retinal/retinol
metabolism, and the photoreceptor outer segment, inner segment and cilium — and
takes every descendant of them over `is_a` + `part_of`, giving 174 terms
(156 BP, 15 CC, 3 MF).

Membership therefore follows the ontology, not the wording of a term name. That
matters in both directions: `9+0 non-motile cilium` is excluded despite being
driven here entirely by photoreceptor genes, because it is not under an eye
root, while `photoreceptor connecting cilium` is included. `sensory organ
development` (GO:0007423) is deliberately *not* a root, since it also covers
ear, nose and lateral line.

**The result is one-sided.** Retina terms were tested in both directions (22
down, 29 up), but none reached significance among down-regulated genes — the
best is retinoid metabolic process at raw p = 0.017, FDR = 0.11. All 19
significant retina terms are up in the cKO.

Both directions are still drawn, because an absent panel cannot distinguish
"tested and found nothing" from "never tested". A direction with no term at
FDR < 0.05 instead shows its strongest nominal candidates (raw p < 0.10 and
fold enrichment > 1) as **hollow grey bubbles**, which are not enrichment
claims. Colour is normalised over the significant terms only, so a pale bubble
never has to double as "not significant". The down panel's cellular-component
band is empty and labelled as such: its three retina CC terms are *depleted*
rather than enriched (fold enrichment 0.73–0.95), so none clears even the
nominal bar.

### Why the redundancy metric is Jaccard, not the overlap coefficient

GO terms nest, and a child's gene set is a subset of its parent's, so the
overlap coefficient |A ∩ B| / min(|A|,|B|) is 1.0 for every parent/child pair
and discards every child of a kept term. On the retina set that collapsed 14
significant BP terms to 3, throwing away exactly the specific terms the figure
is for. Jaccard |A ∩ B| / |A ∪ B| only collapses terms close in size as well as
content, keeping `photoreceptor cell maintenance` beside `visual system
development` while still removing near-duplicates.

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

This matters. Of the terms the hypergeometric test calls significant, 169/390 "down"
and 96/235 "up" also clear the expression-matched test — so matching removes well over
half of the down-regulated hits. (The expression-matched test taken on its own flags 310
down and 179 up terms across everything tested, but those are not the same terms, and
the figure requires both.)

The terms that fail are not random. Every down-regulated *immune* term (immune effector
process, adaptive immune response, …) fails expression matching, while the up-regulated
type I interferon, complement and inflammasome terms survive — i.e. the real immune
signal in this dataset points **up**, and the plain test would have shown it in both
directions.

## Files

| File | Contents |
|---|---|
| `go_enrichment.py` | OBO/GAF parsing, DAG propagation, hypergeometric test |
| `expression_matched_test.py` | expression-matched permutation null |
| `retina_go_set.py` | the retina/vision GO subgraph from curated roots |
| `select_terms.py` | double-significance filter, optional GO allow-list, redundancy pruning, top-N |
| `plot_go_bubble.py` | the bubble plot (PNG + PDF) |
| `run_all.sh` | runs the steps in order, fetching reference data first |
`results/` ships the figure only — `go_bubble_plot.svg` and `go_bubble_plot.pdf`.
The intermediate tables (`go_enrichment_exprmatched.tsv` with every tested term,
its member genes and both sets of p/FDR values; `go_terms_plot.tsv` with the
terms drawn) are regenerated by `run_all.sh` and are not committed. Pass
`GO_BUBBLE_PNG=1` to also write a raster proof.

## Reading the figure

Each bubble is one GO term, labelled with its name and GO accession. x = fold
enrichment (observed / expected fraction of DE genes), bubble area ∝ number of DE
genes in the term, colour = −log10 FDR from the hypergeometric test. The two
columns are the enriched directions (here only one); within each column,
biological-process and
cellular-component terms share one plot and the ontology is carried on the
y-axis as a bracketed group label. The panel count follows the data.

Output is vector only: SVG with live text (`svg.fonttype = "none"`, Liberation
Sans with an Arial/Helvetica fallback stack, no mathtext) so an editor can
restyle the labels, and PDF with subset-embedded TrueType fonts.

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
