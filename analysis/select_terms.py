#!/usr/bin/env python3
"""Pick the GO terms to plot: significant by both tests, de-duplicated by gene overlap."""
import sys
import pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]
TOP_N = int(sys.argv[3]) if len(sys.argv) > 3 else 8
FDR_CUT = 0.05
OVERLAP_CUT = 0.5     # drop a term sharing this much of its DE genes with a stronger one
# GO_ID breaks ties so the selection is reproducible when p-values are identical
# (synonymous terms such as GO:0015267 / GO:0022803 carry the same gene set)

d = pd.read_csv(IN, sep="\t")
d = d[(d.FDR < FDR_CUT) & (d.FDR_expr_matched < FDR_CUT)].copy()
d["gene_set"] = d.genes.str.split(";").apply(set)

kept = []
for (direction, ns), sub in d.groupby(["direction", "ontology"]):
    chosen = []
    for _, row in sub.sort_values(["pvalue", "FDR", "GO_ID"]).iterrows():
        g = row.gene_set
        if any(len(g & c.gene_set) / min(len(g), len(c.gene_set)) > OVERLAP_CUT
               for c in chosen):
            continue
        chosen.append(row)
        if len(chosen) == TOP_N:
            break
    kept.extend(chosen)

out = pd.DataFrame(kept).drop(columns=["gene_set"])
out.to_csv(OUT, sep="\t", index=False)
for (dirn, ns), s in out.groupby(["direction", "ontology"]):
    print(f"\n### {dirn} / {ns}")
    print(s[["term", "n_DE_in_term", "n_term_in_universe", "fold_enrichment",
             "FDR", "FDR_expr_matched"]].to_string(index=False))
