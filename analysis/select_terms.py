#!/usr/bin/env python3
"""
Pick the GO terms to plot: significant by both tests, optionally restricted to a
GO allow-list, then de-duplicated by gene overlap.

Usage: select_terms.py IN.tsv OUT.tsv "BP:8,CC:5" [allowlist.tsv]

Adds a `significant` column.  A direction in which nothing clears FDR < 0.05
would otherwise contribute an empty panel, which hides whether the terms were
tested at all; such a direction instead contributes its strongest nominal
candidates (raw p < 0.10 and fold enrichment > 1), flagged significant = False
so the figure can draw them as hollow bubbles.  They are shown to say how far
short the direction falls, and are not enrichment claims.

Redundancy metric
-----------------
Jaccard, |A and B| / |A or B|, not the overlap coefficient |A and B| / min(|A|,|B|).
GO terms nest, and a child's genes are a subset of its parent's, so the overlap
coefficient is 1.0 for every parent/child pair and discards every child of a
kept term. Jaccard only collapses terms whose gene sets are genuinely close in
size as well as content, which keeps a specific term beside the broad one it
sits under while still removing near-duplicates.
"""
import sys

import pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]
# how many terms to keep per ontology, e.g. "BP:8,CC:5" -- ontologies left out
# of the spec are dropped entirely
SPEC = sys.argv[3] if len(sys.argv) > 3 else "BP:8,CC:5"
TOP_N = {k: int(v) for k, v in (p.split(":") for p in SPEC.split(","))}
ALLOWLIST = sys.argv[4] if len(sys.argv) > 4 else None

FDR_CUT = 0.05
RELAXED_P = 0.10        # display-only fallback for a direction with no hits
JACCARD_CUT = 0.5

d = pd.read_csv(IN, sep="\t")
if ALLOWLIST:
    allow = set(pd.read_csv(ALLOWLIST, sep="\t").GO_ID)
    before = len(d)
    d = d[d.GO_ID.isin(allow)]
    print(f"allow-list {ALLOWLIST}: {before} -> {len(d)} terms tested",
          file=sys.stderr)
d = d.copy()
d["gene_set"] = d.genes.str.split(";").apply(set)


def pick(sub, n):
    """Strongest terms by p, dropping any that duplicate a stronger one."""
    chosen = []
    # GO_ID breaks ties so the selection is reproducible when p-values are
    # identical (synonymous terms such as GO:0015267 / GO:0022803 share a gene set)
    for _, row in sub.sort_values(["pvalue", "FDR", "GO_ID"]).iterrows():
        g = row.gene_set
        if any(len(g & c.gene_set) / len(g | c.gene_set) > JACCARD_CUT
               for c in chosen):
            continue
        chosen.append(row)
        if len(chosen) == n:
            break
    return chosen


kept = []
for direction, per_dir in d.groupby("direction"):
    hits = per_dir[(per_dir.FDR < FDR_CUT) & (per_dir.FDR_expr_matched < FDR_CUT)]
    significant = len(hits) > 0
    if significant:
        pool = hits
    else:
        pool = per_dir[(per_dir.pvalue < RELAXED_P) & (per_dir.fold_enrichment > 1)]
        print(f"{direction}: nothing at FDR < {FDR_CUT}; falling back to "
              f"{len(pool)} nominal candidates (p < {RELAXED_P})", file=sys.stderr)
    for ns, sub in pool.groupby("ontology"):
        if ns not in TOP_N:
            continue
        for row in pick(sub, TOP_N[ns]):
            row = row.copy()
            row["significant"] = significant
            kept.append(row)

if not kept:
    sys.exit("no terms survived the filters")

out = pd.DataFrame(kept).drop(columns=["gene_set"])
out.to_csv(OUT, sep="\t", index=False)
for (dirn, ns), s in out.groupby(["direction", "ontology"]):
    print(f"\n### {dirn} / {ns}  (significant={bool(s.significant.iloc[0])})")
    print(s[["GO_ID", "term", "n_DE_in_term", "fold_enrichment",
             "pvalue", "FDR", "FDR_expr_matched"]].to_string(index=False))
