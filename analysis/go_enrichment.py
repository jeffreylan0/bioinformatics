#!/usr/bin/env python3
"""
GO over-representation analysis for the 1-month CHX10-Cre; Ddx41 fl/fl retina RNA-seq.

Inputs : 1m-deseq.tabular  (DESeq2 table, no header)
         GO ontology (go-basic.obo) + mouse GO annotations (MGI GAF)
Method : annotations propagated over the GO DAG (is_a + part_of), one-sided
         hypergeometric test against the set of genes DESeq2 actually tested,
         Benjamini-Hochberg FDR within each ontology/direction.
Outputs: go_enrichment_<direction>.tsv  and a combined table.
"""
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

DESEQ = sys.argv[1]
OBO = sys.argv[2]
GAF = sys.argv[3]
OUTDIR = sys.argv[4]

PADJ_CUT = 0.05
MIN_TERM, MAX_TERM = 10, 500   # term size in the tested universe
MIN_HITS = 5                   # DE genes needed in a term to report it

# ---------------------------------------------------------------- ontology
def parse_obo(path):
    terms, alt = {}, {}
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {"id": None, "name": None, "ns": None,
                       "parents": set(), "alt": [], "obsolete": False}
                continue
            if line.startswith("["):          # [Typedef] etc.
                if cur and cur["id"]:
                    terms[cur["id"]] = cur
                cur = None
                continue
            if cur is None or not line:
                if cur and not line and cur["id"]:
                    terms[cur["id"]] = cur
                    cur = None
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                cur["id"] = val
            elif key == "name":
                cur["name"] = val
            elif key == "namespace":
                cur["ns"] = val
            elif key == "alt_id":
                cur["alt"].append(val)
            elif key == "is_a":
                cur["parents"].add(val.split(" ! ")[0].strip())
            elif key == "relationship":
                rel, _, target = val.partition(" ")
                if rel == "part_of":
                    cur["parents"].add(target.split(" ! ")[0].strip())
            elif key == "is_obsolete" and val == "true":
                cur["obsolete"] = True
    if cur and cur["id"]:
        terms[cur["id"]] = cur
    terms = {k: v for k, v in terms.items() if not v["obsolete"]}
    for t in terms.values():
        for a in t["alt"]:
            alt[a] = t["id"]
    return terms, alt


def ancestors(go_id, terms, cache):
    """All ancestors of go_id, inclusive, over is_a/part_of."""
    if go_id in cache:
        return cache[go_id]
    cache[go_id] = {go_id}            # guard against cycles
    out = {go_id}
    for p in terms[go_id]["parents"]:
        if p in terms:
            out |= ancestors(p, terms, cache)
    cache[go_id] = out
    return out


# -------------------------------------------------------------- annotation
def parse_gaf(path, terms, alt):
    """symbol -> set of directly annotated GO ids (mouse, non-NOT, non-ND)."""
    gene2go = defaultdict(set)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            symbol, qualifier, go_id, evidence, taxon = f[2], f[3], f[4], f[6], f[12]
            if "NOT" in qualifier or evidence == "ND":
                continue
            if not taxon.startswith("taxon:10090"):
                continue
            go_id = alt.get(go_id, go_id)
            if go_id in terms:
                gene2go[symbol].add(go_id)
    return gene2go


# -------------------------------------------------------------------- main
sys.setrecursionlimit(100000)
terms, alt = parse_obo(OBO)
print(f"ontology terms (non-obsolete): {len(terms)}", file=sys.stderr)
gene2go_direct = parse_gaf(GAF, terms, alt)
print(f"annotated mouse genes: {len(gene2go_direct)}", file=sys.stderr)

cache = {}
gene2go = {g: set().union(*(ancestors(t, terms, cache) for t in gos))
           for g, gos in gene2go_direct.items()}

cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
res = pd.read_csv(DESEQ, sep="\t", header=None, names=cols,
                  na_values=["NA", "", "NaN"])
tested = res.dropna(subset=["padj"]).copy()

# universe = genes DESeq2 tested AND that carry at least one GO annotation
universe = sorted(set(tested.gene) & set(gene2go))
uset = set(universe)
tested = tested[tested.gene.isin(uset)]
print(f"DESeq2-tested genes: {len(res.dropna(subset=['padj']))}; "
      f"with GO annotation (universe): {len(universe)}", file=sys.stderr)

sig = tested[tested.padj < PADJ_CUT]
sets = {
    "Down in cKO": sorted(sig.loc[sig.log2FoldChange < 0, "gene"]),
    "Up in cKO":   sorted(sig.loc[sig.log2FoldChange > 0, "gene"]),
}

# term -> genes in universe
term2genes = defaultdict(set)
for g in universe:
    for t in gene2go[g]:
        term2genes[t].add(g)

NS_SHORT = {"biological_process": "BP",
            "molecular_function": "MF",
            "cellular_component": "CC"}

M = len(universe)
all_out = []
for label, genes in sets.items():
    gset = set(genes)
    N = len(gset)
    rows = []
    for t, tg in term2genes.items():
        n = len(tg)
        if n < MIN_TERM or n > MAX_TERM:
            continue
        hits = tg & gset
        k = len(hits)
        if k < MIN_HITS:
            continue
        p = hypergeom.sf(k - 1, M, n, N)
        rows.append({
            "GO_ID": t,
            "term": terms[t]["name"],
            "ontology": NS_SHORT.get(terms[t]["ns"], terms[t]["ns"]),
            "direction": label,
            "n_DE_in_term": k,
            "n_term_in_universe": n,
            "n_DE_total": N,
            "n_universe": M,
            "gene_ratio": k / N,
            "bg_ratio": n / M,
            "fold_enrichment": (k / N) / (n / M),
            "pvalue": p,
            "genes": ";".join(sorted(hits)),
        })
    df = pd.DataFrame(rows)
    # Benjamini-Hochberg within each ontology for this direction
    df["FDR"] = np.nan
    for ns, idx in df.groupby("ontology").groups.items():
        sub = df.loc[idx].sort_values("pvalue")
        m = len(sub)
        q = sub.pvalue.values * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1].clip(max=1.0)
        df.loc[sub.index, "FDR"] = q
    df = df.sort_values(["ontology", "FDR", "pvalue"])
    tag = label.split()[0].lower()
    df.to_csv(f"{OUTDIR}/go_enrichment_{tag}.tsv", sep="\t", index=False)
    print(f"{label}: {N} genes, {len(df)} terms tested, "
          f"{(df.FDR < 0.05).sum()} with FDR<0.05", file=sys.stderr)
    all_out.append(df)

pd.concat(all_out).to_csv(f"{OUTDIR}/go_enrichment_all.tsv", sep="\t", index=False)

# gene x term binary matrix, for the expression-matched permutation test
from scipy import sparse
keep_terms = sorted(t for t, g in term2genes.items() if MIN_TERM <= len(g) <= MAX_TERM)
col = {t: j for j, t in enumerate(keep_terms)}
rows, cols_ = [], []
for i, g in enumerate(universe):
    for t in gene2go[g]:
        j = col.get(t)
        if j is not None:
            rows.append(i)
            cols_.append(j)
mat = sparse.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols_)),
                        shape=(len(universe), len(keep_terms)))
np.savez_compressed(f"{OUTDIR}/gene2go_matrix.npz", data=mat.data, indices=mat.indices,
                    indptr=mat.indptr, shape=np.array(mat.shape),
                    universe=np.array(universe, dtype=object),
                    terms=np.array(keep_terms, dtype=object))
