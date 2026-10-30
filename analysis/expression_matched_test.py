#!/usr/bin/env python3
"""
Expression-matched sensitivity test for the GO over-representation results.

In this dataset the DE direction is strongly confounded with expression level
(down-regulated genes are mostly low-expressed, up-regulated genes mostly
high-expressed), so a plain hypergeometric universe can report terms that are
merely enriched for lowly- or highly-expressed genes.  This test rebuilds the
null by sampling background gene sets that match the DE set's baseMean decile
distribution, and reports an empirical p-value and BH-FDR per GO term.
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import sparse

ENRICH_TSV, DESEQ, GENE2GO_NPZ, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
B = int(sys.argv[5]) if len(sys.argv) > 5 else 5000
RNG = np.random.default_rng(0)
N_DECILES = 10

store = np.load(GENE2GO_NPZ, allow_pickle=True)
universe = list(store["universe"])
terms = list(store["terms"])
Mmat = sparse.csr_matrix((store["data"], store["indices"], store["indptr"]),
                         shape=tuple(store["shape"]))      # genes x terms, binary
gidx = {g: i for i, g in enumerate(universe)}
tidx = {t: i for i, t in enumerate(terms)}

cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
res = pd.read_csv(DESEQ, sep="\t", header=None, names=cols,
                  na_values=["NA", "", "NaN"]).dropna(subset=["padj"])
res = res[res.gene.isin(gidx)].set_index("gene")
bm = res.loc[universe, "baseMean"].values

# decile strata on log expression
ranks = pd.qcut(pd.Series(bm).rank(method="first"), N_DECILES, labels=False).values
strata = [np.where(ranks == d)[0] for d in range(N_DECILES)]

enr = pd.read_csv(ENRICH_TSV, sep="\t")
out = []
for direction, sub in enr.groupby("direction"):
    sig_genes = set()
    for s in sub.genes:
        sig_genes.update(str(s).split(";"))
    de_idx = np.array(sorted(gidx[g] for g in sig_genes if g in gidx))
    per_stratum = np.array([np.isin(de_idx, s).sum() for s in strata])

    keep = [tidx[t] for t in sub.GO_ID]
    sub_mat = Mmat[:, keep]                                   # genes x terms
    obs = np.asarray(sub_mat[de_idx].sum(axis=0)).ravel()

    ge = np.zeros(len(keep), dtype=np.int64)
    for b in range(B):
        pick = np.concatenate([RNG.choice(s, size=k, replace=False)
                               for s, k in zip(strata, per_stratum) if k])
        null = np.asarray(sub_mat[pick].sum(axis=0)).ravel()
        ge += (null >= obs)
    p = (ge + 1) / (B + 1)

    s2 = sub.copy()
    s2["p_expr_matched"] = p
    s2["FDR_expr_matched"] = np.nan
    for ns, idx in s2.groupby("ontology").groups.items():
        t = s2.loc[idx].sort_values("p_expr_matched")
        m = len(t)
        q = t.p_expr_matched.values * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1].clip(max=1.0)
        s2.loc[t.index, "FDR_expr_matched"] = q
    out.append(s2)
    print(f"{direction}: {len(s2)} terms, "
          f"{(s2.FDR_expr_matched < 0.05).sum()} keep FDR<0.05 after expression matching",
          file=sys.stderr)

pd.concat(out).to_csv(OUT, sep="\t", index=False)
