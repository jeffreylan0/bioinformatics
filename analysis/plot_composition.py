#!/usr/bin/env python3
"""
Production figure: the 1-month Ddx41 cKO retina loses cones and inner neurons
while rods and glia gain share, so photoreceptor gene sets rise in *relative*
terms without being transcriptionally induced.

Shows the log2 fold change of curated cell-type marker panels, one dot per gene
plus the panel median, against the shift a highly expressed gene gets from
library scaling alone. Gene-level detail is a separate figure
(plot_marker_genes.py).

Usage: plot_composition.py counts.tabular deseq.tabular OUT_BASE
Writes vector only (SVG + PDF); set FIGURE_PNG=1 for a raster proof.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import figure_style as fs
from marker_panels import PANELS

COUNTS, DESEQ, OUT_BASE = sys.argv[1], sys.argv[2], sys.argv[3]
HI_EXPR = 5000        # "highly expressed" cut for the library-wide reference
JITTER = 0.13
RNG = np.random.default_rng(0)


fs.apply_rcparams()

cnt = pd.read_csv(COUNTS, sep="\t", index_col=0)
cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
res = (pd.read_csv(DESEQ, sep="\t", header=None, names=cols,
                   na_values=["NA", "", "NaN"])
       .set_index("gene").dropna(subset=["log2FoldChange"]))

wt = [c for c in cnt.columns if "WT" in c]
ko = [c for c in cnt.columns if "-DDX41-" in c and "WT" not in c]
tot_wt, tot_ko = cnt[wt].sum().sum(), cnt[ko].sum().sum()

# Reference: where the bulk of highly expressed genes sits. This is the shift a
# gene gets for free from library scaling, and the bar a panel must clear to be
# called a real gain.
ref = res[res.baseMean >= HI_EXPR].log2FoldChange.median()

rows = []
for name, genes in PANELS.items():
    g = [x for x in genes if x in res.index and x in cnt.index]
    lfc = res.loc[g, "log2FoldChange"]
    rows.append({
        "panel": name, "n": len(g), "genes": g, "lfc": lfc.values,
        "median": lfc.median(),
        "wt_share": cnt.loc[g, wt].sum().sum() / tot_wt,
        "ko_share": cnt.loc[g, ko].sum().sum() / tot_ko,
    })
t = pd.DataFrame(rows).sort_values("median").reset_index(drop=True)  # gains on top
t["pct"] = t.ko_share / t.wt_share - 1
y = np.arange(len(t))

# ---------------------------------------------------------------- geometry ---
ROW_IN = 0.56
HEADER_IN, FOOTER_IN = 1.75, 0.95
LEFT_IN, RIGHT_IN = 2.35, 0.55
PANEL_W = 6.60
FIG_H = ROW_IN * len(t) + HEADER_IN + FOOTER_IN
FIG_W = LEFT_IN + PANEL_W + RIGHT_IN

fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(1, 1,
                      left=LEFT_IN / FIG_W, right=(FIG_W - RIGHT_IN) / FIG_W,
                      top=1 - HEADER_IN / FIG_H, bottom=FOOTER_IN / FIG_H)

axA = fig.add_subplot(gs[0, 0])
axA.axvline(0, color=fs.BASELINE, lw=1.0, zorder=1)
axA.axvline(ref, color=fs.MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)

for i, r in t.iterrows():
    colour = fs.POS if r["median"] > 0 else fs.NEG
    axA.scatter(r["lfc"], np.full(len(r["lfc"]), i)
                + RNG.uniform(-JITTER, JITTER, len(r["lfc"])),
                s=26, facecolor=colour, alpha=0.30, edgecolor="none", zorder=2)
    axA.scatter([r["median"]], [i], s=150, marker="D", facecolor=colour,
                edgecolor=fs.SURFACE, linewidth=1.1, zorder=4)

axA.set_yticks(y)
axA.set_yticklabels([])
axA.set_ylim(-0.7, len(t) - 0.3)
axA.set_xlim(-1.95, 4.75)
fs.tidy(axA)
axA.set_xlabel("log2 fold change (cKO / control)", fontsize=9.5,
               color=fs.INK_2, labelpad=7)
axA.set_title("Marker gene fold change by cell population", fontsize=11,
              color=fs.INK, fontweight="bold", pad=10)

# row labels: panel name with its gene count underneath
trans = axA.get_yaxis_transform()
for i, r in t.iterrows():
    axA.annotate(r["panel"], xy=(-0.018, i), xycoords=trans, xytext=(0, 1.5),
                 textcoords="offset points", ha="right", va="bottom",
                 fontsize=9.5, color=fs.INK)
    axA.annotate(f"{r['n']} genes", xy=(-0.018, i), xycoords=trans,
                 xytext=(0, -2.5), textcoords="offset points", ha="right",
                 va="top", fontsize=7.4, color=fs.MUTED)

axA.annotate(f"library-wide shift of highly\nexpressed genes ({ref:+.2f})",
             xy=(ref, 1.35), xytext=(0.80, 1.55),
             fontsize=7.8, color=fs.MUTED, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.8,
                             shrinkA=0, shrinkB=2))

# the decomposition, stated where it applies
rod = t.index[t.panel == "Rod"][0]
axA.annotate(
    f"rod {t.loc[rod, 'median']:+.2f}  =  {ref:+.2f} library-wide scaling\n"
    f"  +  {t.loc[rod, 'median'] - ref:+.2f} genuine gain in rod share",
    xy=(t.loc[rod, "median"], rod), xytext=(1.35, rod + 0.60),
    fontsize=8.2, color=fs.INK_2, ha="left", va="center",
    arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.8, shrinkA=0,
                    shrinkB=9))

# ------------------------------------------------------------- titles ---------
TX = 0.50 / FIG_W
fig.suptitle("Cell-type composition shifts in the 1-month Ddx41 "
             "conditional-knockout retina",
             x=TX, y=1 - 0.42 / FIG_H, ha="left", fontsize=14,
             color=fs.INK, fontweight="bold")
fig.text(TX, 1 - 0.72 / FIG_H,
         "CHX10-Cre; Ddx41 fl/fl vs. littermate control, n = 3 per genotype. "
         "Each dot is one marker gene; the diamond is the panel median. The "
         "dashed line is the median\nshift of all genes with baseMean >= 5,000, "
         "i.e. what a highly expressed gene gains from library scaling alone. "
         "Rods clear that line by only 0.15 log2 units\n(about 11%), so the "
         "photoreceptor gene sets that score as enriched are rising in share as "
         "cones and inner neurons are lost, not being induced.",
         ha="left", va="top", fontsize=8.8, color=fs.MUTED, linespacing=1.6)

fig.savefig(f"{OUT_BASE}.svg")
fig.savefig(f"{OUT_BASE}.pdf")
out = [f"{OUT_BASE}.svg", f"{OUT_BASE}.pdf"]
if os.environ.get("FIGURE_PNG"):
    fig.savefig(f"{OUT_BASE}.png", dpi=300)
    out.append(f"{OUT_BASE}.png")
print("wrote " + ", ".join(out))
