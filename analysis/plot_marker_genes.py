#!/usr/bin/env python3
"""
Production figure: every marker gene behind the cell-population summary, with
its own log2 fold change.

One small panel per population, genes sorted by fold change. The stem runs from
zero, so the bar length is the fold change itself; the dashed line is the shift
a highly expressed gene gets from library scaling alone, which is the bar a gene
has to clear before "up" means anything about composition. Hollow dots are genes
that are not significant on their own (padj >= 0.05).

Usage: plot_marker_genes.py deseq.tabular OUT_BASE
Writes vector only (SVG + PDF); set FIGURE_PNG=1 for a raster proof.
"""
import os
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import figure_style as fs
from marker_panels import PANELS

DESEQ, OUT_BASE = sys.argv[1], sys.argv[2]
HI_EXPR = 5000
PADJ_CUT = 0.05

fs.apply_rcparams()

cols = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
res = (pd.read_csv(DESEQ, sep="\t", header=None, names=cols,
                   na_values=["NA", "", "NaN"])
       .set_index("gene").dropna(subset=["log2FoldChange"]))
ref = res[res.baseMean >= HI_EXPR].log2FoldChange.median()

# panels ordered by median fold change, strongest gain first
order = sorted(PANELS,
               key=lambda k: -res.loc[[g for g in PANELS[k] if g in res.index],
                                      "log2FoldChange"].median())
data = {k: res.loc[[g for g in PANELS[k] if g in res.index]]
             .sort_values("log2FoldChange", ascending=False) for k in order}

# ---------------------------------------------------------------- geometry ---
N_COL = 3
ROW_IN = 0.190                       # one gene
TITLE_IN = 0.30                      # panel heading
GAP_IN = 0.34                        # between panels
HEADER_IN, FOOTER_IN = 1.60, 0.85
LABEL_IN, PLOT_IN, COLGAP_IN = 1.05, 2.90, 0.65
MARGIN_IN = 0.50


def column_height(cats):
    return (sum(len(data[k]) * ROW_IN + TITLE_IN for k in cats)
            + GAP_IN * max(len(cats) - 1, 0))


# Split the ordered panels into contiguous columns, choosing the split points
# that minimise the tallest column. Contiguous keeps the median-fold-change
# ordering readable down one column and on to the next; one panel (rod, 22
# genes) is much larger than the rest, so a naive halfway split is lopsided.
best, best_h = None, float("inf")
for cuts in combinations(range(1, len(order)), N_COL - 1):
    bounds = (0,) + cuts + (len(order),)
    cols = [order[a:b] for a, b in zip(bounds[:-1], bounds[1:])]
    if any(not c for c in cols):
        continue
    h = max(column_height(c) for c in cols)
    if h < best_h:
        best, best_h = cols, h
columns = best

col_h = [column_height(c) for c in columns]
FIG_H = max(col_h) + HEADER_IN + FOOTER_IN
FIG_W = (2 * MARGIN_IN + N_COL * (LABEL_IN + PLOT_IN)
         + (N_COL - 1) * COLGAP_IN)

lo = min(d.log2FoldChange.min() for d in data.values())
hi = max(d.log2FoldChange.max() for d in data.values())
pad = 0.08 * (hi - lo)
XLIM = (lo - pad, hi + pad)

fig = plt.figure(figsize=(FIG_W, FIG_H))
bottom_axes = []

for ci, cat_list in enumerate(columns):
    x0 = (MARGIN_IN + ci * (LABEL_IN + PLOT_IN + COLGAP_IN) + LABEL_IN) / FIG_W
    w = PLOT_IN / FIG_W
    y_cursor = FIG_H - HEADER_IN       # inches from the bottom, walking down

    for ki, cat in enumerate(cat_list):
        sub = data[cat]
        h = len(sub) * ROW_IN
        y_cursor -= TITLE_IN
        ax = fig.add_axes([x0, (y_cursor - h) / FIG_H, w, h / FIG_H])
        y = np.arange(len(sub))[::-1]          # first gene at the top

        ax.axvline(0, color=fs.BASELINE, lw=1.0, zorder=1)
        ax.axvline(ref, color=fs.MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)

        for yi, (_, r) in zip(y, sub.iterrows()):
            colour = fs.POS if r.log2FoldChange > 0 else fs.NEG
            sig = (not np.isnan(r.padj)) and r.padj < PADJ_CUT
            ax.plot([0, r.log2FoldChange], [yi, yi], color=colour, lw=1.8,
                    solid_capstyle="butt", alpha=0.85, zorder=2)
            ax.scatter([r.log2FoldChange], [yi], s=34,
                       facecolor=colour if sig else fs.SURFACE,
                       edgecolor=fs.SURFACE if sig else colour,
                       linewidth=1.0 if sig else 1.3, zorder=3)

        ax.set_yticks(y)
        ax.set_yticklabels(sub.index, fontsize=7.6, color=fs.INK,
                           fontstyle="italic")
        ax.set_ylim(-0.75, len(sub) - 0.25)
        ax.set_xlim(*XLIM)
        fs.tidy(ax)
        ax.set_title(f"{cat}   ({len(sub)} genes)", fontsize=9.6, color=fs.INK,
                     fontweight="bold", loc="left", pad=5)

        last = ki == len(cat_list) - 1
        if last:
            bottom_axes.append(ax)
            ax.set_xlabel("log2 fold change (cKO / control)", fontsize=9,
                          color=fs.INK_2, labelpad=6)
        else:
            ax.set_xticklabels([])
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", length=0)
        y_cursor -= h + GAP_IN

# ------------------------------------------------------------------ legend ---
handles = [
    Line2D([], [], marker="o", linestyle="none", markersize=7,
           markerfacecolor=fs.MUTED, markeredgecolor=fs.SURFACE,
           markeredgewidth=1.0, label=f"padj < {PADJ_CUT}"),
    Line2D([], [], marker="o", linestyle="none", markersize=7,
           markerfacecolor=fs.SURFACE, markeredgecolor=fs.MUTED,
           markeredgewidth=1.3, label="not significant"),
    Line2D([], [], color=fs.MUTED, lw=0.9, ls=(0, (4, 3)),
           label=f"library-wide shift ({ref:+.2f})"),
]
leg = fig.legend(handles=handles, loc="lower center",
                 bbox_to_anchor=(0.5, 0.004), ncol=3, frameon=False,
                 fontsize=8.8, labelcolor=fs.INK_2, handletextpad=0.6,
                 columnspacing=2.8)

TX = MARGIN_IN / FIG_W
fig.suptitle("Marker genes behind each cell population, 1-month Ddx41 "
             "conditional-knockout retina",
             x=TX, y=1 - 0.40 / FIG_H, ha="left", fontsize=14, color=fs.INK,
             fontweight="bold")
fig.text(TX, 1 - 0.70 / FIG_H,
         "CHX10-Cre; Ddx41 fl/fl vs. littermate control, n = 3 per genotype. "
         "Every gene in the curated panels, sorted by fold change within its "
         "population.\nCones, horizontal cells and retinal ganglion cells fall "
         "coherently; rod genes sit just above the library-scaling line, which "
         "is what a gain in share\nlooks like rather than induction. Reactive "
         "gliosis is the one population whose genes move by more than two-fold "
         "on their own.",
         ha="left", va="top", fontsize=8.8, color=fs.MUTED, linespacing=1.6)

fig.savefig(f"{OUT_BASE}.svg")
fig.savefig(f"{OUT_BASE}.pdf")
out = [f"{OUT_BASE}.svg", f"{OUT_BASE}.pdf"]
if os.environ.get("FIGURE_PNG"):
    fig.savefig(f"{OUT_BASE}.png", dpi=300)
    out.append(f"{OUT_BASE}.png")
print("wrote " + ", ".join(out))
