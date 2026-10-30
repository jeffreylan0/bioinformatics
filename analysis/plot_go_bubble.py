#!/usr/bin/env python3
"""Bubble plot of the significantly enriched GO terms."""
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

TSV, OUT_BASE = sys.argv[1], sys.argv[2]

# --- design tokens (light surface) -------------------------------------------
SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK_2     = "#52514e"
MUTED     = "#898781"
GRID      = "#e1e0d9"
BASELINE  = "#c3c2b7"
# sequential blue ramp, steps 250-700 (step 250 is the lightest that clears 2:1
# on the light surface, so the palest bubble still reads as a mark)
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue", RAMP)

AREA_PER_GENE = 5.5     # marker area in pt^2 per DE gene -> area-proportional
WRAP = 46

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

d = pd.read_csv(TSV, sep="\t")
d["neglog10FDR"] = -np.log10(d.FDR)

DIRS = ["Down in cKO", "Up in cKO"]
NS_ORDER = ["BP", "CC", "MF"]
NS_FULL = {"BP": "Biological process", "CC": "Cellular component",
           "MF": "Molecular function"}

# constant row pitch everywhere: each grid row is sized by its fullest panel
row_n = {ns: max(len(d[(d.direction == x) & (d.ontology == ns)]) for x in DIRS)
         for ns in NS_ORDER}

norm = Normalize(vmin=d.neglog10FDR.min(), vmax=d.neglog10FDR.max())
xmax = d.fold_enrichment.max() * 1.18

FIG_H = 0.46 * sum(row_n.values()) + 2.6
HEADER_IN, FOOTER_IN = 1.85, 1.05      # inches reserved for titles / legend
fig = plt.figure(figsize=(14.5, FIG_H))
gs = fig.add_gridspec(
    len(NS_ORDER), 2,
    height_ratios=[row_n[ns] for ns in NS_ORDER],
    left=0.255, right=0.9, top=1 - HEADER_IN / FIG_H,
    bottom=FOOTER_IN / FIG_H, hspace=0.20, wspace=0.66,
)

for r, ns in enumerate(NS_ORDER):
    for c, direction in enumerate(DIRS):
        ax = fig.add_subplot(gs[r, c])
        sub = (d[(d.direction == direction) & (d.ontology == ns)]
               .sort_values("fold_enrichment"))
        y = np.arange(len(sub)) + (row_n[ns] - len(sub))   # top-align
        ax.scatter(sub.fold_enrichment, y,
                   s=sub.n_DE_in_term * AREA_PER_GENE,
                   c=sub.neglog10FDR, cmap=CMAP, norm=norm,
                   linewidths=0.9, edgecolors=SURFACE, zorder=3)

        ax.set_yticks(y)
        ax.set_yticklabels([textwrap.fill(t, WRAP) for t in sub.term],
                           fontsize=8.2, color=INK_2)
        ax.set_ylim(-0.7, row_n[ns] - 0.3)
        ax.set_xlim(1.0, xmax)
        ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(BASELINE)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.tick_params(axis="y", length=0, pad=4)
        ax.tick_params(axis="x", colors=MUTED, labelsize=8.2, length=3,
                       color=BASELINE)

        if r == 0:
            ax.set_title(direction, fontsize=11.5, color=INK,
                         fontweight="bold", pad=10)
        if r == len(NS_ORDER) - 1:
            ax.set_xlabel("Fold enrichment", fontsize=9.5, color=INK_2,
                          labelpad=7)
        else:
            ax.set_xticklabels([])
        if c == 0:
            ax.annotate(NS_FULL[ns], xy=(0, 1), xycoords="axes fraction",
                        xytext=(-168, 16), textcoords="offset points",
                        fontsize=10, color=INK, fontweight="bold",
                        ha="left", va="bottom")

# --- colour bar ---------------------------------------------------------------
cax = fig.add_axes([0.922, 0.40, 0.013, 0.26])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax)
cb.set_label("$-\\log_{10}$ FDR", fontsize=9, color=INK_2, labelpad=9)
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8.2, colors=MUTED, length=3, color=BASELINE)

# --- size legend --------------------------------------------------------------
ticks = [10, 30, 60, 100]
handles = [Line2D([], [], marker="o", linestyle="none",
                  markersize=np.sqrt(t * AREA_PER_GENE),
                  markerfacecolor=MUTED, markeredgecolor=SURFACE,
                  markeredgewidth=0.9, label=str(t)) for t in ticks]
leg = fig.legend(handles=handles, title="DE genes in term",
                 loc="lower center", bbox_to_anchor=(0.56, 0.004),
                 ncol=len(ticks), frameon=False, handletextpad=0.6,
                 columnspacing=2.6, labelspacing=0.9, borderpad=0.8,
                 fontsize=9, labelcolor=INK_2)
leg.get_title().set_fontsize(9)
leg.get_title().set_color(INK_2)

fig.suptitle("GO enrichment in the 1-month CHX10-Cre; Ddx41$^{fl/fl}$ retina",
             x=0.075, y=1 - 0.42 / FIG_H, ha="left", fontsize=15, color=INK,
             fontweight="bold")
fig.text(0.075, 1 - 0.72 / FIG_H,
         "Over-representation among genes differentially expressed vs. littermate control "
         "(DESeq2 padj < 0.05; n = 3 per genotype).\n"
         "Terms shown pass FDR < 0.05 in both a hypergeometric test against all "
         "17,444 tested GO-annotated genes and an expression-matched\npermutation test, "
         "and are de-duplicated by gene overlap; top 8 per ontology per direction.",
         ha="left", va="top", fontsize=8.8, color=MUTED, linespacing=1.55)

fig.savefig(f"{OUT_BASE}.png", dpi=300)
fig.savefig(f"{OUT_BASE}.pdf")
print(f"wrote {OUT_BASE}.png / .pdf")
