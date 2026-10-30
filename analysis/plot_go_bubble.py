#!/usr/bin/env python3
"""
Production figure: GO over-representation in the 1-month Ddx41 cKO retina.

Writes vector output only (SVG with live, editable text + PDF with embedded
TrueType fonts).  Set GO_BUBBLE_PNG=1 to also drop a raster copy for proofing.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

TSV, OUT_BASE = sys.argv[1], sys.argv[2]

# --- design tokens (light surface) -------------------------------------------
SURFACE  = "#fcfcfb"
INK      = "#0b0b0b"
INK_2    = "#52514e"
MUTED    = "#898781"
GRID     = "#e1e0d9"
BASELINE = "#c3c2b7"
# sequential blue ramp, steps 250-700: step 250 is the lightest that still
# clears 2:1 on the light surface, so the palest bubble reads as a mark
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("seq_blue", RAMP)

AREA_PER_GENE = 5.5          # marker area in pt^2 per DE gene -> area-proportional

# Live text in the SVG and embedded TrueType in the PDF, so the figure stays
# editable in Illustrator/Inkscape.  Liberation Sans is metric-compatible with
# Arial; no mathtext anywhere, since it would not survive svg.fonttype="none".
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

d = pd.read_csv(TSV, sep="\t")
d["neglog10FDR"] = -np.log10(d.FDR)

DIRS = ["Down in cKO", "Up in cKO"]
NS_ORDER = [ns for ns in ("BP", "CC") if ns in set(d.ontology)]
NS_FULL = {"BP": "Biological process", "CC": "Cellular component"}

# constant row pitch everywhere: each grid row is sized by its fullest panel
row_n = {ns: max(len(d[(d.direction == x) & (d.ontology == ns)]) for x in DIRS)
         for ns in NS_ORDER}

norm = Normalize(vmin=d.neglog10FDR.min(), vmax=d.neglog10FDR.max())
xmax = d.fold_enrichment.max() * 1.18

ROW_IN = 0.52                      # inches per term row (two-line labels)
HEADER_IN, FOOTER_IN = 1.80, 1.00
FIG_H = ROW_IN * sum(row_n.values()) + HEADER_IN + FOOTER_IN
FIG_W = 13.6

fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(
    len(NS_ORDER), 2,
    height_ratios=[row_n[ns] for ns in NS_ORDER],
    left=0.235, right=0.90, top=1 - HEADER_IN / FIG_H,
    bottom=FOOTER_IN / FIG_H, hspace=0.22, wspace=0.60,
)

for r, ns in enumerate(NS_ORDER):
    for c, direction in enumerate(DIRS):
        ax = fig.add_subplot(gs[r, c])
        sub = (d[(d.direction == direction) & (d.ontology == ns)]
               .sort_values("fold_enrichment"))
        y = np.arange(len(sub)) + (row_n[ns] - len(sub))      # top-align

        ax.scatter(sub.fold_enrichment, y,
                   s=sub.n_DE_in_term * AREA_PER_GENE,
                   c=sub.neglog10FDR, cmap=CMAP, norm=norm,
                   linewidths=0.9, edgecolors=SURFACE, zorder=3)

        # term name with its GO accession set underneath, both right-aligned
        trans = ax.get_yaxis_transform()
        for yi, term, go_id in zip(y, sub.term, sub.GO_ID):
            ax.annotate(term, xy=(-0.018, yi), xycoords=trans,
                        xytext=(0, 1.5), textcoords="offset points",
                        ha="right", va="bottom", fontsize=9, color=INK)
            ax.annotate(go_id, xy=(-0.018, yi), xycoords=trans,
                        xytext=(0, -2.5), textcoords="offset points",
                        ha="right", va="top", fontsize=7.2, color=MUTED)

        ax.set_yticks(y)
        ax.set_yticklabels([])
        ax.set_ylim(-0.6, row_n[ns] - 0.4)
        ax.set_xlim(1.0, xmax)
        ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(BASELINE)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=3,
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
                        xytext=(-158, 15), textcoords="offset points",
                        fontsize=10.5, color=INK, fontweight="bold",
                        ha="left", va="bottom")

# --- colour bar ---------------------------------------------------------------
# Drawn by hand rather than with fig.colorbar: matplotlib's colorbar is a quad
# mesh that is rasterised by default, and un-rasterising it leaves hairline
# seams between the quads in PDF/SVG. Overlapping rectangles cannot seam, and
# keep the figure fully vector.
cax = fig.add_axes([0.922, 0.38, 0.0125, 0.26])
n_steps = 256
edges = np.linspace(norm.vmin, norm.vmax, n_steps + 1)
for lo, hi in zip(edges[:-1], edges[1:]):
    h = hi - lo
    cax.add_patch(Rectangle((0, lo - 0.25 * h), 1, 1.5 * h,
                            facecolor=CMAP(norm(0.5 * (lo + hi))),
                            edgecolor="none", linewidth=0))
cax.set_xlim(0, 1)
cax.set_ylim(norm.vmin, norm.vmax)
cax.set_xticks([])
cax.yaxis.tick_right()
cax.yaxis.set_label_position("right")
cax.set_ylabel("-log10 FDR", fontsize=9, color=INK_2, labelpad=9)
for sp in cax.spines.values():
    sp.set_visible(False)
cax.tick_params(axis="y", labelsize=8.5, colors=MUTED, length=3, color=BASELINE)

# --- size legend --------------------------------------------------------------
ticks = [10, 30, 60, 100]
handles = [Line2D([], [], marker="o", linestyle="none",
                  markersize=np.sqrt(t * AREA_PER_GENE),
                  markerfacecolor=MUTED, markeredgecolor=SURFACE,
                  markeredgewidth=0.9, label=str(t)) for t in ticks]
leg = fig.legend(handles=handles, title="DE genes in term",
                 loc="lower center", bbox_to_anchor=(0.55, 0.004),
                 ncol=len(ticks), frameon=False, handletextpad=0.6,
                 columnspacing=2.6, labelspacing=0.9, borderpad=0.8,
                 fontsize=9, labelcolor=INK_2)
leg.get_title().set_fontsize(9)
leg.get_title().set_color(INK_2)

fig.suptitle("GO enrichment in the 1-month Ddx41 conditional-knockout retina",
             x=0.062, y=1 - 0.40 / FIG_H, ha="left", fontsize=15, color=INK,
             fontweight="bold")
fig.text(0.062, 1 - 0.70 / FIG_H,
         "CHX10-Cre; Ddx41 fl/fl vs. littermate control, n = 3 per genotype. "
         "Over-representation among DESeq2 padj < 0.05 genes, against the 17,444 "
         "tested GO-annotated genes.\nTerms pass FDR < 0.05 in both a hypergeometric "
         "and an expression-matched permutation test, are de-duplicated by gene "
         "overlap, and are the top 8 (BP) and 5 (CC) per direction.",
         ha="left", va="top", fontsize=8.8, color=MUTED, linespacing=1.6)

fig.savefig(f"{OUT_BASE}.svg")
fig.savefig(f"{OUT_BASE}.pdf")
written = [f"{OUT_BASE}.svg", f"{OUT_BASE}.pdf"]
if os.environ.get("GO_BUBBLE_PNG"):
    fig.savefig(f"{OUT_BASE}.png", dpi=300)
    written.append(f"{OUT_BASE}.png")
print("wrote " + ", ".join(written))
