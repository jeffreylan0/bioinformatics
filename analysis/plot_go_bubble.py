#!/usr/bin/env python3
"""
Production figure: GO over-representation in the 1-month Ddx41 cKO retina.

Biological-process and cellular-component terms share one plot per direction;
the ontology is carried on the y-axis as a bracketed group label.

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
BAND_GAP = 0.9               # blank rows between the two ontology bands

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
NS_ORDER = [ns for ns in ("BP", "CC") if ns in set(d.ontology)]   # top to bottom
NS_FULL = {"BP": "Biological process", "CC": "Cellular component"}

# Rows are laid out once for every panel so the ontology bands line up across
# directions: the fullest panel sets each band's height.
band_n = {ns: max(len(d[(d.direction == x) & (d.ontology == ns)]) for x in DIRS)
          for ns in NS_ORDER}
n_rows = sum(band_n.values())

# y ranges per band, counting down from the top
band_top = {}
cursor = n_rows + BAND_GAP * (len(NS_ORDER) - 1)
for ns in NS_ORDER:
    band_top[ns] = cursor - 1          # y of the band's first (highest) row
    cursor -= band_n[ns] + BAND_GAP
Y_MAX = n_rows + BAND_GAP * (len(NS_ORDER) - 1) - 1

norm = Normalize(vmin=d.neglog10FDR.min(), vmax=d.neglog10FDR.max())
xmax = d.fold_enrichment.max() * 1.18

ROW_IN = 0.52                      # inches per term row (two-line labels)
HEADER_IN, FOOTER_IN = 1.70, 1.00
FIG_H = ROW_IN * (n_rows + BAND_GAP) + HEADER_IN + FOOTER_IN
FIG_W = 14.5
LEFT, RIGHT, WSPACE = 0.215, 0.90, 0.95

fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(1, 2, left=LEFT, right=RIGHT,
                      top=1 - HEADER_IN / FIG_H, bottom=FOOTER_IN / FIG_H,
                      wspace=WSPACE)

axes, label_artists = [], []
for c, direction in enumerate(DIRS):
    ax = fig.add_subplot(gs[0, c])
    axes.append(ax)
    trans = ax.get_yaxis_transform()
    panel_labels = []

    for ns in NS_ORDER:
        sub = (d[(d.direction == direction) & (d.ontology == ns)]
               .sort_values("fold_enrichment", ascending=False))
        y = band_top[ns] - np.arange(len(sub))     # highest enrichment on top

        ax.scatter(sub.fold_enrichment, y,
                   s=sub.n_DE_in_term * AREA_PER_GENE,
                   c=sub.neglog10FDR, cmap=CMAP, norm=norm,
                   linewidths=0.9, edgecolors=SURFACE, zorder=3)

        # term name with its GO accession set underneath, both right-aligned
        for yi, term, go_id in zip(y, sub.term, sub.GO_ID):
            panel_labels.append(
                ax.annotate(term, xy=(-0.018, yi), xycoords=trans,
                            xytext=(0, 1.5), textcoords="offset points",
                            ha="right", va="bottom", fontsize=9, color=INK))
            ax.annotate(go_id, xy=(-0.018, yi), xycoords=trans,
                        xytext=(0, -2.5), textcoords="offset points",
                        ha="right", va="top", fontsize=7.2, color=MUTED)

    label_artists.append(panel_labels)

    ax.set_yticks([])
    ax.set_ylim(-0.65, Y_MAX + 0.65)
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
    ax.set_title(direction, fontsize=11.5, color=INK, fontweight="bold", pad=10)
    ax.set_xlabel("Fold enrichment", fontsize=9.5, color=INK_2, labelpad=7)

# --- ontology brackets on the y-axis -----------------------------------------
# Placed just clear of the widest term label, which is measured from the drawn
# figure rather than guessed, so the bracket never lands on top of a name.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for ax, labels in zip(axes, label_artists):
    inv = ax.transAxes.inverted()
    left_edge = min(inv.transform((t.get_window_extent(renderer).x0, 0))[0]
                    for t in labels)
    rule_x = left_edge - 0.035
    trans = ax.get_yaxis_transform()
    for ns in NS_ORDER:
        hi, lo = band_top[ns], band_top[ns] - band_n[ns] + 1
        ax.plot([rule_x, rule_x], [lo - 0.42, hi + 0.42], transform=trans,
                color=BASELINE, lw=1.2, solid_capstyle="butt",
                clip_on=False, zorder=2)
        ax.annotate(NS_FULL[ns], xy=(rule_x - 0.018, 0.5 * (hi + lo)),
                    xycoords=trans, ha="center", va="center", rotation=90,
                    fontsize=10, color=INK, fontweight="bold")

# --- colour bar ---------------------------------------------------------------
# Drawn by hand rather than with fig.colorbar: matplotlib's colorbar is a quad
# mesh that is rasterised by default, and un-rasterising it leaves hairline
# seams between the quads in PDF/SVG. Overlapping rectangles cannot seam, and
# keep the figure fully vector.
cax = fig.add_axes([0.925, 0.34, 0.0125, 0.30])
edges = np.linspace(norm.vmin, norm.vmax, 257)
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
             x=0.052, y=1 - 0.40 / FIG_H, ha="left", fontsize=15, color=INK,
             fontweight="bold")
fig.text(0.052, 1 - 0.70 / FIG_H,
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
