"""Shared design tokens and vector-output settings for the figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- light-surface tokens -----------------------------------------------------
SURFACE  = "#fcfcfb"
INK      = "#0b0b0b"
INK_2    = "#52514e"
MUTED    = "#898781"
GRID     = "#e1e0d9"
BASELINE = "#c3c2b7"

# sequential blue ramp, steps 250-700: step 250 is the lightest that still
# clears 2:1 on the light surface, so the palest mark still reads as a mark
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQ_BLUE = LinearSegmentedColormap.from_list("seq_blue", RAMP)

# diverging pair for polarity (gain vs loss); saturated at both ends so a
# near-zero value stays visible instead of fading into the surface
POS = "#2a78d6"
NEG = "#e34948"


def apply_rcparams():
    """Live text in SVG, embedded TrueType in PDF, no mathtext anywhere.

    svg.fonttype="none" keeps labels editable in Illustrator/Inkscape; mathtext
    would not survive it, so every figure spells out "log2" rather than using
    a math expression. Liberation Sans is metric-compatible with Arial.
    """
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
    })


def tidy(ax, keep_left=False):
    """Recessive chrome: hairline x-grid, one baseline, muted ticks."""
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    if keep_left:
        ax.spines["left"].set_color(BASELINE)
        ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_color(BASELINE)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=3, color=BASELINE)
