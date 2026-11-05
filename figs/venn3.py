# Minimal three-set Venn (equal circles) for matplotlib; region counts are drawn as given (not area-proportional).
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

POS = {"A_only": (-0.62, 0.38), "B_only": (0.62, 0.38), "C_only": (0, -0.78), "AB": (0, 0.52), "AC": (-0.38, -0.2),
       "BC": (0.38, -0.2), "ABC": (0, 0.06)}
CEN = [(-0.35, 0.2), (0.35, 0.2), (0, -0.4)]


def venn3(ax, counts, labels, colors, totals=None, r=0.62, fs=10):
    for (x, y), c in zip(CEN, colors):
        ax.add_patch(Circle((x, y), r, facecolor=c, alpha=0.22, edgecolor=c, lw=1.4))
    for k, (x, y) in POS.items():
        ax.text(x, y, f"{counts[k]:,}", ha="center", va="center", fontsize=fs + 1,
                fontweight="bold" if k == "ABC" else "normal", color="#222222")
    lp = [(-0.78, 0.98, "right"), (0.78, 0.98, "left"), (0, -1.16, "center")]
    for i, ((x, y, ha), lab, c) in enumerate(zip(lp, labels, colors)):
        t = lab if totals is None else f"{lab}\n(n = {totals[i]:,})"
        ax.text(x, y, t, ha=ha, va="center", fontsize=fs, color=c)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.25); ax.set_aspect("equal"); ax.axis("off")
