# Option A (2026-10-31): redraw main Fig. 5 in the R0 design from the real species-level amino-acid composition.
# R0 layout kept (Main_Figures_1-6.pdf page 5): 2 x 2 panels with framed titles, residues grouped and shaded by
# property (Small/special, Hydrophobic, Polar, Charged), split half-violins (seawater left, freshwater right) with
# box plots and points, shared y label, legend on top, q note at the bottom. Residue labels sit below each axis
# (R0 clipped "Trp (W)", "Cys (C)", "Asn (N)" in panel d).
# Data: aa/aa_species_clr.tsv (species-level CLR; ALL_ANY n = 42, ALL_CYT/ALL_EXT/ALL_MEM n = 40).
# q labels: main-model PGLS (iqtree_rerun/out/rerun_fixed/aa_pgls_results.tsv), BH across 20 residues per set.
# Output: SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig5/Fig5_R1.{png,pdf} + source data.
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

R = os.environ.get("IZ_ROOT", ".")
A = R + "/_R1_work/analysis"
OUT = R + "/SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig5"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "font.weight": "bold", "pdf.fonttype": 42,
                     "ps.fonttype": 42, "axes.linewidth": 0.8, "axes.labelweight": "bold"})
RED, BLUE, DK = "#ca2a28", "#2371a9", "#222222"
GROUPS = [("Small/special", "GAP", "#FBE4D5"), ("Hydrophobic", "VLIMFW", "#E6EEF9"),
          ("Polar", "STCNQY", "#ECF6E1"), ("Charged", "DEKRH", "#EEE7F5")]
NAME = {"G": "Gly", "A": "Ala", "P": "Pro", "V": "Val", "L": "Leu", "I": "Ile", "M": "Met", "F": "Phe", "W": "Trp",
        "S": "Ser", "T": "Thr", "C": "Cys", "N": "Asn", "Q": "Gln", "Y": "Tyr", "D": "Asp", "E": "Glu", "K": "Lys",
        "R": "Arg", "H": "His"}
AA = "".join(g[1] for g in GROUPS)
PANELS = [("a", "ALL_ANY", "All protein models"), ("b", "ALL_CYT", "Cytosolic proteins"),
          ("c", "ALL_EXT", "Extracellular proteins"), ("d", "ALL_MEM", "Membrane proteins")]

d = pd.read_csv(A + "/aa/aa_species_clr.tsv", sep="\t")
res = pd.read_csv(A + "/iqtree_rerun/out/rerun_fixed/aa_pgls_results.tsv", sep="\t")
res = res[res.model.str.startswith("1 main")]

FW_IN, FH_IN = 11.69, 8.27
fig = plt.figure(figsize=(FW_IN, FH_IN))
L0, GAP, RM = 0.075, 0.05, 0.012
CW = (1 - L0 - RM - GAP) / 2
# vertical layout in inches: legend 0.38, per row title 0.30 + 0.06 gap + axes 2.42 + label zone 0.86, 0.20 between
# rows, note below
TT, TG, AH, LZ = 0.30, 0.06, 2.42, 0.86
y1 = FH_IN - 0.38 - TT - TG - AH
y2 = y1 - LZ - 0.20 - TT - TG - AH
ROWS = [(y1 / FH_IN, AH / FH_IN), (y2 / FH_IN, AH / FH_IN)]   # (axes bottom, axes height) per row
rng = np.random.default_rng(20261031)
src, qrows, axes = [], [], []

for i, (tag, ds, title) in enumerate(PANELS):
    x = d[d.dataset == ds]
    fw = x[x.salinity3 == "FRESHWATER"]
    sw = x[x.salinity3.isin(["MARINE", "BRACKISH"])]
    assert len(fw) == 8 and len(x) == (42 if ds == "ALL_ANY" else 40), (ds, len(fw), len(x))
    r_ = res[res.dataset == ds].set_index("aa")
    assert set(r_.index) == set(AA) and (r_.n == len(x)).all()
    col, row = i % 2, i // 2
    bot, hgt = ROWS[row]
    left = L0 + col * (CW + GAP)
    ax = fig.add_axes([left, bot, CW, hgt])
    vals = x[["clr_" + a for a in AA]].values
    lo, hi = vals.min(), vals.max()
    ax.set_ylim(lo - 0.12, hi + 0.30)
    ax.set_xlim(-0.6, len(AA) - 0.4)
    # property shading (extends behind the residue and group labels below the axis)
    xs, spans = 0, []
    for gname, gaa, gcol in GROUPS:
        x0, x1 = xs - 0.5 + 0.03, xs + len(gaa) - 0.5 - 0.03
        ax.add_patch(Rectangle((x0, -LZ / AH), x1 - x0, 1 + LZ / AH, transform=ax.get_xaxis_transform(),
                               facecolor=gcol, edgecolor="none", clip_on=False, zorder=0))
        spans.append((gname, x0, x1))
        xs += len(gaa)
    labs = []
    for j, a in enumerate(AA):
        sv, fv = sw["clr_" + a].values, fw["clr_" + a].values
        for v, side, c, mk in ((sv, -1, BLUE, "^"), (fv, 1, RED, "o")):
            g = np.linspace(v.min(), v.max(), 120)
            k = gaussian_kde(v)(g)
            w = 0.40 * k / k.max()
            ax.fill_betweenx(g, j, j + side * w, facecolor=c, alpha=0.18, edgecolor="none", zorder=2)
            ax.plot(j + side * w, g, color=c, lw=0.8, zorder=3)
            ax.plot([j, j], [v.min(), v.max()], color=c, lw=0.6, zorder=3)
            jx = j + side * rng.uniform(0.07, 0.36, len(v))
            ax.scatter(jx, v, marker=mk, s=11, facecolor=c, edgecolor="none", alpha=0.55, zorder=4)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            iqr = q3 - q1
            wl, wh = v[v >= q1 - 1.5 * iqr].min(), v[v <= q3 + 1.5 * iqr].max()
            bx = j + side * 0.10
            ax.plot([bx, bx], [wl, q1], color=DK, lw=0.7, zorder=5)
            ax.plot([bx, bx], [q3, wh], color=DK, lw=0.7, zorder=5)
            ax.add_patch(Rectangle((bx - 0.045, q1), 0.09, max(iqr, 1e-4), facecolor="white", edgecolor=DK,
                                   lw=0.7, zorder=5))
            ax.plot([bx - 0.045, bx + 0.045], [med, med], color=DK, lw=1.1, zorder=6)
        rr = r_.loc[a]
        sig = rr.q_BH < 0.05
        if sig:
            own = max(sv.max(), fv.max())
            nb = [x["clr_" + AA[k]].max() for k in range(max(0, j - 1), min(len(AA), j + 2))]
            top = max(nb) + 0.06
            ax.plot([j - 0.22, j - 0.22, j + 0.22, j + 0.22], [own + 0.03, top, top, own + 0.03], color=DK, lw=0.8,
                    zorder=6)
            ax.text(j, top + 0.015, "q=%.4f" % rr.q_BH, ha="center", va="bottom", fontsize=10, zorder=6)
            qrows.append((tag, ds, a))
        labs.append(ax.text(j, -0.02, "%s (%s)" % (NAME[a], a), transform=ax.get_xaxis_transform(), rotation=90,
                            ha="center", va="top", fontsize=11, color=RED if sig else DK))
        for sp, hab, v in [(s, "Seawater", y) for s, y in zip(sw.species, sv)] + \
                          [(s, "Freshwater", y) for s, y in zip(fw.species, fv)]:
            src.append((tag, ds, title, a, sp, x.set_index("species").loc[sp, "salinity3"].lower(), hab, v))
    axes.append((ax, spans, labs))
    ax.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", labelsize=11, width=0.8, length=3)
    ax.text(0.012, 0.97, tag, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top", ha="left")
    # framed title strip
    t = fig.add_axes([left, bot + hgt + TG / FH_IN, CW, TT / FH_IN])
    t.set_xticks([]); t.set_yticks([])
    for s in t.spines.values():
        s.set_linewidth(0.8)
    t.text(0.5, 0.5, "%s (n = %d)" % (title, len(x)), ha="center", va="center", fontsize=12)

fig.canvas.draw()
rend = fig.canvas.get_renderer()
for ax, spans, labs in axes:
    ymin = min(l.get_window_extent(rend).y0 for l in labs)            # display units
    yl = ax.transAxes.inverted().transform((0, ymin - 4))[1]           # 4 px below the lowest residue label
    for gname, x0, x1 in spans:
        ax.plot([x0 + 0.05, x1 - 0.05], [yl, yl], transform=ax.get_xaxis_transform(), color=DK, lw=0.8,
                clip_on=False)
        ax.text((x0 + x1) / 2, yl - 0.012, gname, transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=11)
fig.text(0.018, 0.47, "Species-level CLR amino-acid frequency", rotation=90, ha="center", va="center", fontsize=12)
H = [Line2D([], [], ls="none", label="Habitat"),
     Line2D([], [], ls="none", marker="^", ms=7, mfc=BLUE, mec=BLUE, alpha=0.8, label="Seawater"),
     Line2D([], [], ls="none", marker="o", ms=7, mfc=RED, mec=RED, alpha=0.8, label="Freshwater"),
     Line2D([], [], ls="none", label="   AA property")] + \
    [Patch(facecolor=c, edgecolor="none", label=g) for g, _, c in sorted(GROUPS)]
lg = fig.legend(handles=H, loc="upper center", bbox_to_anchor=(0.53, 1.0), ncol=len(H), frameon=False,
                fontsize=12, handlelength=1.2, handletextpad=0.4, columnspacing=1.1)
for tx in lg.get_texts():
    tx.set_fontweight("bold")
fig.text(L0, 0.010, "q labels are shown only for PGLS BH-FDR q < 0.05 (main model; 20 residues per protein set) "
                    "on CLR-transformed amino-acid frequencies.", ha="left", va="bottom", fontsize=10)
for ext in ("png", "pdf"):
    fig.savefig(OUT + "/Fig5_R1." + ext, dpi=300)
plt.close(fig)

pd.DataFrame(src, columns=["panel", "dataset", "protein_set", "residue", "species", "salinity_class", "habitat_plotted",
                           "clr"]).to_csv(OUT + "/Fig5_source_data.tsv", sep="\t", index=False)
keep = ["dataset", "aa", "n", "n_FW", "mean_FW", "mean_comp", "lambda", "beta_FW", "se", "ci_lo", "ci_hi", "t", "P",
        "q_BH"]
res[res.dataset.isin([p[1] for p in PANELS])][keep].to_csv(OUT + "/Fig5_pgls_main_model.tsv", sep="\t", index=False)
print("q-labelled:", qrows)
