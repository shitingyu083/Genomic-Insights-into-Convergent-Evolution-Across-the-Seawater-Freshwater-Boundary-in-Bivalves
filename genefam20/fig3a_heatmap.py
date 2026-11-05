# Option A (2026-10-31): redraw Fig. 3a in the R0 design from the real 20-family copy numbers.
# Layout, row order, row labels, module bar and colours are copied from the R0 panel as drawn
# (Main_Figures_1-6.pdf page 3, sampled at 4x); values = row z-scores of ln(L2 gene count + 1)
# across the 40 counted species; stars = main-model PGLS q (BH over the 20 families, L2).
# Output: SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig3/Fig3a_heatmap_R1.{png,pdf} + source data.
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

R = os.environ.get("IZ_ROOT", ".")
W = R + "/_R1_work/analysis/genefam20"
OUT = R + "/SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig3"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "font.weight": "bold", "pdf.fonttype": 42})

FAMS = ["ATG8", "BECN1", "AQP", "ABC", "P450", "CTS", "VATP", "CA", "VWA", "CTLD",
        "SLC", "NHE", "HSP70", "NLR", "CBM14", "BIRC", "FREP", "HSP90", "GST", "PDI"]
LABEL = ["ATG8/LC3/GABARAP", "Beclin 1", "Aquaporin/AQP", "ABC transporter", "Cytochrome P450",
         "Cysteine cathepsin", "V-type H$^{+}$-ATPase", "Carbonic anhydrase", "VWA/vWFA-domain protein",
         "CTLD-containing protein", "SLC transporter", "Na$^{+}$/H$^{+}$ exchanger", "HSP70/HSPA",
         "NLR/NACHT-domain protein", "Peritrophin-A/chitin-binding protein", "IAP/BIRC", "FREP/FReD protein",
         "HSP90/HSPC", "GST", "PDI/P4HB"]
MODULES = [("Autophagy and lysosomal turnover", "#00A087"), ("Transmembrane transport and detoxification", "#F39B7F"),
           ("Osmoregulation and ion homeostasis", "#8491B4"), ("Acid-base balance and biomineralization", "#91D1C2"),
           ("Extracellular matrix and adhesion", "#DCAFF1"), ("Innate immunity and host defense", "#B09C85"),
           ("Protein homeostasis and ER proteostasis", "#7E6148"), ("Apoptosis and cellular stress response", "#3C5488"),
           ("Oxidative stress defense", "#EFC000")]
# module of each row as drawn in R0 (index into MODULES)
ROWMOD = [0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 2, 2, 6, 5, 4, 7, 5, 6, 8, 6]
FW_COLS = ["Pst", "Cfl", "Dpo", "Dro", "Lfo", "Sol", "Awo", "Ude"]
SW_COLS = ["Sbr", "Mye", "Cgi", "Sco", "Tna", "Mco", "Csi", "Bpl", "Mve", "Mme", "Sgl", "Cho", "Rph", "Aja", "Apl",
           "Ama", "Air", "Apu", "Cfa", "Cbi", "Can", "Cni", "Cvi", "Mar", "Mvi", "Mca", "Mga", "Oed", "Pmax", "Pfu",
           "Pim", "Pmag", "Sso", "Yba"]
HAB_FW, HAB_SW, NA_COL = "#B5222F", "#286EAE", "#CCCCCC"
# colour ramp anchors sampled from the R0 colour bar (value: colour)
RAMP = [(-1.8, "#2C71AF"), (-1.35, "#5B96C5"), (-0.9, "#A2D1E0"), (-0.45, "#E0EDF1"), (0, "#FCFCFC"),
        (0.45, "#EFD2C1"), (0.9, "#E6A391"), (1.35, "#DC6952"), (1.8, "#C63B36")]
cmap = LinearSegmentedColormap.from_list("r0", [((v + 1.8) / 3.6, c) for v, c in RAMP])

# ---- data
cn = pd.read_csv(W + "/genefam20_counts.tsv", sep="\t")
cn = cn[cn.family.isin(FAMS)]
cols = FW_COLS + SW_COLS
assert set(cols) == set(cn[cn.salinity != "OUTGROUP"].species), "column set differs from counts table"
assert set(cn[cn.salinity == "FRESHWATER"].species) == set(FW_COLS)
use = cn[(cn.salinity != "OUTGROUP") & (cn.pfam_route_complete == 1)]
mat = use.pivot(index="family", columns="species", values="L2").loc[FAMS]
assert mat.shape[1] == 40 and not mat.isna().any().any(), mat.shape
lx = np.log(mat + 1)
z = lx.sub(lx.mean(axis=1), axis=0).div(lx.std(axis=1, ddof=1), axis=0)
zc = z.clip(-1.8, 1.8).reindex(columns=cols)          # Sco, Mvi -> NaN (not counted)
notc = [c for c in cols if c not in mat.columns]
assert notc == ["Sco", "Mvi"], notc

res = pd.read_csv(W + "/out/genefam20_pgls_results.tsv", sep="\t")
res = res[(res.level == "L2") & res.model.str.startswith("1 main")].set_index("family").loc[FAMS]
assert (res.n == 40).all() and (res.n_FW == 8).all()
star = lambda q: "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""
stars = [star(q) for q in res.q_BH]

src = (use[use.family.isin(FAMS)][["species", "salinity", "family", "L2"]].copy())
src["ln_L2_plus_1"] = np.log(src.L2 + 1)
src["row_z"] = [z.loc[f, s] for f, s in zip(src.family, src.species)]
src["row_z_plotted"] = src.row_z.clip(-1.8, 1.8)
src = src.merge(res[["ratio_FW", "ratio_lo", "ratio_hi", "P", "q_BH"]].reset_index(), on="family")
src["star"] = src.q_BH.map(star)
src.sort_values(["family", "species"]).to_csv(OUT + "/Fig3a_source_data.tsv", sep="\t", index=False)

# ---- layout in inches (one axes spanning the figure, data units = inches)
cw, ch, gap = 0.155, 0.19, 0.08
x_mod = 0.55
x0 = x_mod + cw + 0.10                                   # first heatmap column
xs = [x0 + i * cw for i in range(8)] + [x0 + 8 * cw + gap + i * cw for i in range(34)]
x_end = xs[-1] + cw
y_top = 4.55                                             # top of first heatmap row
y_hab = y_top + 0.10
FW_, H_ = 14.6, 5.25
fig = plt.figure(figsize=(FW_, H_))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, FW_); ax.set_ylim(0, H_); ax.axis("off")

for i, f in enumerate(FAMS):
    yb = y_top - (i + 1) * ch
    ax.add_patch(Rectangle((x_mod, yb), cw, ch, facecolor=MODULES[ROWMOD[i]][1], edgecolor="white", lw=0.4))
    for j, c in enumerate(cols):
        v = zc.loc[f, c]
        fc = NA_COL if np.isnan(v) else cmap((v + 1.8) / 3.6)
        ax.add_patch(Rectangle((xs[j], yb), cw, ch, facecolor=fc, edgecolor="white", lw=0.4))
    ax.text(x_end + 0.10, yb + ch / 2, LABEL[i] + ("  " + stars[i] if stars[i] else ""), va="center", ha="left")
ax.add_patch(Rectangle((xs[0], y_hab), 8 * cw, 0.16, facecolor=HAB_FW, edgecolor="white", lw=0.4))
ax.add_patch(Rectangle((xs[8], y_hab), 34 * cw, 0.16, facecolor=HAB_SW, edgecolor="white", lw=0.4))
ax.text(x_end + 0.10, y_hab + 0.08, "Habitat", va="center", ha="left")
y_bot = y_top - 20 * ch
for j, c in enumerate(cols):
    ax.text(xs[j] + cw * 0.75, y_bot - 0.06, c, rotation=45, ha="right", va="top", rotation_mode="anchor")
ax.text(x_mod + cw * 0.75, y_bot - 0.06, "Module", rotation=45, ha="right", va="top", rotation_mode="anchor")
ax.text(0.12, H_ - 0.12, "a", fontsize=14, fontweight="bold", va="top")

# colour bar
xc, cbw, cb0, cb1 = 10.15, 0.2, y_bot + 0.45, y_top + 0.12
vals = np.linspace(-1.8, 1.8, 181)
hh = (cb1 - cb0) / len(vals)
for k, v in enumerate(vals):
    ax.add_patch(Rectangle((xc, cb0 + k * hh), cbw, hh * 1.02, facecolor=cmap((v + 1.8) / 3.6), lw=0))
for v, lab in ((1.8, "≥ 1.8"), (0.9, "0.9"), (0, "0"), (-0.9, "−0.9"), (-1.8, "≤ −1.8")):
    ax.text(xc + cbw + 0.08, cb0 + (v + 1.8) / 3.6 * (cb1 - cb0), lab, va="center", ha="left")
ax.add_patch(Rectangle((xc, cb0 - 0.42), cbw, 0.19, facecolor=NA_COL, edgecolor="white", lw=0.4))
ax.text(xc + cbw + 0.08, cb0 - 0.325, "Not counted", va="center", ha="left")

# legends
xl, ps = 11.15, 0.19
ax.text(xl, y_top + 0.12, "Habitat", va="center", ha="left")
for k, (lab, col) in enumerate((("Freshwater", HAB_FW), ("Seawater", HAB_SW))):
    yy = y_top - 0.12 - k * (ps + 0.02)
    ax.add_patch(Rectangle((xl, yy - ps / 2), ps, ps, facecolor=col, lw=0))
    ax.text(xl + ps + 0.07, yy, lab, va="center", ha="left")
ym = y_top - 0.95
ax.text(xl, ym, "Module", va="center", ha="left")
for k, (lab, col) in enumerate(MODULES):
    yy = ym - 0.24 - k * (ps + 0.02)
    ax.add_patch(Rectangle((xl, yy - ps / 2), ps, ps, facecolor=col, lw=0))
    ax.text(xl + ps + 0.07, yy, lab, va="center", ha="left")

for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/Fig3a_heatmap_R1.{ext}", dpi=300)
plt.close(fig)
print("stars:", {f: s for f, s in zip(FAMS, stars) if s})
print("q:", {f: float(f"{q:.3g}") for f, q in zip(FAMS, res.q_BH)})
print("FW-block mean z:", {f: round(float(z.loc[f, FW_COLS].mean()), 2) for f in FAMS})
