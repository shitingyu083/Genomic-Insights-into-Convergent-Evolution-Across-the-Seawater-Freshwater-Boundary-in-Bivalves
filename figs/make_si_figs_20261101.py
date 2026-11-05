# R1 revision, option A (2026-11-01): supplementary figures with the R1 checks (main Figs 2-5 keep the R0 layout).
#   Fig. S2  copy numbers of the 20 gene families of Fig. 3a by salinity class, with the main-model PGLS
#            (genefam20/genefam20_counts.tsv; genefam20/out/genefam20_pgls_results.tsv, level L2, model "1 main").
#   Fig. S3  assembly-span contrast: species values by salinity class and the PGLS estimate as a function of
#            Pagel's lambda (pgls_span/span_lambda_curve.tsv, span_reference_tests.tsv; revised Table S1).
#   Fig. S4  amino-acid composition checks: euhaline-only, BUSCO >= 85% and one-to-one orthologue estimates,
#            GC4 vs proteome GC index and GC4-adjusted estimates (iqtree_rerun/out/rerun_fixed/aa_*.tsv).
# All numbers are read from the analysis outputs (PGLS fits on the revised Table S1 values, 2026-11-01).
# Output: SUBMIT_R1_revised_20261024/Figures_R1/supplementary/FigS2-S4_R1.{png,pdf} + source data.
import os, textwrap
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

R = os.environ.get("IZ_ROOT", ".")
W = R + "/_R1_work"
AN = W + "/analysis"
PG = AN + "/iqtree_rerun/out/rerun_fixed"
OUT = R + "/SUBMIT_R1_revised_20261024/Figures_R1/supplementary"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.labelsize": 12, "axes.titlesize": 12,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
                     "axes.edgecolor": "#222222", "xtick.color": "#222222", "ytick.color": "#222222",
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold", "axes.unicode_minus": True})
DK, MD, LT = "#222222", "#888888", "#CCCCCC"
RED, BLUE, LBLUE, TEAL, PURPLE = "#ca2a28", "#2371a9", "#8cc2ea", "#00a79d", "#8a92c2"
PBLUE, PPINK = "#c8e2eb", "#f6dcd9"
HAB = {"Freshwater": RED, "Brackish": LBLUE, "Euhaline": BLUE}
TICK = {"Freshwater": RED, "Brackish": "#4a8fc4", "Euhaline": BLUE}
GRP = {"FRESHWATER": "Freshwater", "BRACKISH": "Brackish", "MARINE": "Euhaline"}
MINUS = "\u2212"


def tag(ax, s, x=-0.14):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def habitat_strip(ax, d, i, col, rng, s=34):
    q1, med, q3 = np.percentile(d, [25, 50, 75])
    ax.add_patch(plt.Rectangle((i - 0.3, q1), 0.6, q3 - q1, facecolor=col, alpha=0.18, edgecolor=col, lw=0.8, zorder=1))
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(d)), d, s=s, facecolor=col, edgecolor="white", lw=0.6, zorder=3, alpha=0.95)
    ax.hlines(med, i - 0.3, i + 0.3, color=col, lw=2.2, zorder=4)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

def _sig2(x):
    """x to 2 significant figures, keeping trailing zeros (0.0010, 0.010, 0.25)."""
    x = float(f"{x:.2g}")
    k = int(np.floor(np.log10(abs(x))))
    return f"{x:.{max(0, 1 - k)}f}"

def fnum(x):
    """P or q value for figure text: 2 significant figures, a.b x 10^-k below 0.001."""
    if x < 1e-3:
        m, e = f"{x:.1e}".split("e")
        return f"{m} × 10$^{{-{int(e[1:])}}}$"
    return _sig2(x)

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------------------------------------------------------------- Fig. S2
FAMS = ["ATG8", "BECN1", "AQP", "ABC", "P450", "CTS", "VATP", "CA", "VWA", "CTLD",
        "SLC", "NHE", "HSP70", "NLR", "CBM14", "BIRC", "FREP", "HSP90", "GST", "PDI"]
LABEL = ["ATG8/LC3/GABARAP", "Beclin 1", "Aquaporin/AQP", "ABC transporter", "Cytochrome P450",
         "Cysteine cathepsin", "V-type H$^{+}$-ATPase", "Carbonic anhydrase", "VWA/vWFA-domain protein",
         "CTLD-containing protein", "SLC transporter", "Na$^{+}$/H$^{+}$ exchanger", "HSP70/HSPA",
         "NLR/NACHT-domain protein", "Peritrophin-A/chitin-binding protein", "IAP/BIRC", "FREP/FReD protein",
         "HSP90/HSPC", "GST", "PDI/P4HB"]
cn = pd.read_csv(AN + "/genefam20/genefam20_counts.tsv", sep="\t")
use = cn[(cn.salinity != "OUTGROUP") & (cn.pfam_route_complete == 1)]
assert set(use.salinity) <= set(GRP), set(use.salinity)
assert not set(use.species) & {"Sco", "Mvi"}
gp = pd.read_csv(AN + "/genefam20/out/genefam20_pgls_results.tsv", sep="\t")
gp = gp[(gp.level == "L2") & gp.model.str.startswith("1 main")].set_index("family").loc[FAMS]
assert (gp.n == 40).all() and (gp.n_FW == 8).all(), gp[["n", "n_FW"]]
G3 = ("Freshwater", "Brackish", "Euhaline")
fig, axs = plt.subplots(4, 5, figsize=(17.5, 15))
rng = np.random.default_rng(7)
src = []
for k, (f, lab) in enumerate(zip(FAMS, LABEL)):
    ax = axs.flat[k]
    d = use[(use.family == f) & use.L2.notna()]
    assert len(d) == int(gp.loc[f, "n"]), (f, len(d))
    cnt = d.salinity.map(GRP).value_counts().to_dict()
    assert cnt == {"Freshwater": 8, "Brackish": 10, "Euhaline": 22}, (f, cnt)
    for i, g in enumerate(G3):
        habitat_strip(ax, d.loc[d.salinity.map(GRP) == g, "L2"].to_numpy(float), i, HAB[g], rng, s=30)
    r = gp.loc[f]
    sig = r.q_BH < 0.05
    if sig:
        ax.set_facecolor("#fbeeec")
    ax.set_title("\n".join(textwrap.wrap(lab, 26)), fontsize=12, fontweight="bold", color=RED if sig else DK, pad=22)
    ax.text(0.5, 1.02, f"ratio {r.ratio_FW:.2f} ({r.ratio_lo:.2f}–{r.ratio_hi:.2f}), q = {fnum(r.q_BH)}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color=RED if sig else DK)
    ax.set_xlim(-0.6, 2.6)
    ax.set_xticks(range(3))
    if k >= 15:
        ax.set_xticklabels([f"{g}\n(n = {cnt[g]})" for g in G3])
        for t, g in zip(ax.get_xticklabels(), G3):
            t.set_color(TICK[g])
    else:
        ax.set_xticklabels([])
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    if k % 5 == 0:
        ax.set_ylabel("Copy number (genes)")
    for _, row in d.iterrows():
        src.append({"family": f, "family_label": lab.replace("$^{+}$", "+"), "species": row.species,
                    "salinity_class": GRP[row.salinity], "copy_number_L2": int(row.L2)})
hd = [Line2D([], [], ls="", marker="o", ms=8, mfc=HAB[g], mec="white", label=f"{g} (n = {n})")
      for g, n in zip(G3, (8, 10, 22))]
hd.append(Patch(facecolor="#fbeeec", edgecolor=LT, label="main-model q < 0.05"))
fig.legend(handles=hd, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0), fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.965))
save(fig, "FigS2_R1")
pd.DataFrame(src).to_csv(OUT + "/FigS2_source_data.tsv", sep="\t", index=False)
gp.reset_index()[["family", "level", "model", "n", "n_FW", "lambda", "beta_FW", "se", "P", "ratio_FW", "ratio_lo",
                  "ratio_hi", "q_BH"]].to_csv(OUT + "/FigS2_pgls_main_model.tsv", sep="\t", index=False)
print("Fig S2: main-model q < 0.05 in", [f for f in FAMS if gp.loc[f, "q_BH"] < 0.05])

# ---------------------------------------------------------------- Fig. S3
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
s1 = pd.read_csv(R + "/SUBMIT_R1_revised_20261024/Supplementary_Table_S1_retained_bivalve_genomes_R1.tsv", sep="\t")
s1 = s1[s1["Habitat group in analyses"] != "Outgroup (not in habitat contrasts)"].copy()
s1["cls"] = s1["Salinity class (sensitivity analyses)"]
assert s1.cls.value_counts().to_dict() == {"Euhaline": 23, "Brackish": 11, "Freshwater": 8}
cv = pd.read_csv(AN + "/pgls_span/span_lambda_curve.tsv", sep="\t").sort_values("lambda")
rt = pd.read_csv(AN + "/pgls_span/span_reference_tests.tsv", sep="\t").set_index("test")
pm = pd.read_csv(PG + "/pgls_span_results.tsv", sep="\t").iloc[0]
assert pm.model.startswith("main") and pm.n == 42 and pm.n_FW == 8
mk = cv[cv.mark.notna()].set_index("mark")
ml = mk.loc["ML estimate (main model)"]
assert abs(ml["lambda"] - pm["lambda"]) < 1e-3 and abs(ml.P - pm.P) < 1e-4 and abs(ml.ratio - pm.ratio_FW) < 1e-3
assert (rt.n == 42).all() and (rt.n_FW == 8).all() and np.allclose(cv.ratio, 10 ** cv.beta)
inside = cv[cv.logLik_ML >= cv.logLik_ML.max() - stats.chi2.ppf(0.95, 1) / 2]
lam_lo, lam_hi = inside["lambda"].min(), inside["lambda"].max()
print(f"Fig S3: 95% likelihood interval for lambda {lam_lo:.2f}-{lam_hi:.2f}; "
      f"P inside {inside.P.min():.4f}-{inside.P.max():.3f}; ratio inside {inside.ratio.min():.2f}-{inside.ratio.max():.2f}")
tcrit = stats.t.ppf(0.975, int(pm.n) - 2)
cv["ratio_lo"], cv["ratio_hi"] = 10 ** (cv.beta - tcrit * cv.se), 10 ** (cv.beta + tcrit * cv.se)
assert abs(cv.loc[cv.mark == "ML estimate (main model)", "ratio_lo"].iloc[0] - pm.ratio_lo) < 1e-3

fig = plt.figure(figsize=(16.5, 5.4))
gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.25, 1.2], wspace=0.34)
ax = fig.add_subplot(gs[0])
rng = np.random.default_rng(11)
for i, g in enumerate(G3):
    habitat_strip(ax, s1.loc[s1.cls == g, "Assembly span (Mb)"].to_numpy(float), i, HAB[g], rng)
ax.set_yscale("log")
v_all = s1["Assembly span (Mb)"].astype(float)
ax.yaxis.set_major_locator(FixedLocator([t for t in (200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000)
                                         if v_all.min() / 1.15 <= t <= v_all.max() * 1.15]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_xlim(-0.6, 2.6)
ax.set_xticks(range(3), [f"{g}\n(n = {(s1.cls == g).sum()})" for g in G3])
for t, g in zip(ax.get_xticklabels(), G3):
    t.set_color(TICK[g])
ax.set_ylabel("Assembly span (Mb)")
tag(ax, "a", x=-0.26)

MARK = [("no phylogenetic correction", "s", MD, "λ = 0"),
        ("original PGLS run (0.442)", "D", PURPLE, "λ = 0.442 (original PGLS run)"),
        ("ML estimate (main model)", "o", DK, f"λ = {ml['lambda']:.2f} (ML estimate, main model)"),
        ("Brownian motion", "^", BLUE, "λ = 1 (Brownian motion)")]
for j, key in ((1, "P"), (2, "ratio")):
    ax = fig.add_subplot(gs[j])
    ax.axvspan(lam_lo, lam_hi, color=PBLUE, alpha=0.6, lw=0, zorder=0)
    if key == "ratio":
        ax.fill_between(cv["lambda"], cv.ratio_lo, cv.ratio_hi, color=MD, alpha=0.25, lw=0, zorder=1)
        ax.axhline(1, color=MD, ls="--", lw=1)
    else:
        ax.axhline(0.05, color=RED, ls="--", lw=1, zorder=1)
        ax.text(0.012, 0.044, "P = 0.05", color=RED, fontsize=10, va="top")
    ax.plot(cv["lambda"], cv[key], color=DK, lw=1.6, zorder=2)
    for mkey, m, c, lab in MARK:
        r = cv[cv.mark == mkey].iloc[0]
        if key == "ratio":
            ax.errorbar(r["lambda"], r.ratio, yerr=[[r.ratio - r.ratio_lo], [r.ratio_hi - r.ratio]], fmt=m, ms=8,
                        mfc=c, mec="white", mew=0.8, ecolor=c, elinewidth=1.2, capsize=3, zorder=4)
        else:
            ax.scatter(r["lambda"], r.P, marker=m, s=75, facecolor=c, edgecolor="white", lw=0.8, zorder=4,
                       label=f"{lab}: P = {fnum(r.P)}")
    ax.set_xlim(-0.03, 1.03)
    ax.set_xlabel("Pagel's λ (fixed in the PGLS fit)")
    if key == "P":
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 3)
        ax.set_ylabel("P (freshwater effect)")
        hb = ax.get_legend_handles_labels()[0] + [Patch(facecolor=PBLUE, alpha=0.6, label=(
            f"95% likelihood interval for λ ({lam_lo:.2f}–{lam_hi:.2f})"))]
        ax.legend(handles=hb, frameon=False, loc="upper left", fontsize=10, handletextpad=0.4, borderaxespad=0.2)
        ax.text(0.97, 0.04, "Without phylogenetic correction\n(8 freshwater vs 34 comparators):\n"
                f"Wilcoxon rank-sum P = {fnum(rt.loc['Wilcoxon rank-sum (span)', 'P'])}\n"
                f"Welch t (log$_{{10}}$ span) P = {fnum(rt.loc['Welch t (log10 span)', 'P'])}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=10, linespacing=1.3)
        tag(ax, "b", x=-0.16)
    else:
        ax.set_ylabel("Freshwater/comparator span ratio\n(PGLS estimate, 95% CI)")
        tag(ax, "c", x=-0.16)
save(fig, "FigS3_R1")
s1.rename(columns={"cls": "Salinity class"})[["ID", "Species", "Habitat group in analyses", "Salinity class",
                                              "Assembly span (Mb)"]].to_csv(OUT + "/FigS3a_source_data.tsv",
                                                                            sep="\t", index=False)
cv.assign(in_95pct_likelihood_interval=cv.logLik_ML >= cv.logLik_ML.max() - stats.chi2.ppf(0.95, 1) / 2).to_csv(
    OUT + "/FigS3bc_source_data.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- Fig. S4
GROUPS = [("Small/special", "GAP", "#FBE4D5"), ("Hydrophobic", "VLIMFW", "#E6EEF9"),
          ("Polar", "STCNQY", "#ECF6E1"), ("Charged", "DEKRH", "#EEE7F5")]  # residue order and groups of Fig. 5
NAME = {"G": "Gly", "A": "Ala", "P": "Pro", "V": "Val", "L": "Leu", "I": "Ile", "M": "Met", "F": "Phe", "W": "Trp",
        "S": "Ser", "T": "Thr", "C": "Cys", "N": "Asn", "Q": "Gln", "Y": "Tyr", "D": "Asp", "E": "Glu", "K": "Lys",
        "R": "Arg", "H": "His"}
AA = "".join(g[1] for g in GROUPS)
ap = pd.read_csv(PG + "/aa_pgls_results.tsv", sep="\t")
SETS = [("ALL_ANY", "1 main", "Main model (n = 42)", "o", RED, 42, 8),
        ("ALL_ANY", "2 comparator", "Euhaline comparators only (n = 31)", "s", BLUE, 31, 8),
        ("ALL_ANY", "3 BUSCO", "BUSCO C ≥ 85% (n = 36)", "^", TEAL, 36, 6),
        ("SC8_OGC", "1 main", "One-to-one orthologues (n = 42)", "D", PURPLE, 42, 8)]
fig = plt.figure(figsize=(16, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[2.3, 1, 1], wspace=0.3)
ax = fig.add_subplot(gs[0])
x = np.arange(20)
srcA = []
for j, (ds, mod, lab, m, c, n, nfw) in enumerate(SETS):
    d = ap[(ap.dataset == ds) & ap.model.str.startswith(mod)].set_index("aa").loc[list(AA)]
    assert (d.n == n).all() and (d.n_FW == nfw).all(), (ds, mod)
    if j == 0:
        assert d.q_BH.min() >= 0.097 - 1e-3, d.q_BH.min()
    ax.errorbar(x + (j - 1.5) * 0.19, d.beta_FW, yerr=[d.beta_FW - d.ci_lo, d.ci_hi - d.beta_FW], fmt=m, ms=5.5,
                mfc=c, mec="white", mew=0.5, ecolor=c, elinewidth=0.9, capsize=0, label=lab, zorder=3)
    srcA.append(d.reset_index()[["dataset", "model", "aa", "n", "n_FW", "lambda", "beta_FW", "se", "ci_lo", "ci_hi",
                                 "P", "q_BH"]])
srcA = pd.concat(srcA)
for g, s, col in GROUPS:
    i0, i1 = AA.index(s[0]), AA.index(s[-1])
    ax.axvspan(i0 - 0.5, i1 + 0.5, color=col, lw=0, zorder=0)
    ax.text((i0 + i1) / 2, 1.01, g, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=10)
ax.axhline(0, color=MD, ls="--", lw=1, zorder=1)
ax.set_xlim(-0.5, 19.5)
ax.set_xticks(x, [NAME[a] for a in AA])
for t, a in zip(ax.get_xticklabels(), AA):
    if a in "GAPQ":
        t.set_fontweight("bold")
span_y = srcA.ci_hi.max() - srcA.ci_lo.min()
ax.set_ylim(srcA.ci_lo.min() - 0.08 * span_y, srcA.ci_hi.max() + 0.08 * span_y)
ax.set_ylabel("β$_{FW}$ (CLR units, 95% CI)")
ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.07), fontsize=10)
tag(ax, "a", x=-0.09)

sp = pd.read_csv(W + "/server_data/r1rev/aa_gc/aa_gc_species.tsv", sep="\t")
sp = sp[sp.matches_proteome_count == 1]
cl = pd.read_csv(AN + "/aa/aa_species_clr.tsv", sep="\t")
cl = cl[cl.dataset == "ALL_ANY"][["species", "salinity3", "gcp"]]
d = sp.merge(cl, on="species")
assert len(d) == 10 and (d.salinity3 == "FRESHWATER").sum() == 5 and set(d.salinity3) <= set(GRP)
rho = pd.read_csv(PG + "/aa_gc_proxy_cor.tsv", sep="\t")
rho = rho[(rho.set == "S10") & (rho.x == "GC4") & (rho.y == "gcp")].iloc[0]
assert abs(stats.spearmanr(d.GC4, d.gcp).correlation - rho.rho) < 1e-3
ax = fig.add_subplot(gs[1])
fit = np.polyfit(d.GC4, d.gcp, 1)
xx = np.linspace(d.GC4.min(), d.GC4.max(), 50)
ax.plot(xx, np.polyval(fit, xx), color=MD, lw=1, zorder=1)  # visual guide only; the statistic is Spearman's rho
for g in G3:
    e = d[d.salinity3.map(GRP) == g]
    if len(e):
        ax.scatter(e.GC4, e.gcp, s=62, facecolor=HAB[g], edgecolor="white", lw=0.8, label=f"{g} (n = {len(e)})", zorder=3)
        for _, r in e.iterrows():
            off = {"Lfo": (6, -4), "Sbr": (-4, 7)}.get(r.species, (4, 3))  # Lfo and Sbr lie close together
            ax.annotate(r.species, (r.GC4, r.gcp), xytext=off, textcoords="offset points", fontsize=10, color=TICK[g],
                        ha="right" if off[0] < 0 else "left", va="center" if r.species == "Lfo" else "baseline")
ax.set(xlabel="GC4", ylabel="Proteome GC index ln[GARP/FYMINK]")
ax.set_xlim(d.GC4.min() - 0.025, d.GC4.max() + 0.03)
ax.text(0.97, 0.05, f"Spearman ρ = {rho.rho:.3f}\nP = {fnum(rho.P)}, n = {int(rho.n)}", transform=ax.transAxes,
        ha="right", fontsize=10)
ax.legend(frameon=False, loc="upper left")
tag(ax, "b", x=-0.22)

gt = pd.read_csv(PG + "/aa_gc_tests.tsv", sep="\t")
gt = gt[gt.set == "S10"]
ax = fig.add_subplot(gs[2])
srcC = []
for j, (mod, lab, c) in enumerate((("FW", "Without GC4", RED), ("FW+GC4", "With GC4 as covariate", BLUE))):
    e = gt[gt.model == mod].set_index("response").loc[["clr_G", "clr_A", "clr_P", "clr_Q"]]
    assert (e.n == 10).all() and (e.n_FW == 5).all(), mod
    ax.errorbar(np.arange(4) + (j - 0.5) * 0.25, e.beta_FW, yerr=[e.beta_FW - e.ci_lo, e.ci_hi - e.beta_FW], fmt="o",
                ms=7, mfc=c, mec="white", mew=0.7, ecolor=c, elinewidth=1.4, capsize=3, label=lab, zorder=3)
    srcC.append(e.reset_index())
srcC = pd.concat(srcC)
ax.axhline(0, color=MD, ls="--", lw=1)
ax.set_xticks(range(4), ["Gly", "Ala", "Pro", "Gln"])
ax.set(ylabel="β$_{FW}$ (CLR units, 95% CI)", xlabel="10 species (5 freshwater)")
ax.set_ylim(min(-0.105, srcC.ci_lo.min() - 0.01), max(0.13, srcC.ci_hi.max() + 0.045))
ax.legend(frameon=False, loc="upper right")
tag(ax, "c", x=-0.22)
save(fig, "FigS4_R1")
srcA.to_csv(OUT + "/FigS4a_source_data.tsv", sep="\t", index=False)
d[["species", "salinity3", "GC4", "gcp"]].to_csv(OUT + "/FigS4b_source_data.tsv", sep="\t", index=False)
srcC.to_csv(OUT + "/FigS4c_source_data.tsv", sep="\t", index=False)
print("done", sorted(os.listdir(OUT)))
