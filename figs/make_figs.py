# R1 revision: redraw Figs 2-5 from the real reanalysis tables (no simulated values).
# Output: SUBMIT_R1_revised_20261024/Figures_R1/Fig{2,3,4,5}.png (300 dpi) + .pdf (vector).
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.environ.get("IZ_ROOT", ".")
W = R + "/_R1_work"
# 2026-10-29: PGLS results refitted on the IQ-TREE branch lengths (OrthoFinder topology)
PG = W + "/analysis/iqtree_rerun/out/rerun_fixed"
A = W + "/reanalysis/A_convergence_v4"
OUT = R + "/SUBMIT_R1_revised_20261024/Figures_R1"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.labelsize": 12, "axes.titlesize": 12,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
                     "axes.edgecolor": "#222222", "xtick.color": "#222222", "ytick.color": "#222222"})
# Palette sampled from the original R0 figures (Main_Figures_1-6.pdf vector fills):
# freshwater red #ca2a28, seawater blue #2371a9, brackish light blue #8cc2ea, salmon #dd6751,
# BUSCO purple #8a92c2 / teal #00a79d, pale panel blue #c8e2eb / pale pink #f6dcd9.
DK, MD, LT = "#222222", "#888888", "#CCCCCC"
RED, BLUE, LBLUE, SALMON, PURPLE, TEAL = "#ca2a28", "#2371a9", "#8cc2ea", "#dd6751", "#8a92c2", "#00a79d"
PBLUE, PPINK = "#c8e2eb", "#f6dcd9"
HAB = {"Freshwater": ("o", RED), "Brackish": ("o", LBLUE), "Euhaline": ("o", BLUE)}


def tag(ax, s):
    ax.text(-0.14, 1.04, s, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


def habitat_strip(ax, d, i, col, rng, s=34):
    # light box (IQR) + jittered points with white edges, as in the original habitat dots
    q1, med, q3 = np.percentile(d, [25, 50, 75])
    ax.add_patch(plt.Rectangle((i - 0.3, q1), 0.6, q3 - q1, facecolor=col, alpha=0.18, edgecolor=col, lw=0.8, zorder=1))
    ax.scatter(i + rng.uniform(-0.17, 0.17, len(d)), d, s=s, facecolor=col, edgecolor="white", lw=0.6, zorder=3, alpha=0.95)
    ax.hlines(med, i - 0.3, i + 0.3, color=col, lw=2.2, zorder=4)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def short(p):
    return p.replace("(sister Lfo)", "").replace("(sister Cfl)", "").replace("(sister Dreissenidae)", "").replace("|", "–")


# ---------- Fig 2 ----------
a3 = pd.read_csv(A + "/A3_pair_summary.tsv", sep="\t")
a6 = pd.read_csv(A + "/A6_castoe_branch_pairs.tsv", sep="\t")
a7 = json.load(open(A + "/A7_castoe_regression.json"))["best_state"]
a5 = pd.read_csv(A + "/A5_unit_recurrent_candidates.tsv", sep="\t")
assert len(a6[a6.cls == "background"]) == 2846 and len(a5) == 127

fig, axs = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw={"width_ratios": [1, 1, 1]})
ax = axs[0]
sty = {"FW": ("o", RED, "Freshwater pairs (n = 6)"), "FW-control": ("s", PURPLE, "Freshwater–sister pairs (n = 9)"),
       "control": ("^", BLUE, "Sister–sister pairs (n = 3)")}
lim = [0, max(a3.obs_C.max(), a3.exp_C_fgene.max()) * 1.12]
ax.fill_between(lim, lim, [lim[1]] * 2, color=PPINK, alpha=0.45, lw=0, zorder=0)
ax.text(lim[1] * 0.04, lim[1] * 0.93, "Observed > expected", fontsize=10, color=RED, va="top")
ax.plot(lim, lim, ls="--", color=MD, lw=1)
for k, (m, c, lab) in sty.items():
    d = a3[a3.kind == k]
    ax.scatter(d.exp_C_fgene, d.obs_C, marker=m, s=58, facecolor=c, edgecolor="white", lw=0.8, label=lab, zorder=3)
    for _, r in d[d.q_C_fgene < 0.10].iterrows():
        ax.annotate(short(r.pair), (r.exp_C_fgene, r.obs_C), xytext=(-6, 0), textcoords="offset points", fontsize=10,
                    ha="right", va="center", color=c)
ax.set(xlim=lim, ylim=lim, xlabel="Expected C (JTT-fgene)", ylabel="Observed C (parallel + convergent)")
ax.legend(frameon=False, loc="lower right")
tag(ax, "a")

ax = axs[1]
bg = a6[a6.cls == "background"]
ax.scatter(bg.D, bg.C, s=7, color=LBLUE, alpha=0.55, edgecolor="none", label="Background pairs (n = 2,846)", zorder=1)
xs = np.linspace(0, bg.D.max(), 200)
b0, b1 = a7["intercept"], a7["slope"]
foc = pd.DataFrame(a7["focal"])
# OLS 95% prediction interval from the background pairs; checked against the A7 pi95_hi of the six focal pairs
from scipy import stats
n = len(bg); xb = bg.D.mean(); sxx = ((bg.D - xb) ** 2).sum()
res = bg.C - (b0 + b1 * bg.D); s = np.sqrt((res ** 2).sum() / (n - 2)); tq = stats.t.ppf(0.975, n - 2)
pi = lambda v: tq * s * np.sqrt(1 + 1 / n + (v - xb) ** 2 / sxx)
hw = np.abs(b0 + b1 * foc.D + pi(foc.D) - foc.pi95_hi).max()
assert hw < 0.05, hw
ax.fill_between(xs, b0 + b1 * xs - pi(xs), b0 + b1 * xs + pi(xs), color=BLUE, alpha=0.15, lw=0, label="95% prediction interval")
ax.plot(xs, b0 + b1 * xs, color=BLUE, lw=1.4, label=f"Fit (slope {b1:.3f}, r² = {a7['r2']:.2f})")
ax.scatter(foc.D, foc.C, s=64, color=RED, edgecolor="white", lw=0.8, zorder=3, label="Freshwater pairs (n = 6)")
ax.set(xlabel="Divergent substitutions (D)", ylabel="Parallel + convergent (C)")
ax.legend(frameon=False, loc="upper left")
tag(ax, "b")

ax = axs[2]
t = a5.groupby(["units", "type"]).size().unstack(fill_value=0)
for c in ("parallel", "convergent", "mixed"):
    if c not in t: t[c] = 0
t = t.loc[t.sum(axis=1).sort_values().index]
lab = [u.replace(";", " + ") if u.count(";") < 2 else u.replace(";", " +" + chr(10), 1).replace(";", " + ") for u in t.index]
left = np.zeros(len(t))
for c, col in (("parallel", BLUE), ("convergent", SALMON), ("mixed", PURPLE)):
    ax.barh(lab, t[c], left=left, color=col, edgecolor="white", lw=0.6, height=0.7, label=c.capitalize())
    left += t[c].values
for i, v in enumerate(left):
    ax.text(v + 0.5, i, str(int(v)), va="center", fontsize=10)
ax.set(xlabel="Candidate sites (n = 127)")
ax.legend(frameon=False, loc="lower right")
tag(ax, "c")
fig.tight_layout()
save(fig, "Fig2")
print("Fig2 PI max deviation from A7:", round(hw, 4))

# ---------- Fig 3 ----------
gc = pd.read_csv(W + "/analysis/genefam/genefam_counts.tsv", sep="\t")
gp = pd.read_csv(PG + "/genefam_pgls_results.tsv", sep="\t")
gp = gp[(gp.level == "L2") & gp.model.str.startswith("1 main")].set_index("family")
gc = gc[(gc.salinity != "OUTGROUP") & (gc.pfam_route_complete == 1)].dropna(subset=["L2"])
grp = {"FRESHWATER": "Freshwater", "BRACKISH": "Brackish", "MARINE": "Euhaline"}
fams = ["ATG8", "BECN1", "AQP", "ABC", "P450"]
fig, axs = plt.subplots(1, 5, figsize=(15, 4))
rng = np.random.default_rng(1)
for ax, f in zip(axs, fams):
    d = gc[gc.family == f]
    assert len(d) == int(gp.loc[f, "n"]), (f, len(d))
    for i, (s, g) in enumerate(grp.items()):
        y = d[d.salinity == s].L2.values
        habitat_strip(ax, y, i, HAB[g][1], rng, s=30)
    ax.set_xlim(-0.6, 2.6)
    ax.set_xticks(range(3), [f"{g}\n(n = {(d.salinity == s).sum()})" for s, g in grp.items()], rotation=45, ha="right")
    for lbl, (s, g) in zip(ax.get_xticklabels(), grp.items()):
        lbl.set_color(HAB[g][1] if g != "Brackish" else "#4a8fc4")
    r = gp.loc[f]
    sig = r.q_BH < 0.05
    ax.set_title(f"{f}", fontsize=12, fontweight="bold", color=RED if sig else DK, pad=18)
    ax.text(0.5, 1.02, f"ratio {r.ratio_FW:.2f} ({r.ratio_lo:.2f}–{r.ratio_hi:.2f}), q = {r.q_BH:.2g}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, color=RED if sig else DK)
    if sig:
        ax.set_facecolor("#fbeeec")
axs[0].set_ylabel("Copy number (genes)")
fig.tight_layout()
save(fig, "Fig3")

# ---------- Fig 4 ----------
s1 = pd.read_csv(R + "/SUBMIT_R1_revised_20261024/Supplementary_Table_S1_retained_bivalve_genomes_R1.tsv", sep="\t")
s1["span"] = pd.to_numeric(s1["Assembly span (Mb)"], errors="coerce")
s1 = s1[s1["Habitat group in analyses"] != "Outgroup (not in habitat contrasts)"].dropna(subset=["span"])
cls = s1["Salinity class (sensitivity analyses)"].str.lower()
s1["g"] = np.where(cls.str.startswith("fresh"), "Freshwater", np.where(cls.str.startswith("brack"), "Brackish", "Euhaline"))
cnt = s1.g.value_counts()
assert (cnt["Freshwater"], cnt["Brackish"], cnt["Euhaline"]) == (8, 11, 23), cnt  # 2026-11-02: Nov-1 Table S1, Teredo navalis in brackish
ps = pd.read_csv(W + "/analysis/pgls_span/pgls_span_results_R1.tsv", sep="\t")  # 2026-11-02: refit on Nov-1 Table S1 (IQ-TREE tree)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.5]})
ax = axs[0]
for i, g in enumerate(["Freshwater", "Brackish", "Euhaline"]):
    y = s1[s1.g == g].span.values
    habitat_strip(ax, y, i, HAB[g][1], rng, s=44)
ax.set_xlim(-0.6, 2.6)
ax.set_yscale("log")
ax.set_xticks(range(3), [f"{g}\n(n = {cnt[g]})" for g in ["Freshwater", "Brackish", "Euhaline"]])
ax.set_ylabel("Assembly span (Mb)")
tag(ax, "a")
ax = axs[1]
y = np.arange(len(ps))[::-1]
for yi, (_, r) in zip(y, ps.iterrows()):
    c = RED if r.model.startswith("main") else BLUE
    if yi % 2 == 0:
        ax.axhspan(yi - 0.5, yi + 0.5, color=PBLUE, alpha=0.35, lw=0, zorder=0)
    ax.errorbar(r.ratio_FW, yi, xerr=[[r.ratio_FW - r.ratio_lo], [r.ratio_hi - r.ratio_FW]], fmt="o", color=c,
                ecolor=c, capsize=3, ms=7 if c == RED else 6, mec="white", mew=0.7, elinewidth=1.4, zorder=3)
ax.set_ylim(-0.5, len(ps) - 0.5)
ax.axvline(1, ls="--", color=MD)
ax.set_yticks(y, [m.replace("main: FW vs brackish+euhaline", "Main model")[:1].upper() + m.replace("main: FW vs brackish+euhaline", "Main model")[1:] for m in ps.model])
ax.set_xticks([0.5, 1, 1.5, 2, 2.5])
xmax = ps.ratio_hi.max() * 1.05
for yi, (_, r) in zip(y, ps.iterrows()):
    ax.text(xmax, yi, f"n = {r.n} ({r.n_FW} FW), P = {r.P:.3f}", va="center", fontsize=10)
ax.set_xlabel("Freshwater / comparator span ratio (95% CI)")
tag(ax, "b")
fig.tight_layout()
save(fig, "Fig4")

# ---------- Fig 5 ----------
ap = pd.read_csv(PG + "/aa_pgls_results.tsv", sep="\t")
AA = list("ACDEFGHIKLMNPQRSTVWY")
sets = [("ALL_ANY", "1 main", "Main (n = 42)", "o", RED), ("ALL_ANY", "2 comparator", "Euhaline only (n = 31)", "s", BLUE),
        ("ALL_ANY", "3 BUSCO", "BUSCO ≥ 85% (n = 36)", "^", TEAL), ("SC8_OGC", "1 main", "One-to-one orthologues (n = 42)", "D", PURPLE)]
fig = plt.figure(figsize=(15, 5))
gs = fig.add_gridspec(1, 3, width_ratios=[2.2, 1, 1])
ax = fig.add_subplot(gs[0])
x = np.arange(20)
for j, (ds, mod, lab, m, c) in enumerate(sets):
    d = ap[(ap.dataset == ds) & ap.model.str.startswith(mod)].set_index("aa").loc[AA]
    if j == 0:
        assert (d.n == 42).all() and d.q_BH.min() >= 0.097 - 1e-3, d.q_BH.min()
    xo = x + (j - 1.5) * 0.18
    ax.errorbar(xo, d.beta_FW, yerr=[d.beta_FW - d.ci_lo, d.ci_hi - d.beta_FW], fmt=m, ms=5.5, mfc=c, mec="white",
                mew=0.5, ecolor=c, elinewidth=0.9, capsize=0, label=lab, zorder=3)
for i, a in enumerate(AA):
    if a in "GAPQ":
        ax.axvspan(i - 0.45, i + 0.45, color=PPINK, alpha=0.6, lw=0, zorder=0)
ax.axhline(0, color=MD, ls="--", lw=1)
ax.set_xticks(x, AA)
for lbl in ax.get_xticklabels():
    if lbl.get_text() in "GAPQ": lbl.set_fontweight("bold"); lbl.set_color(RED)
ax.set(xlabel="Amino acid", ylabel="β$_{FW}$ (CLR, 95% CI)")
ax.set_ylim(-0.085, 0.115)
ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=10)
tag(ax, "a")

sp = pd.read_csv(W + "/server_data/r1rev/aa_gc/aa_gc_species.tsv", sep="\t")
sp = sp[sp.matches_proteome_count == 1]
cl = pd.read_csv(W + "/analysis/aa/aa_species_clr.tsv", sep="\t")
cl = cl[cl.dataset == "ALL_ANY"][["species", "salinity3", "gcp"]]
d = sp.merge(cl, on="species")
assert len(d) == 10 and (d.salinity3 == "FRESHWATER").sum() == 5
rho = pd.read_csv(PG + "/aa_gc_proxy_cor.tsv", sep="\t")
rho = rho[(rho.set == "S10") & (rho.x == "GC4") & (rho.y == "gcp")].iloc[0]
ax = fig.add_subplot(gs[1])
fit = np.polyfit(d.GC4, d.gcp, 1)
xx = np.linspace(d.GC4.min(), d.GC4.max(), 50)
ax.plot(xx, np.polyval(fit, xx), color=MD, lw=1, ls="-", zorder=1)  # visual guide only; statistic is Spearman
for fw, (m, c, lab) in {True: ("o", RED, "Freshwater"), False: ("o", BLUE, "Comparator")}.items():
    e = d[(d.salinity3 == "FRESHWATER") == fw]
    ax.scatter(e.GC4, e.gcp, marker=m, s=62, facecolor=c, edgecolor="white", lw=0.8, label=lab, zorder=3)
    for _, r in e.iterrows():
        ax.annotate(r.species, (r.GC4, r.gcp), xytext=(4, 3), textcoords="offset points", fontsize=10, color=c)
ax.set(xlabel="GC4", ylabel="GC index ln[GARP/FYMINK]")
ax.text(0.97, 0.05, f"Spearman ρ = {rho.rho:.3f}\nP = {rho.P:.1e}, n = {int(rho.n)}", transform=ax.transAxes,
        ha="right", fontsize=10)
ax.legend(frameon=False, loc="upper left")
tag(ax, "b")

gt = pd.read_csv(PG + "/aa_gc_tests.tsv", sep="\t")
gt = gt[gt.set == "S10"]
ax = fig.add_subplot(gs[2])
for j, (mod, lab, c) in enumerate((("FW", "Without GC4", RED), ("FW+GC4", "With GC4 covariate", BLUE))):
    e = gt[gt.model == mod].set_index("response").loc[["clr_G", "clr_A", "clr_P", "clr_Q"]]
    xo = np.arange(4) + (j - 0.5) * 0.25
    ax.errorbar(xo, e.beta_FW, yerr=[e.beta_FW - e.ci_lo, e.ci_hi - e.beta_FW], fmt="o", ms=7, mfc=c, mec="white",
                mew=0.7, ecolor=c, elinewidth=1.4, capsize=3, label=lab, zorder=3)
ax.axhline(0, color=MD, ls="--", lw=1)
ax.set_xticks(range(4), ["Gly", "Ala", "Pro", "Gln"])
ax.set(ylabel="β$_{FW}$ (CLR, 95% CI)", xlabel="10 species (5 FW)")
ax.set_ylim(-0.105, 0.13)
ax.legend(frameon=False, loc="upper right")
tag(ax, "c")
fig.tight_layout()
save(fig, "Fig5")
print("done", sorted(os.listdir(OUT)))
