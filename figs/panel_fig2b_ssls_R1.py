"""Fig. 2b (R1): dSSLS from the real server run (tables written by ../ssls/ssls_stats.py).
R0 raincloud layout kept: jittered points left, half-violin right, IQR bar with median; P and BH q on the bracket."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
ST = os.path.join(HERE, "..", "ssls", "stats")
PKG = os.path.join(HERE, "..", "..", "..", "SUBMIT_R1_revised_20261024", "Figures_R1", "panels_for_Fig2")
QA = os.path.join(HERE, "..", "..", "qa_figs_20261031")
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "font.weight": "bold", "axes.labelweight": "bold",
                     "axes.linewidth": 1.0, "pdf.fonttype": 42, "svg.fonttype": "none"})
G = pd.read_csv(os.path.join(ST, "ssls_gene_table.tsv"), sep="\t")
S = pd.read_csv(os.path.join(ST, "ssls_site_table.tsv"), sep="\t")
T = pd.read_csv(os.path.join(ST, "ssls_tests.tsv"), sep="\t")
pick = lambda lv: T[(T.level == lv) & (T.subset == "recurrent") & (T.role == "primary")].iloc[0]
tg, ts = pick("gene"), pick("site")
gr = G[G.rand1000]; gc = G[(G.set == "cand") & (G.status == "ok") & G.has_recurrent]
sub = S[S.og.isin(set(S.og[S.is_recurrent]))]
so, sc = sub[~sub.is_recurrent], sub[sub.is_recurrent]
assert (len(gr), len(gc), len(so), len(sc)) == (tg.n_comp, tg.n_focal, ts.n_comp, ts.n_focal)
W = 0.5                                   # site panel: inverse hyperbolic sine axis, linear within about +-0.5
fwd = lambda x: np.arcsinh(np.asarray(x, float) / W)
MINUS = "−"

def fmt(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return rf"{m}$\times$10$^{{{int(e)}}}$"

def rain(ax, x, v, fill, edge, size, rng):
    d = gaussian_kde(v); y = np.linspace(v.min(), v.max(), 400); w = d(y); w = 0.34 * w / w.max()
    ax.fill_betweenx(y, x + 0.03, x + 0.03 + w, color=fill, lw=0, zorder=1)
    ax.plot(x + 0.03 + w, y, color=edge, lw=0.8, zorder=2)
    q1, md, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([x + 0.09] * 2, [q1, q3], color="black", lw=3.4, solid_capstyle="butt", zorder=3)
    ax.scatter([x + 0.09], [md], s=18, color="white", edgecolor="black", lw=0.8, zorder=4)
    ax.scatter(x - 0.06 - rng.uniform(0, 0.24, len(v)), v, s=size, color=edge, alpha=0.45, lw=0, rasterized=True, zorder=2)

def bracket(ax, y, h, t):
    ax.plot([0.09, 0.09, 1.09, 1.09], [y, y + h, y + h, y], color="black", lw=1.0)
    ax.text(0.59, y + 1.6 * h, f"P = {fmt(t.P_MW)}\nq = {fmt(t.q_MW)}", ha="center", va="bottom", linespacing=1.15)

rng = np.random.default_rng(20261101)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(5.8, 3.3), gridspec_kw={"wspace": 0.45})
rain(a1, 0, gr.dssls_mean.values, "#f2a7a5", "#ca2a28", 5, rng)
rain(a1, 1, gc.dssls_mean.values, "#f2a7a5", "#ca2a28", 14, rng)
top = max(gr.dssls_mean.max(), gc.dssls_mean.max()); bracket(a1, top + 0.10, 0.08, tg)
a1.set_ylim(min(gr.dssls_mean.min(), gc.dssls_mean.min()) - 0.15, top + 0.62)
a1.set_ylabel("Gene-level ΔSSLS (mean per site)")
a1.set_xticks([0, 1], [f"Random\ngenes\nn = {len(gr):,}", f"Genes with\nrecurrent\nsubstitutions\nn = {len(gc)}"])
rain(a2, 0, fwd(so.dssls), "#8fb6dd", "#2371a9", 3, rng)
rain(a2, 1, fwd(sc.dssls), "#8fb6dd", "#2371a9", 5, rng)
top = max(fwd(so.dssls).max(), fwd(sc.dssls).max()); bracket(a2, top + 0.25, 0.20, ts)
a2.set_ylim(min(fwd(so.dssls).min(), fwd(sc.dssls).min()) - 0.4, top + 1.9)
tk = [-5, -2, -1, 0, 1, 2, 5, 10]
a2.set_yticks(fwd(tk), [f"{MINUS}{-t}" if t < 0 else f"{t}" for t in tk])
a2.set_ylabel("Site-level ΔSSLS (asinh scale)")
a2.set_xticks([0, 1], [f"Other\nsites\nn = {len(so):,}", f"Recurrent-\nsubstitution\nsites\nn = {len(sc)}"])
for ax in (a1, a2):
    ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-0.42, 1.52); ax.tick_params(axis="x", length=0, pad=4)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(PKG, f"Fig2b_ssls_R1.{ext}"), dpi=300, bbox_inches="tight", pad_inches=0.04)
fig.savefig(os.path.join(QA, "Fig2b_ssls_R1_qa.png"), dpi=250, bbox_inches="tight", pad_inches=0.04)
src = pd.concat([pd.DataFrame({"level": "gene", "group": "random_gene", "og": gr.og, "site": "", "n_taxa": gr.n_taxa,
                               "n_fw": gr.n_fw, "dssls": gr.dssls_mean}),
                 pd.DataFrame({"level": "gene", "group": "recurrent_substitution_gene", "og": gc.og, "site": "",
                               "n_taxa": gc.n_taxa, "n_fw": gc.n_fw, "dssls": gc.dssls_mean}),
                 pd.DataFrame({"level": "site", "group": "other_site", "og": so.og, "site": so.site, "dssls": so.dssls}),
                 pd.DataFrame({"level": "site", "group": "recurrent_substitution_site", "og": sc.og, "site": sc.site,
                               "dssls": sc.dssls})])
src.to_csv(os.path.join(PKG, "Fig2b_source_data.tsv"), sep="\t", index=False)
print({k: len(v) for k, v in src.groupby("group")}, "| gene P/q", fmt(tg.P_MW), fmt(tg.q_MW), "| site P/q", fmt(ts.P_MW), fmt(ts.q_MW))
