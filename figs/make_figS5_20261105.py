# R1 revision (2026-11-05), audit item (4) / R2-M5: per-locus DESeq2 contrasts for the Fig. 6f genes.
#   Fig. S5   dot panel, 17 RefSeq loci of the nine genes in Fig. 6f x six treatment-vs-NS contrasts;
#             dot area = -log10(BH q), colour = log2 fold change, black outline = q < 0.05.
#   Table S9  the same statistics for every locus and contrast, read counts of the 17 loci, sample list.
# Input: server ~/R1_revision_20261024/fig6f_de (salmon 1.4.0 on GCF_014805675.1 RefSeq transcripts,
# gene-level sums, DESeq2 1.42.1 ~ trt, NS reference, filter >= 10 reads in >= 3 libraries, BH within
# contrast over all tested genes); copied to SUBMIT_R1_revised_20261024/_R1_work/fig6f_de (sha256 checked).
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

R = (os.environ.get("IZ_ROOT", ".") + "/SUBMIT_R1_revised_20261024")
IN = R + "/_R1_work/fig6f_de"
OUT = R + "/Figures_R1/supplementary"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.labelsize": 12, "pdf.fonttype": 42,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.linewidth": 0.8,
                     "axes.edgecolor": "#222222", "axes.unicode_minus": True})
MINUS = "\u2212"

LOCI = [("ATG7", "LOC123529011"), ("ATG5", "LOC123552033"), ("ATG5", "LOC123552930"),
        ("BECN1", "LOC123538342"), ("SQSTM1", "LOC123556581"), ("SQSTM1", "LOC123556582"),
        ("TFEB", "LOC123528641"), ("PRKN", "LOC123552401"), ("ATG4", "LOC123554457"),
        ("ATG4", "LOC123524438"), ("ATG9", "LOC123532473"),
        ("DRAM1-like", "LOC123536557"), ("DRAM1-like", "LOC123536558"), ("DRAM1-like", "LOC123536931"),
        ("DRAM2-like", "LOC123536158"), ("DRAM2-like", "LOC123536560"), ("DRAM2-like", "LOC123537516")]
CON = ["AH1", "AH5", "CH", "AL1", "AL5", "CL"]
CLAB = {"AH1": "Acute high\nday 1", "AH5": "Acute high\nday 5", "CH": "Chronic high\nday 50",
        "AL1": "Acute low\nday 1", "AL5": "Acute low\nday 5", "CL": "Chronic low\nday 50"}

de = pd.read_csv(IN + "/deseq2_all_vs_NS.tsv", sep="\t")
cnt = pd.read_csv(IN + "/counts17.tsv", sep="\t", index_col=0)
smp = pd.read_csv(IN + "/samples.tsv", sep="\t", header=None, names=["run", "treatment"])
gp = pd.read_csv(IN + "/gene_product.tsv", sep="\t", header=None, names=["gene", "product"])
gp["gene"] = gp["gene"].str.replace("gene=", "", regex=False)
gp["product"] = gp["product"].str.replace("%2C", ",", regex=False).str.replace(r",? transcript variant X\d+", "", regex=True)
prod = gp.drop_duplicates("gene").set_index("gene")["product"]
n_tested = de["gene"].nunique()
assert len(de) == 6 * n_tested

rows = []
for g, loc in LOCI:
    c = cnt.loc[loc, smp["run"]].astype(int)
    ok = int((c >= 10).sum())
    for t in CON:
        r = de[(de["gene"] == loc) & (de["contrast"] == t + "_vs_NS")]
        base = dict(gene=g, locus=loc, product=prod.get(loc, ""), contrast=t + " vs NS",
                    libraries_with_10_reads=ok, total_reads=int(c.sum()))
        if len(r) == 1:
            r = r.iloc[0]
            base.update(tested="yes", baseMean=r.baseMean, log2FC=r.log2FoldChange, lfcSE=r.lfcSE,
                        Wald_stat=r.stat, P=r.pvalue, BH_q=r.padj)
        else:
            assert len(r) == 0
            base.update(tested="no (fewer than 3 libraries with >= 10 reads)", baseMean=np.nan, log2FC=np.nan,
                        lfcSE=np.nan, Wald_stat=np.nan, P=np.nan, BH_q=np.nan)
        rows.append(base)
T = pd.DataFrame(rows)

# ---------- figure
nl = len(LOCI)
fig, ax = plt.subplots(figsize=(8.6, 8.2))
lim = np.nanmax(np.abs(T["log2FC"]))
lim = np.ceil(lim * 2) / 2
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
cmap = plt.get_cmap("RdBu_r")
SZ = lambda q: 18 + 70 * np.clip(-np.log10(q), 0, 6)
for i, (g, loc) in enumerate(LOCI):
    y = nl - 1 - i
    for j, t in enumerate(CON):
        r = T[(T.locus == loc) & (T.contrast == t + " vs NS")].iloc[0]
        if r.tested != "yes":
            ax.scatter(j, y, marker="x", s=22, color="#888888", lw=0.9, zorder=3)
        elif np.isnan(r.BH_q):
            ax.scatter(j, y, s=18, facecolor="none", edgecolor="#888888", lw=0.8, zorder=3)
        else:
            ax.scatter(j, y, s=SZ(r.BH_q), color=cmap(norm(r.log2FC)), zorder=3,
                       edgecolor="#000000" if r.BH_q < 0.05 else "none", lw=1.0)
ax.set_xlim(-0.6, len(CON) - 0.4); ax.set_ylim(-0.7, nl - 0.3)
ax.set_xticks(range(len(CON))); ax.set_xticklabels([CLAB[t] for t in CON], fontsize=10)
ax.xaxis.tick_top()
ax.set_yticks(range(nl))
ax.set_yticklabels([f"{g}  {loc}" for g, loc in LOCI][::-1], fontsize=10)
for y in range(nl):
    ax.axhline(y, color="#eeeeee", lw=0.6, zorder=0)
ax.axvline(2.5, color="#888888", lw=0.8, ls=(0, (3, 2)))
# separate genes
prev = None
for i, (g, _) in enumerate(LOCI):
    if prev is not None and g != prev:
        ax.axhline(nl - 1 - i + 0.5, color="#bbbbbb", lw=0.6)
    prev = g
for s in ("top", "right", "bottom", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(length=0)
ax.set_xlabel("Treatment vs normal seawater (NS)", fontsize=12, labelpad=8)
ax.xaxis.set_label_position("top")

sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cax = fig.add_axes([0.93, 0.55, 0.022, 0.25])
cb = fig.colorbar(sm, cax=cax)
cb.set_label("log$_2$ fold change", fontsize=11)
cb.ax.tick_params(labelsize=10)
qs = [0.5, 0.05, 0.001, 1e-6]
h = [Line2D([], [], ls="", marker="o", markersize=np.sqrt(SZ(q)), markerfacecolor="#999999",
            markeredgecolor="none") for q in qs]
h.append(Line2D([], [], ls="", marker="o", markersize=np.sqrt(SZ(0.01)), markerfacecolor="white",
                markeredgecolor="#000000"))
h.append(Line2D([], [], ls="", marker="o", markersize=np.sqrt(18), markerfacecolor="none",
                markeredgecolor="#888888"))
h.append(Line2D([], [], ls="", marker="x", markersize=5, color="#888888"))
lab = ["q = 0.5", "q = 0.05", "q = 0.001", "q ≤ 10$^{" + MINUS + "6}$", "outline: q < 0.05",
       "q not computed\n(DESeq2 independent\nfiltering)", "not tested\n(< 3 libraries\nwith ≥ 10 reads)"]
fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.905, 0.53), frameon=False, fontsize=10,
           title="BH-adjusted q", title_fontsize=10, labelspacing=0.9, handletextpad=0.8)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/FigS5_R1.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------- source data and Table S9
T.to_csv(f"{OUT}/FigS5_source_data.tsv", sep="\t", index=False, float_format="%.6g")
C = cnt.loc[[l for _, l in LOCI], smp["run"]].copy()
C.insert(0, "gene", [g for g, _ in LOCI]); C.index.name = "locus"
S = smp.copy()
readme = pd.DataFrame({"Table S9": [
    "Per-locus differential expression of the nine autophagy-related genes shown in Fig. 6f (BioProject PRJNA975821, 21 libraries, three per treatment).",
    "Reads were quantified with salmon 1.4.0 (selective alignment, automatic library type; IU detected in all libraries) against the RefSeq transcripts of GCF_014805675.1 and summed to RefSeq genes.",
    "DESeq2 1.42.1, design ~ treatment, normal seawater (NS) as reference; genes with >= 10 reads in >= 3 libraries were tested (%d genes); Wald tests; Benjamini-Hochberg q within each contrast over all tested genes." % n_tested,
    "Treatments: AH1/AH5/CH = acute high salinity day 1 / day 5 / chronic high day 50; AL1/AL5/CL = acute low salinity day 1 / day 5 / chronic low day 50.",
    "All libraries come from one tank per treatment, so q values describe differences between treatment groups in this experiment and are exploratory.",
    "Loci with fewer than three libraries at >= 10 reads were not tested; their read counts are given in sheet 'read_counts'.",
]})
with pd.ExcelWriter(R + "/Supplementary_Table_S9_fig6f_locus_DESeq2.xlsx") as xw:
    readme.to_excel(xw, sheet_name="README", index=False)
    T.to_excel(xw, sheet_name="DESeq2_vs_NS", index=False)
    C.to_excel(xw, sheet_name="read_counts")
    S.to_excel(xw, sheet_name="samples", index=False)
T.to_csv(R + "/Supplementary_Table_S9_fig6f_locus_DESeq2.tsv", sep="\t", index=False, float_format="%.6g")

print("tested genes:", n_tested)
pd.set_option("display.width", 220)
print(T[["gene", "locus", "contrast", "libraries_with_10_reads", "log2FC", "BH_q"]].to_string())
print("q<0.05:", (T.BH_q < 0.05).sum(), "of", T.BH_q.notna().sum())
