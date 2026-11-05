# Stand-alone replacement panels for the author's R0 Fig. 2 (author pastes them in; the main figure is not redrawn).
# panel_a_gene_venn      -> replaces R0 Fig. 2a (old Fsite 1,931 / 894 / 263 / 109 Venn)
# panel_a_pathway_venn   -> alternative/companion for Fig. 2a at KEGG pathway level
# Input: bsite/three_lines.py output dir (venn.json). Output: Figures_R1/panels_for_Fig2/*.png (300 dpi) + .pdf
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from venn3 import venn3

R = os.environ.get("IZ_ROOT", ".")
TL = os.environ.get("TL_DIR", R + "/_R1_work/analysis/bsite/final_3l")
OUT = os.environ.get("FIG_OUT", R + "/SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig2")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42})
RED, BLUE, TEAL = "#ca2a28", "#2371a9", "#00a79d"
vj = json.load(open(TL + "/venn.json"))
tot = vj["totals"]
LAB = ["Positive selection", "Parallel/convergent", "Expanded families"]


def one(name, counts, totals, title):
    fig, ax = plt.subplots(figsize=(3.4, 3.4))
    venn3(ax, counts, LAB, [RED, BLUE, TEAL], totals=totals, fs=9.5)
    ax.set_title(title, fontsize=10.5)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


one("panel_a_gene_venn", vj["gene_main"], [tot["PSG"], tot["PCG"], tot["EXP_hog"]],
    f"Genes (tested orthogroups, n = {tot['tested']:,})")
one("panel_a_pathway_venn", vj["pathway_represented"], [tot["pw_PSG"], tot["pw_PCG"], tot["pw_EXP"]],
    "KEGG pathways")
print("written", OUT, tot)
