"""Fig. 3b (R1 revision, 2026-11-01): gene trees of the Pfam-defined ATG8/LC3/GABARAP, Beclin 1 and
aquaporin families, drawn from the IQ-TREE 2.4.0 maximum-likelihood trees (<family>.treefile; MAFFT
v7.490 L-INS-i alignments, ModelFinder BIC models, 1,000 ultrafast bootstrap replicates; run_fig3b.sh).
One tip per gene (the longest protein); tip habitat from fig3b_tip_labels.tsv (the 40 counted bivalves
plus Lottia gigantea copies). Unrooted equal-angle layout (Felsenstein 2004) started at the tip centroid.
Edge rule, independent of rooting: an edge is drawn red when every tip on one side of it is a freshwater
gene copy. Black dots on branches mark splits with UFBoot >= 95; branches that only join copies of one
identical sequence (set aside by IQ-TREE during the search and re-attached at zero length; <family>.log)
carry no bootstrap value and no dot.
Outputs (Figures_R1/panels_for_Fig3b): Fig3b_trees_R1.{pdf,png}, Fig3b_source_data.tsv,
Fig3b_freshwater_clades.tsv, Fig3b_tree_summary.tsv.
"""
import math, pathlib, re, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

sys.setrecursionlimit(20000)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic"})
ROOT = pathlib.Path(__file__).resolve().parents[3]
TREES = ROOT / "_R1_work/server_data/r1rev/fig3b_trees"
OUT = ROOT / "SUBMIT_R1_revised_20261024/Figures_R1/panels_for_Fig3b"
HAB_FW, HAB_SW, OUTG, EDGE = "#B5222F", "#286EAE", "#888888", "#A9C6E3"
W, H = 14.6, 5.6
FAMS = [("ATG8", "ATG8/LC3/GABARAP", 3.6, 16), ("BECN1", "Beclin 1", 3.0, 40), ("AQP", "Aquaporin/AQP", 4.6, 6)]
UNIT = {"Pst": "Unionida", "Sol": "Unionida", "Awo": "Unionida", "Ude": "Unionida",
        "Dpo": "Dreissenidae", "Dro": "Dreissenidae", "Cfl": "Corbicula", "Lfo": "Limnoperna"}
key = lambda a, b: (a, b) if a < b else (b, a)


def parse_newick(s):
    """Minimal Newick reader -> parent, children, label, branch length per node (node 0 = root)."""
    s = s.strip().rstrip(";")
    par, kids, lab, bl = [-1], [[]], [""], [0.0]
    cur, i, n = 0, 0, len(s)
    while i < n:
        c = s[i]
        if c == "(" or c == ",":
            p = cur if c == "(" else par[cur]
            par.append(p); kids.append([]); lab.append(""); bl.append(0.0)
            cur = len(par) - 1; kids[p].append(cur); i += 1
        elif c == ")":
            cur = par[cur]; i += 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in ",();":
                j += 1
            bl[cur] = float(s[i + 1:j]); i = j
        else:
            j = i
            while j < n and s[j] not in ",():;":
                j += 1
            lab[cur] = s[i:j].strip(); i = j
    return par, kids, lab, bl


def orient(nb, start):
    """Re-root the unrooted tree at `start`; children keep the cyclic (Newick) order."""
    par, ch, order, k = {start: -1}, {}, [start], 0
    while k < len(order):
        v = order[k]; k += 1
        lst = nb[v]
        if par[v] < 0:
            c = list(lst)
        else:
            i = lst.index(par[v]); c = lst[i + 1:] + lst[:i]
        ch[v] = c
        for w in c:
            par[w] = v; order.append(w)
    return par, ch, order
def counts(ch, order, isleaf, isfw):
    tip, nfw = {}, {}
    for v in reversed(order):
        if isleaf[v]:
            tip[v], nfw[v] = 1, int(isfw[v])
        else:
            tip[v] = sum(tip[w] for w in ch[v]); nfw[v] = sum(nfw[w] for w in ch[v])
    return tip, nfw


def layout(ch, order, tip, blen, start):
    """Equal-angle unrooted layout: each subtree gets a wedge proportional to its tip count."""
    xy, wedge = {start: (0.0, 0.0)}, {start: (0.0, 2 * math.pi)}
    for v in order:
        a0, a1 = wedge[v]
        tot, a = sum(tip[w] for w in ch[v]), a0
        for w in ch[v]:
            span = (a1 - a0) * tip[w] / tot
            wedge[w] = (a, a + span); mid = a + span / 2; l = blen[key(v, w)]
            xy[w] = (xy[v][0] + l * math.cos(mid), xy[v][1] + l * math.sin(mid))
            a += span
    return xy


def centroid(nb, isleaf):
    """Internal node whose removal leaves the most even split of the tips."""
    par, ch, order = orient(nb, 0)
    tip, _ = counts(ch, order, isleaf, isleaf)
    T, best, bv = tip[0], None, None
    for v in order:
        if not isleaf[v]:
            m = max([tip[w] for w in ch[v]] + ([T - tip[v]] if par[v] >= 0 else []))
            if best is None or m < best:
                best, bv = m, v
    return bv


def fit(xy, bw, bh):
    """Rotation (5-degree steps) that lets the drawing fill a bw x bh box at the largest scale."""
    pts, best = list(xy.values()), None
    for d in range(0, 180, 5):
        c, s = math.cos(math.radians(d)), math.sin(math.radians(d))
        xs = [x * c - y * s for x, y in pts]; ys = [x * s + y * c for x, y in pts]
        sc = min(bw / (max(xs) - min(xs)), bh / (max(ys) - min(ys)))
        if best is None or sc > best[0]:
            best = (sc, d)
    c, s = math.cos(math.radians(best[1])), math.sin(math.radians(best[1]))
    return {v: (x * c - y * s, x * s + y * c) for v, (x, y) in xy.items()}, best[1]
def leaves_below(ch, v):
    out, st = [], [v]
    while st:
        u = st.pop()
        if ch[u]:
            st.extend(ch[u])
        else:
            out.append(u)
    return out


tl = pd.read_csv(OUT / "fig3b_tip_labels.tsv", sep="\t")
fig = plt.figure(figsize=(W, H))
BOT, TOP, X0 = 1.15, 0.35, 0.35          # inches: label strip below the trees, top margin, left margin
src_rows, clade_rows, summ, x_in = [], [], [], X0
for fam, title, bw, msize in FAMS:
    par0, kids0, lab, bl0 = parse_newick((TREES / f"{fam}.treefile").read_text())
    N = len(par0)
    nb = {v: ([] if par0[v] < 0 else [par0[v]]) + kids0[v] for v in range(N)}
    blen, sup = {}, {}
    for v in range(1, N):
        blen[key(v, par0[v])] = max(bl0[v], 0.0)
        if kids0[v] and lab[v]:
            sup[key(v, par0[v])] = float(lab[v])
    isleaf = {v: not kids0[v] for v in range(N)}
    t = tl[tl.family == fam].set_index("label")
    hab = {v: t.at[lab[v], "habitat"] for v in range(N) if isleaf[v]}     # KeyError if a tip is unmatched
    spp = {v: t.at[lab[v], "species"] for v in range(N) if isleaf[v]}
    assert len(hab) == len(t) and t.index.is_unique, (fam, len(hab), len(t))
    lg, ident = (TREES / f"{fam}.log").read_text(), {}       # IQ-TREE identical-sequence groups (log NOTEs)
    for a, b in re.findall(r"NOTE: (\S+) is identical to (\S+) but kept", lg) + \
            re.findall(r"NOTE: (\S+) \(identical to (\S+)\) is ignored but added at the end", lg):
        ident.setdefault(b, {b}).add(a)
    grp = {x: r for r, g in ident.items() for x in g}
    hab_l = {lab[v]: hab[v] for v in hab}
    mixed = sorted(r for r, g in ident.items() if len({hab_l[x] for x in g}) > 1)
    if mixed:
        print(fam, "identical-sequence groups with mixed habitats:", mixed)
    isfw = {v: isleaf[v] and hab[v] == "freshwater" for v in range(N)}
    c0 = centroid(nb, isleaf)
    par, ch, order = orient(nb, c0)
    tip, nfw = counts(ch, order, isleaf, isfw)
    T, F = tip[c0], nfw[c0]
    for v in order[1:]:
        assert T - tip[v] != F - nfw[v], (fam, v)          # the side away from v is never all-freshwater
    red = {v: v != c0 and nfw[v] == tip[v] for v in order}  # edge par[v]-v: every tip below v freshwater
    bh = H - BOT - TOP
    xy, rot = fit(layout(ch, order, tip, blen, c0), bw, bh)
    xs = [p[0] for p in xy.values()]; ys = [p[1] for p in xy.values()]
    sc = min(bw / (max(xs) - min(xs)), bh / (max(ys) - min(ys))) * 0.94
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    ax = fig.add_axes([x_in / W, BOT / H, bw / W, bh / H])
    ax.set_xlim(cx - bw / 2 / sc, cx + bw / 2 / sc); ax.set_ylim(cy - bh / 2 / sc, cy + bh / 2 / sc)
    ax.set_aspect("equal"); ax.axis("off")
    assert set(hab.values()) <= {"freshwater", "brackish", "euhaline", "outgroup"}, set(hab.values())
    unl = [v for v in order[1:] if not isleaf[v] and key(v, par[v]) not in sup]
    for v in unl:                       # branches without UFBoot must only join copies of one identical sequence
        g = {grp.get(lab[u], lab[u]) for u in leaves_below(ch, v)}
        assert len(g) == 1 and g <= set(ident), (fam, g)
    seg = {True: [], False: []}
    for v in order[1:]:
        seg[red[v]].append([xy[par[v]], xy[v]])
    ax.add_collection(LineCollection(seg[False], colors=EDGE, linewidths=0.6, zorder=1, capstyle="round"))
    ax.add_collection(LineCollection(seg[True], colors=HAB_FW, linewidths=1.6, zorder=2, capstyle="round"))
    dots = [((xy[par[v]][0] + xy[v][0]) / 2, (xy[par[v]][1] + xy[v][1]) / 2)
            for v in order[1:] if not isleaf[v] and sup.get(key(v, par[v]), 0) >= 95]
    ax.scatter([d[0] for d in dots], [d[1] for d in dots], s=4, c="#222222", linewidths=0, zorder=3)
    leaves = [v for v in order if isleaf[v]]
    for h, col, z in (("euhaline", HAB_SW, 4), ("brackish", HAB_SW, 4), ("outgroup", OUTG, 4), ("freshwater", HAB_FW, 5)):
        L = [v for v in leaves if hab[v] == h]
        ax.scatter([xy[v][0] for v in L], [xy[v][1] for v in L], s=msize, c=col, edgecolors="white",
                   linewidths=0.3, zorder=z)
    raw = (max(xs) - min(xs)) * 0.2; e = math.floor(math.log10(raw)); m = raw / 10 ** e
    sl = (1 if m < 2 else 2 if m < 5 else 5) * 10 ** e
    x_l, yb = cx - bw / 2 / sc + 0.12 / sc, (BOT - 0.1) / H
    tr = blended_transform_factory(ax.transData, fig.transFigure)
    ax.add_line(Line2D([x_l, x_l + sl], [yb, yb], color="black", lw=1.0, transform=tr, clip_on=False))
    ax.text(x_l + sl + 0.06 / sc, yb, f"{sl:g} substitutions per site", transform=tr, ha="left", va="center", fontsize=10)
    xm = (x_in + bw / 2) / W
    nf = sum(hab[v] == "freshwater" for v in leaves); no = sum(hab[v] == "outgroup" for v in leaves); ns = T - nf - no
    fig.text(xm, 0.64 / H, title, ha="center", va="center", fontsize=14, fontweight="bold")
    fig.text(xm, 0.33 / H, f"{nf} freshwater + {ns} seawater genes; {no} $\\it{{Lottia}}$", ha="center", va="center", fontsize=10)
    cl_of, k = {}, 0
    for v in order:
        if red[v] and not red.get(par[v], False):
            k += 1; lv = leaves_below(ch, v); cid = f"{fam}_FW{k:03d}"
            s_ = sorted({spp[u] for u in lv}); un = sorted({UNIT[s] for s in s_})
            cl_of.update({u: cid for u in lv})
            clade_rows.append(dict(family=fam, clade=cid, n_genes=len(lv), n_species=len(s_), species=",".join(s_),
                                   n_units=len(un), units=",".join(un), UFBoot=sup.get(key(v, par[v])),
                                   branch="terminal" if isleaf[v] else "internal" if key(v, par[v]) in sup
                                   else "identical sequences (no support value)",
                                   labels=",".join(sorted(lab[u] for u in lv))))
    for v in leaves:
        src_rows.append(dict(family=fam, label=lab[v], species=spp[v], habitat=hab[v], unit=UNIT.get(spp[v], ""),
                             x=round(xy[v][0], 6), y=round(xy[v][1], 6), terminal_branch_red=red[v],
                             freshwater_clade=cl_of.get(v, ""), identical_to=grp.get(lab[v], "")))
    cr = [r for r in clade_rows if r["family"] == fam]; mu = [r for r in cr if r["n_genes"] > 1]
    print(f"{fam}: tips {T} (FW {nf}, SW {ns}, Lottia {no}); rotated {rot} deg; internal branches "
          f"{sum(not isleaf[v] for v in order[1:])} ({len(unl)} joining identical sequences), UFBoot>=95 "
          f"{len(dots)}; FW clades {len(cr)}: singletons "
          f"{len(cr) - len(mu)}, >=2 genes {len(mu)} (UFBoot>=95 {sum((r['UFBoot'] or 0) >= 95 for r in mu)}, "
          f">1 unit {sum(r['n_units'] > 1 for r in mu)}, largest {max([r['n_genes'] for r in cr], default=0)})")
    q = (TREES / f"{fam}.iqtree").read_text()
    mb = re.search(r"Best-fit model according to BIC: (\S+)", q)
    md = re.search(r"Input data: (\d+) sequences with (\d+) amino-acid sites", q)
    assert mb and md and int(md.group(1)) == T, fam
    summ.append(dict(family=fam, tips=T, freshwater=nf, seawater=ns, lottia=no, BIC_model=mb.group(1),
                     aligned_sites=int(md.group(2)), internal_branches=sum(not isleaf[v] for v in order[1:]),
                     UFBoot95=len(dots), fw_clades=len(cr), fw_singletons=len(cr) - len(mu), fw_multi=len(mu),
                     fw_multi_UFBoot95=sum((r["UFBoot"] or 0) >= 95 for r in mu),
                     identical_groups=len(ident), identical_seqs=sum(len(g) for g in ident.values()),
                     branches_identical_no_support=len(unl), identical_mixed_habitat=",".join(mixed),
                     fw_multi_gt1_unit=sum(r["n_units"] > 1 for r in mu),
                     largest_fw_clade=max([r["n_genes"] for r in cr], default=0), rotation_deg=rot, scale_bar=sl))
    x_in += bw + 0.25

hab_h = [Line2D([], [], marker="o", ls="none", mfc=c, mec="white", mew=0.3, ms=7, label=l) for c, l in
         ((HAB_FW, "Freshwater"), (HAB_SW, "Seawater\n(brackish + euhaline)"),
          (OUTG, "$\\it{Lottia\\ gigantea}$\n(outgroup)"))]
br_h = [Line2D([], [], color=HAB_FW, lw=1.6, label="Branch leading only to\nfreshwater genes"),
        Line2D([], [], color=EDGE, lw=0.6, label="Other branch"),
        Line2D([], [], marker="o", ls="none", mfc="#222222", mec="none", ms=2.5, label="UFBoot ≥ 95")]
lx = (x_in + 0.05) / W
l1 = fig.legend(handles=hab_h, loc="upper left", bbox_to_anchor=(lx, (H - TOP) / H), frameon=False, fontsize=10,
                title="Habitat", title_fontsize=10, alignment="left", handlelength=1.4, labelspacing=0.7)
l1.get_title().set_fontweight("bold")
fig.canvas.draw()
y2 = l1.get_window_extent().transformed(fig.transFigure.inverted()).y0 - 0.15 / H
l2 = fig.legend(handles=br_h, loc="upper left", bbox_to_anchor=(lx, y2), frameon=False, fontsize=10,
                title="Branches", title_fontsize=10, alignment="left", handlelength=1.4, labelspacing=0.7)
l2.get_title().set_fontweight("bold")
fig.text(0.1 / W, (H - 0.1) / H, "b", fontsize=14, fontweight="bold", va="top")
fig.savefig(OUT / "Fig3b_trees_R1.pdf"); fig.savefig(OUT / "Fig3b_trees_R1.png", dpi=300)
pd.DataFrame(src_rows).to_csv(OUT / "Fig3b_source_data.tsv", sep="\t", index=False)
pd.DataFrame(clade_rows).to_csv(OUT / "Fig3b_freshwater_clades.tsv", sep="\t", index=False)
pd.DataFrame(summ).to_csv(OUT / "Fig3b_tree_summary.tsv", sep="\t", index=False)
print("wrote", OUT)
