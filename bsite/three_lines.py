"""R1 revision, Fig. 2a: three independent lines of evidence for freshwater adaptation, following the logic of
Yuan et al. (2020, Mol Ecol 29:3954-3969): positively selected genes (PSGs, branch-site model A on the freshwater
unit stems), parallel/convergent substitution genes (PCGs, ancestral reconstruction) and gene families significantly
expanded on freshwater unit stems (CAFE5), intersected at gene/family level and at KEGG pathway level.

usage: python3 three_lines.py <bsite_run_dir> <kegg_dir> <fw_expansion.tsv> <pathway_list.tsv> <out_dir> <pc_control.tsv>
Definitions:
  PSG  : status ok, LRT = max(0, 2(lnL_alt - lnL_null)), p from chi2(1), BH q < 0.05 over all tested HOGs.
  PCG  : >=1 site where >=3 freshwater unit stems change (pp >= 0.8) to the same derived residue absent from all
         comparators (pc_control.py K=3 PP=0.8, with marine 4-branch control rate). Revised after the >=2-unit rule
         turned out to flag ~97% of tested HOGs; the >=2-unit sets are reported as sensitivity.
  EXP  : orthogroup with family p < 0.05 and a significant (branch p < 0.05) gain on >=1 freshwater unit stem
         (main); strict = on >=2 units. Gene-level overlap uses universe HOGs whose parent OG is expanded.
  KEGG : consensus pathways (>= 50% of annotated members), global maps (map011xx/012xx) excluded. Enrichment by
         one-sided Fisher test against the set's own background (PSG/PCG: annotated universe HOGs; EXP: annotated
         orthogroups tested by CAFE5), BH within set.
"""
import json, math, os, sys
from collections import defaultdict

RUN, KD, FWE, PL, OUT, PCC = sys.argv[1:7]
os.makedirs(OUT, exist_ok=True)


def chi2sf1(x):
    return math.erfc(math.sqrt(x / 2)) if x > 0 else 1.0


def bh(ps):
    n = len(ps)
    o = sorted(range(n), key=lambda i: ps[i])
    q = [0.0] * n
    m = 1.0
    for r in range(n - 1, -1, -1):
        i = o[r]
        m = min(m, ps[i] * n / (r + 1))
        q[i] = m
    return q


def lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_gt(a, K, n, N):
    """P(X >= a), X ~ Hypergeometric(N population, K successes, n draws)."""
    hi = min(K, n)
    return min(1.0, sum(math.exp(lchoose(K, x) + lchoose(N - K, n - x) - lchoose(N, n)) for x in range(a, hi + 1)))


def tsv(p):
    f = open(p)
    h = f.readline().rstrip("\n").split("\t")
    return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in f]


names = dict(l.rstrip("\n").split("\t")[:2] for l in open(PL) if "\t" in l)
glob = lambda t: t.startswith("map011") or t.startswith("map012")

# ---- PSG
res = [r for r in tsv(RUN + "/bsite_results.tsv") if r["status"] == "ok"]
for r in res:
    r["lrt"] = max(0.0, 2 * (float(r["lnL_alt"]) - float(r["lnL_null"])))
    r["p"] = chi2sf1(r["lrt"])
for r, q in zip(res, bh([r["p"] for r in res])):
    r["q"] = q
tested = {r["og"] for r in res}
all4 = {r["og"] for r in res if r["units"] == "4"}
PSG = {r["og"] for r in res if r["q"] < 0.05}
PSGnom = {r["og"] for r in res if r["p"] < 0.05}

# ---- PCG (restricted to tested HOGs so that the three lines share one universe)
# main: >=1 site where >=3 freshwater unit stems change (pp >= 0.8) to the same residue absent from all comparators
# (pc_control.py, K=3, PP=0.8; same rule scored on marine 4-branch control sets). Sensitivity: >=2 units (any event;
# and minpp >= 0.8 + comparator-absent).
ev = [e for e in tsv(RUN + "/pc_events.tsv") if e["og"] in tested]
ctl = {r["hog"]: r for r in tsv(PCC) if r["hog"] in tested}
PCG = {h for h, r in ctl.items() if int(r["fw_sites"]) > 0}
PCG2 = {e["og"] for e in ev}
PCG2s = {e["og"] for e in ev if float(e["minpp"]) >= 0.8 and e["absent_in_comparators"] == "1"}
PCGs = PCG2s

# ---- EXP
fw = tsv(FWE)
EXPog = {r["og"] for r in fw if int(r["n_fw_units_expanded"]) >= 1}
EXPog2 = {r["og"] for r in fw if int(r["n_fw_units_expanded"]) >= 2}
cafe_bg = {r["og"] for r in fw}
hk = {r["hog"]: r for r in tsv(KD + "/hog_kegg.tsv")}
ok_ = {r["og"]: r for r in tsv(KD + "/og_kegg.tsv")}
EXP = {h for h in tested if hk[h]["og"] in EXPog}
EXP2 = {h for h in tested if hk[h]["og"] in EXPog2}


def venn(a, b, c):
    return {"A_only": len(a - b - c), "B_only": len(b - a - c), "C_only": len(c - a - b), "AB": len((a & b) - c),
            "AC": len((a & c) - b), "BC": len((b & c) - a), "ABC": len(a & b & c)}


def pw(r):
    return {t for t in r["pathway_cons"].split(",") if t.startswith("map") and not glob(t)}


hpw = {h: pw(hk[h]) for h in tested}
opw = {o: pw(ok_[o]) if o in ok_ else set() for o in cafe_bg}


def enrich(sel, bgmap, label):
    bg = [x for x in bgmap if bgmap[x]]
    sel = [x for x in sel if x in bgmap and bgmap[x]]
    N, n = len(bg), len(sel)
    cnt_bg, cnt_sel, genes = defaultdict(int), defaultdict(int), defaultdict(list)
    for x in bg:
        for t in bgmap[x]:
            cnt_bg[t] += 1
    for x in sel:
        for t in bgmap[x]:
            cnt_sel[t] += 1
            genes[t].append(x)
    rows = []
    for t, a in cnt_sel.items():
        rows.append([t, names.get(t, "?"), a, n, cnt_bg[t], N, fisher_gt(a, cnt_bg[t], n, N), genes[t]])
    for r, q in zip(rows, bh([r[6] for r in rows])):
        r.append(q)
    rows.sort(key=lambda r: r[6])
    with open(OUT + "/kegg_%s.tsv" % label, "w") as w:
        w.write("pathway\tname\tk_set\tn_set\tK_bg\tN_bg\tp\tq\tmembers\n")
        for r in rows:
            w.write("%s\t%s\t%d\t%d\t%d\t%d\t%.3g\t%.3g\t%s\n" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[8],
                                                                 ",".join(sorted(r[7]))))
    return {r[0] for r in rows}, {r[0] for r in rows if r[8] < 0.05}, {r[0] for r in rows if r[6] < 0.05}


def label(h):
    r = hk[h]
    return r["preferred_name"] if r["preferred_name"] != "-" else r["description"][:40]


L, J = [], {}
pr = lambda *a: (print(*a), L.append(" ".join(str(x) for x in a)))
pr("tested HOGs", len(tested), "| all-4-unit subset", len(all4))
pr("PSG q<0.05", len(PSG), "| nominal p<0.05", len(PSGnom))
fs, cs = sum(int(r["fw_sites"]) for r in ctl.values()), sum(float(r["ctrl_sites_mean"]) for r in ctl.values())
pr("PCG main (>=3 units, pp>=0.8, comparator-absent)", len(PCG), "| sites FW", fs, "| marine 4-branch control mean %.1f"
   % cs, "| FW/control %.2f" % (fs / cs if cs else float("nan")))
pr("PCG sensitivity >=2 units: any", len(PCG2), "| pp>=0.8+absent", len(PCG2s), "| events", len(ev),
   "parallel", sum(e["type"] == "parallel" for e in ev), "convergent", sum(e["type"] == "convergent" for e in ev))
pr("CAFE families", len(cafe_bg), "| expanded >=1 unit", len(EXPog), "| >=2 units", len(EXPog2))
pr("universe HOGs whose OG expanded: >=1 unit", len(EXP), "| >=2 units", len(EXP2))
for tag, sets in [("main", (PSG, PCG, EXP)), ("nominalPSG", (PSGnom, PCG, EXP)), ("strictEXP2", (PSG, PCG, EXP2)),
                  ("PCG2units", (PSG, PCG2s, EXP)),
                  ("all4_main", (PSG & all4, PCG & all4, EXP & all4))]:
    pr("gene Venn", tag, venn(*sets))
    J["gene_" + tag] = venn(*sets)
pr("PSG∩PCG:", ";".join("%s(%s)" % (h, label(h)) for h in sorted(PSG & PCG)))
pr("PSG∩EXP:", ";".join("%s(%s)" % (h, label(h)) for h in sorted(PSG & EXP)))
pr("PCG∩EXP (n=%d) first 30:" % len(PCG & EXP), ";".join("%s(%s)" % (h, label(h)) for h in sorted(PCG & EXP)[:30]))

A = enrich(PSG if PSG else PSGnom, hpw, "PSG" if PSG else "PSGnominal")
B = enrich(PCG, hpw, "PCG")
C = enrich(EXPog, opw, "EXP")
for i, nm in enumerate(["represented", "enriched q<0.05", "enriched p<0.05"]):
    pr("pathway Venn", nm, venn(A[i], B[i], C[i]))
    J["pathway_" + nm.split()[0] + ("" if i == 0 else "_" + nm.split()[1].replace("<", "lt"))] = venn(A[i], B[i], C[i])
for i, nm in enumerate(["represented", "enriched p<0.05"]):
    j = 0 if i == 0 else 2
    s = A[j] & B[j] & C[j]
    pr("pathways shared by all three (%s):" % nm, "; ".join("%s %s" % (t, names.get(t, "?")) for t in sorted(s)))

with open(OUT + "/gene_sets.tsv", "w") as w:
    w.write("hog\tog\tname\tunits\tp\tq\tPSG\tPSG_nominal\tPCG\tPCG_fw_sites\tPCG_ctrl_mean\tPCG2_any\tPCG2_strict"
            "\tEXP\tEXP2\tbeb\n")
    for r in res:
        h = r["og"]
        c = ctl.get(h, {"fw_sites": "NA", "ctrl_sites_mean": "NA"})
        w.write("\t".join([h, hk[h]["og"], label(h), r["units"], "%.3g" % r["p"], "%.3g" % r["q"]] +
                          [str(int(h in s)) for s in (PSG, PSGnom, PCG)] + [c["fw_sites"], c["ctrl_sites_mean"]] +
                          [str(int(h in s)) for s in (PCG2, PCG2s, EXP, EXP2)] + [r["beb"]]) + "\n")
J["totals"] = {"tested": len(tested), "PSG": len(PSG), "PSG_nominal": len(PSGnom), "PCG": len(PCG), "PCG2_any": len(PCG2),
               "PCG2_strict": len(PCG2s), "EXP_hog": len(EXP), "EXP_og": len(EXPog), "EXP2_og": len(EXPog2),
               "cafe_families": len(cafe_bg), "pcg_sites_fw": fs, "pcg_sites_ctrl_mean": cs,
               "pw_PSG": len(A[0]), "pw_PCG": len(B[0]), "pw_EXP": len(C[0])}
json.dump(J, open(OUT + "/venn.json", "w"), indent=1)
open(OUT + "/summary.txt", "w").write("\n".join(L) + "\n")
