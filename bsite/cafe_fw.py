"""R1 revision (Fig. 2a, third line): freshwater-lineage expansions from the authors' existing CAFE5 run
rerun with CAFE5 (the authors' step4_cafe5/input/k2p run used a ladder tree and hit the lambda floor).

Freshwater unit stems (same four units as bsite.py, here with all 42 bivalves):
  Unionida = MRCA(Sol, Ude, Awo, Pst); Dreissenidae = MRCA(Dro, Dpo); Corbicula = Cfl tip; Limnoperna = Lfo tip.
A family counts as expanded on a unit when Gamma_change.tab > 0 on that branch and the branch is flagged significant
(Gamma_branch_probabilities.tab < 0.05) in a family with family-wide p < 0.05 (Gamma_family_results.txt).
Output cafe/fw_expansion.tsv (one row per family) and printed summary.
"""
import os, re, sys

# usage: python3 cafe_fw.py <cafe_result_dir> <model prefix, e.g. Base or Gamma> <out_dir>
K, PFX, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

# node ids from the ASR tree: tips "Xxx<n>" or "Xxx.pep<n>", internal "<n>"
tre = open(K + "/%s_asr.tre" % PFX).read()
first = tre[tre.index("("):tre.index(";", tre.index("(")) + 1]


def parse(s):
    """return parent map and tip-set per node id from a CAFE ASR newick."""
    parent, tips, stack, i, kids = {}, {}, [], 0, {}
    cur = None
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append([])
            i += 1
        elif c in ",)":
            i += 1
            if c == ")":
                ch = stack.pop()
                m = re.match(r"<(\d+)>", s[i:])
                nid = int(m.group(1))
                i += m.end()
                i = skip(s, i)
                tips[nid] = set().union(*[tips[x] for x in ch])
                for x in ch:
                    parent[x] = nid
                if stack:
                    stack[-1].append(nid)
        elif c == ";":
            break
        else:
            m = re.match(r"([A-Za-z]+)(?:\.pep)?<(\d+)>", s[i:])
            nid = int(m.group(2))
            tips[nid] = {m.group(1)}
            i += m.end()
            i = skip(s, i)
            stack[-1].append(nid)
    return parent, tips


def skip(s, i):
    m = re.match(r"\*?(_\d+)?(:[0-9.eE+-]+)?", s[i:])
    return i + m.end()


parent, tips = parse(first)
units = {"Unionida": {"Sol", "Ude", "Awo", "Pst"}, "Dreissenidae": {"Dro", "Dpo"}, "Corbicula": {"Cfl"},
         "Limnoperna": {"Lfo"}}
uid = {}
for u, s in units.items():
    cand = [n for n, t in tips.items() if s <= t]
    uid[u] = min(cand, key=lambda n: len(tips[n]))
    assert tips[uid[u]] == s, (u, tips[uid[u]])

famp = {}
for l in open(K + "/%s_family_results.txt" % PFX):
    if l.startswith("#"):
        continue
    p = l.rstrip("\n").split("\t")
    famp[p[0]] = float(p[1])


def table(path):
    f = open(path)
    h = f.readline().rstrip("\n").split("\t")
    col = {}
    for j, x in enumerate(h[1:], 1):
        m = re.search(r"<(\d+)>", x)
        col[int(m.group(1))] = j
    d = {}
    for l in f:
        p = l.rstrip("\n").split("\t")
        d[p[0]] = p
    return col, d


ccol, chg = table(K + "/%s_change.tab" % PFX)
pcol, bpr = table(K + "/%s_branch_probabilities.tab" % PFX)
rows = []
for og, pv in famp.items():
    r = {"og": og, "family_p": pv}
    nexp = 0
    for u, n in uid.items():
        c = int(chg[og][ccol[n]])
        try:
            bp = float(bpr[og][pcol[n]])
        except (KeyError, ValueError, IndexError):
            bp = float("nan")
        e = int(pv < 0.05 and c > 0 and bp < 0.05)
        r[u + "_change"], r[u + "_p"], r[u + "_exp"] = c, bp, e
        nexp += e
    r["n_fw_units_expanded"] = nexp
    rows.append(r)
keys = list(rows[0])
with open(OUT + "/fw_expansion.tsv", "w") as o:
    o.write("\t".join(keys) + "\n")
    for r in rows:
        o.write("\t".join(str(r[k]) for k in keys) + "\n")
print("families", len(rows), "family p<0.05", sum(r["family_p"] < 0.05 for r in rows))
print("unit node ids", uid)
for u in units:
    print(u, "expanded", sum(r[u + "_exp"] for r in rows))
for k in range(5):
    print(">=%d units" % k, sum(r["n_fw_units_expanded"] >= k for r in rows))
