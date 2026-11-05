"""Control-branch calibration of the parallel/convergent (PCG) criterion.

For every tested HOG, re-reads the aaml ancestral reconstruction (run/og/<HOG>/asr/rst) and scores each site
with the same rule on two branch sets:
  FW   : the four freshwater unit stems (Unionida, Dreissenidae, Corbicula, Limnoperna);
  CTRL : every 4-subset of the marine comparator terminal branches present in the HOG.
A site counts when >= K branches of the set carry a change (ancestor != derived, both pp >= PP) to the same
derived residue, and that residue is absent from every tip outside the set.
Output: per-HOG FW count and mean control count (averaged over 4-subsets), plus genome-wide totals.
Usage: python3 pc_control.py <run> <out.tsv> [K=3] [PP=0.8]
"""
import sys
import itertools
from collections import defaultdict
import bsite as B

RUN, OUT = sys.argv[1], sys.argv[2]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 3
PP = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
ok = [l.split("\t")[0] for l in open(RUN + "/bsite_results.tsv") if "\tok\t" in l]


def order_of(wd):
    ls = open(wd + "/aa.phy").read().split("\n")[1:]
    return [l.split()[0] for l in ls if l.strip()]


def score(branches, sites, st, parent, outside_tips):
    """branches: dict name -> node id; returns list of qualifying site indices."""
    hits = []
    for s, (data, _) in enumerate(sites, 1):
        der = defaultdict(int)
        for n in branches.values():
            a, pa = st(s, parent[n])
            b, pb = st(s, n)
            if a != b and min(pa, pb) >= PP:
                der[b] += 1
        for b, c in der.items():
            if c >= K and b not in [data[i - 1] for i in outside_tips]:
                hits.append(s)
                break
    return hits


fo = open(OUT, "w")
fo.write("hog\tn_ctrl_branches\tn_ctrl_subsets\tfw_sites\tctrl_sites_mean\tctrl_any_subset\n")
tot_fw = tot_ctrl = genes_fw = genes_ctrl = n = 0
for og in ok:
    wd = RUN + "/og/" + og
    try:
        order = order_of(wd)
        parent, root, nt, sites = B.rst_parse(wd + "/asr", order)
    except Exception as e:  # noqa
        print("skip", og, e, file=sys.stderr)
        continue
    idx = {sp: i + 1 for i, sp in enumerate(order)}

    def st(s, node, sites=sites, nt=nt):
        data, states = sites[s - 1]
        if node <= nt:
            return data[node - 1], 1.0
        k = node - nt - 1
        return states[k][0], float(states[k][1])

    def mrca(members):
        anc = []
        for sp in members:
            a, x = [], idx[sp]
            while x in parent:
                a.append(x)
                x = parent[x]
            a.append(x)
            anc.append(a)
        return [m for m in anc[0] if all(m in a for a in anc[1:])][0]

    fw = {u: mrca([x for x in m if x in idx]) for u, m in B.UNITS.items() if any(x in idx for x in m)}
    fw_tips = {idx[sp] for sp in order if sp in B.FW}
    mar = [sp for sp in order if sp not in B.FW]
    fw_hits = score(fw, sites, st, parent, [i for i in idx.values() if i not in fw_tips])
    subs = list(itertools.combinations(mar, 4))
    ch = []
    anyhit = set()
    for sub in subs:
        br = {sp: idx[sp] for sp in sub}
        h = score(br, sites, st, parent, [i for sp, i in idx.items() if sp not in sub])
        ch.append(len(h))
        anyhit.update(h)
    m = sum(ch) / len(ch) if ch else 0.0
    fo.write("%s\t%d\t%d\t%d\t%.4f\t%d\n" % (og, len(mar), len(subs), len(fw_hits), m, len(anyhit)))
    n += 1
    tot_fw += len(fw_hits)
    tot_ctrl += m
    genes_fw += len(fw_hits) > 0
    genes_ctrl += m > 0
fo.close()
print("K=%d PP=%.2f tested=%d | FW sites %d genes %d | control mean sites %.1f genes(mean>0) %d | ratio %.2f"
      % (K, PP, n, tot_fw, genes_fw, tot_ctrl, genes_ctrl, tot_fw / tot_ctrl if tot_ctrl else float("nan")))
