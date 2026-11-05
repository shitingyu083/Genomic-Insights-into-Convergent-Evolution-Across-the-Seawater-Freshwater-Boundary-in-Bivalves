"""Diagnostic: how many tested HOGs carry parallel/convergent sites at increasing stringency (site = og+site;
units = number of freshwater unit stems sharing the derived residue)."""
import sys
from collections import defaultdict

RUN = sys.argv[1]
ok = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(RUN + "/bsite_results.tsv") if "\tok\t" in l}
site = defaultdict(lambda: {"units": set(), "minpp": 1.0, "absent": 1})
for l in open(RUN + "/pc_events.tsv"):
    p = l.rstrip("\n").split("\t")
    if p[0] not in ok:
        continue
    s = site[(p[0], p[1])]
    s["units"].update(p[2].split("-"))
    s["minpp"] = min(s["minpp"], float(p[5]))
    s["absent"] = min(s["absent"], int(p[6]))
print("tested", len(ok), "sites", len(site))
for k in (2, 3, 4):
    for pp in (0, 0.8, 0.95):
        for ab in (0, 1):
            g = {o for (o, _), s in site.items() if len(s["units"]) >= k and s["minpp"] >= pp and s["absent"] >= ab}
            n = sum(1 for s in site.values() if len(s["units"]) >= k and s["minpp"] >= pp and s["absent"] >= ab)
            print("units>=%d minpp>=%.2f absent=%d : sites %d genes %d" % (k, pp, ab, n, len(g)))
