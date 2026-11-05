"""R1 revision: KEGG annotation of (i) branch-site universe HOGs (N1) and (ii) all orthogroups (CAFE families).
Queries of the eggNOG-mapper run are matched to orthogroup members with the same key() normalisation as
ann_extract.py; queries whose key maps to >1 (species, OG) are skipped. A pathway (map IDs only) / KO is assigned to
a HOG/OG when it is carried by >= 50% of its annotated member proteins (consensus), and the union is also kept.
usage: python3 kegg_map.py <universe.tsv> <out_dir>
"""
import csv, os, re, sys
from collections import defaultdict, Counter

HOME = os.path.expanduser("~")
R = os.environ.get("R0_DATA", HOME + "/###qutongjinhua_0710cxcf") + ""
RES = R + "/ortho_clean/orthofinder_output/Results_Dec31"
EM = R + "/orthofinder1229/workflow/convergence_analysis_workflow/step4_cafe5/input/eggnog_output.xlsx.emapper.annotations"
UNI, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)


def key(s):
    s = re.sub(r"[^\x20-\x7e]", "", s)
    return re.sub(r"[\s:,()]", "_", s).upper()


def unq(g):
    g = g.strip()
    if len(g) > 1 and g[0] == g[-1] == '"':
        g = g[1:-1].replace('""', '"')
    return g


uni = {l.split("\t")[0] for l in open(UNI) if l.startswith("HOG")}

# protein key -> (species, OG) from Orthogroups.tsv
g2 = defaultdict(set)
with open(RES + "/Orthogroups/Orthogroups.tsv", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
    sps = [re.sub(r"\.fa$", "", h).replace(".pep", "") for h in next(rd)[1:]]
    for row in rd:
        for sname, cell in zip(sps, row[1:]):
            for g in cell.split(", "):
                g = unq(g)
                if g:
                    g2[key(g)].add((sname, row[0]))

# protein key -> universe HOG from N1.tsv
k2hog, hog2og = {}, {}
with open(RES + "/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    for l in f:
        x = l.rstrip("\n").split("\t")
        h = x[0].replace("N1.", "")
        if h not in uni:
            continue
        hog2og[h] = x[1]
        for cell in x[3:]:
            for g in cell.split(", "):
                g = unq(g)
                if g:
                    k2hog[key(g)] = h

ogn, ogp, ogk, ogname = Counter(), defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
hn, hp, hk, hname, hdesc = Counter(), defaultdict(Counter), defaultdict(Counter), defaultdict(Counter), defaultdict(Counter)
n = amb = none = 0
with open(EM, encoding="utf-8", errors="replace") as f:
    for l in f:
        if l.startswith("#"):
            continue
        p = l.rstrip("\n").split("\t")
        if len(p) < 13:
            continue
        n += 1
        k = key(p[0])
        hits = g2.get(k)
        if not hits:
            none += 1
            continue
        if len(hits) > 1:
            amb += 1
            continue
        og = next(iter(hits))[1]
        pw = {t for t in p[12].split(",") if t.startswith("map")}
        ko = {t for t in p[11].split(",") if t.startswith("ko:")}
        nm = p[8] if p[8] not in ("", "-") else ""
        ogn[og] += 1
        for t in pw: ogp[og][t] += 1
        for t in ko: ogk[og][t] += 1
        if nm: ogname[og][nm] += 1
        h = k2hog.get(k)
        if h:
            hn[h] += 1
            for t in pw: hp[h][t] += 1
            for t in ko: hk[h][t] += 1
            if nm: hname[h][nm] += 1
            if p[7] not in ("", "-"): hdesc[h][p[7]] += 1


def cons(c, m):
    return ",".join(sorted(t for t, v in c.items() if v >= 0.5 * m)) or "-"


def union(c):
    return ",".join(sorted(c)) or "-"


def top(c):
    return c.most_common(1)[0][0] if c else "-"


with open(OUT + "/hog_kegg.tsv", "w") as w:
    w.write("hog\tog\tn_annot\tpreferred_name\tdescription\tko_cons\tpathway_cons\tpathway_union\n")
    for h in sorted(uni):
        m = hn[h]
        w.write("\t".join([h, hog2og.get(h, "-"), str(m), top(hname[h]), top(hdesc[h]), cons(hk[h], m),
                           cons(hp[h], m), union(hp[h])]) + "\n")
with open(OUT + "/og_kegg.tsv", "w") as w:
    w.write("og\tn_annot\tpreferred_name\tko_cons\tpathway_cons\tpathway_union\n")
    for og in sorted(ogn):
        m = ogn[og]
        w.write("\t".join([og, str(m), top(ogname[og]), cons(ogk[og], m), cons(ogp[og], m), union(ogp[og])]) + "\n")
print("emapper rows", n, "ambiguous", amb, "no OG", none, "OGs annotated", len(ogn),
      "universe HOGs", len(uni), "annotated", sum(hn[h] > 0 for h in uni),
      "with consensus pathway", sum(cons(hp[h], hn[h]) != "-" for h in uni))
