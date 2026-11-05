# Repository copy (2026-11-05): section 4 of the working script, read-only checks against files of the
# original submission, is omitted; sections 1-3 are unchanged.
# Candidate-pool functional prioritisation for R1-4 / R2-16 on the v4 convergence reanalysis.
# Test set   : orthogroups (OGs) with >= 1 recurrent-substitution candidate site (A5, 127 sites).
# Background : the 99 OGs analysed in the v4 reanalysis (A1).
# Annotation : eggNOG-mapper 2.1.13 rows extracted read-only on the server (ann_targets.tsv); rows flagged
#              AMBIGUOUS (protein key assigned to more than one orthogroup) are dropped.
# OG terms   : 'majority' = term carried by >= 50 % of the members whose field is non-empty;
#              'any' = term carried by any member (sensitivity).
# Test       : exact one-sided hypergeometric over-representation (checked against scipy); Benjamini-Hochberg
#              within each annotation type x rule over every term with >= MIN_K background OGs.
# Background per type = OGs with >= 1 assigned term of that type (clusterProfiler convention).
import os, csv, json, math, hashlib, collections
from fractions import Fraction
from scipy.stats import hypergeom

W = (os.environ.get("IZ_ROOT", ".") + "/_R1_work")
A = os.path.join(W, "reanalysis/A_convergence_v4")
ANN = os.path.join(W, "server_data/r1rev/ann/ann/ann_targets.tsv")
OUT = os.path.join(W, "analysis/enrich")
MIN_K = 2
KINDS = ["KEGG_Pathway", "GOs", "PFAMs", "KEGG_ko"]


def md5(p):
    with open(p, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def read_tsv(p):
    with open(p, encoding="utf8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_tsv(name, rows, cols):
    with open(os.path.join(OUT, name), "w", encoding="utf8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])


U = sorted({r["og"] for r in read_tsv(os.path.join(A, "A1_unit_presence_by_OG.tsv"))})
A5 = read_tsv(os.path.join(A, "A5_unit_recurrent_candidates.tsv"))
C = sorted({r["og"] for r in A5})
assert len(U) == 99 and len(C) == 37 and set(C) <= set(U)
sites = collections.defaultdict(collections.Counter)
for r in A5:
    sites[r["og"]][r["type"]] += 1

# Reference names: KEGG REST list/pathway (downloaded 2026-10-25) and GO go-basic.obo (release 2026-07-26).
# The 32 MB obo is a temporary download: names for every GO ID used are kept in go_terms_used.tsv, which is
# read instead when the obo has been removed.
KNAME = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(OUT, "kegg_list_pathway_map.txt"), encoding="utf8") if "\t" in l)
OBO, GOUSED = os.path.join(OUT, "go-basic.obo"), os.path.join(OUT, "go_terms_used.tsv")
GONAME, GONS, ALT, cur = {}, {}, {}, None
if os.path.exists(OBO):
    with open(OBO, encoding="utf8") as fh:
        for line in fh:
            if line.startswith("["):
                cur = None
            elif line.startswith("data-version: "):
                OBO_VERSION = line[14:].strip()
            elif line.startswith("id: GO:"):
                cur = line[4:].strip()
            elif cur and line.startswith("name: "):
                GONAME[cur] = line[6:].strip()
            elif cur and line.startswith("namespace: "):
                GONS[cur] = line[11:].strip()
            elif cur and line.startswith("alt_id: "):
                ALT[line[8:].strip()] = cur
    OBO_MD5 = md5(OBO)
else:
    rows = read_tsv(GOUSED)
    OBO_VERSION, OBO_MD5 = rows[0]["obo_version"], rows[0]["obo_md5"]
    for r in rows:
        ALT[r["go_id"]] = r["primary_id"]
        if r["name"]:  # IDs absent from the obo stay unnamed, as in the obo branch
            GONAME[r["primary_id"]], GONS[r["primary_id"]] = r["name"], r["namespace"]

with open(ANN, encoding="utf8", newline="") as fh:
    hdr = fh.readline().rstrip("\r\n").split("\t")
    raw = [l.rstrip("\r\n").split("\t") for l in fh]
ix = {h: i for i, h in enumerate(hdr)}
assert all(len(f) in (23, 24) for f in raw) and all(f[23] == "AMBIGUOUS" for f in raw if len(f) == 24)
Uset, Cset = set(U), set(C)
mem = collections.defaultdict(list)
for f in raw:
    if len(f) == 23 and f[1] in Uset:
        mem[f[1]].append(f)
n_amb_U = sum(1 for f in raw if len(f) == 24 and f[1] in Uset)


def terms(f, kind):
    v = f[ix[kind]]
    if v in ("", "-"):
        return set()
    items = set(v.split(","))
    if kind == "KEGG_Pathway":
        items = {t for t in items if t.startswith("map")}
    elif kind == "GOs":
        items = {t for t in items if t.startswith("GO:")}
    return items


assign = {}
for kind in KINDS:
    for rule in ("majority", "any"):
        d = {}
        for og, fs in mem.items():
            ts = [s for s in (terms(f, kind) for f in fs) if s]
            if not ts:
                continue
            cnt = collections.Counter(t for s in ts for t in s)
            keep = {t for t, c in cnt.items() if rule == "any" or 2 * c >= len(ts)}
            if keep:
                d[og] = keep
        assign[(kind, rule)] = d


def top_value(og, col):
    vals = collections.Counter(f[ix[col]] for f in mem.get(og, []) if f[ix[col]] not in ("", "-"))
    if not vals:
        return "-", 0, 0
    best = max(vals.values())
    return "|".join(sorted(v for v, c in vals.items() if c == best)), best, sum(vals.values())


gene = {og: top_value(og, "Preferred_name")[0] for og in U}


def label(kind, t):
    if kind == "KEGG_Pathway":
        return KNAME.get(t, "")
    if kind == "GOs":
        g = ALT.get(t, t)
        return GONAME.get(g, "") + (" [" + GONS.get(g, "") + "]" if g in GONS else "")
    return ""


def p_exact(k, M, K, n):
    num = sum(math.comb(K, x) * math.comb(M - K, n - x) for x in range(k, min(K, n) + 1))
    return float(Fraction(num, math.comb(M, n)))


def bh(ps):
    m, run, q = len(ps), 1.0, [1.0] * len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        run = min(run, ps[i] * m / rank)
        q[i] = run
    return q


def fmt(x, d=4):
    if isinstance(x, float):
        return "NA" if math.isnan(x) else f"{x:.{d}g}"
    return "NA" if x is None else str(x)


# ---- 1. Over-representation of candidate OGs within the analysed pool ----
def enrich(kind, rule, min_k):
    d = assign[(kind, rule)]
    bg = set(d)
    cand = bg & Cset
    M, n = len(bg), len(cand)
    cnt_bg = collections.Counter(t for og in bg for t in d[og])
    cnt_c = collections.Counter(t for og in cand for t in d[og])
    rows = []
    for t in sorted(t for t, K in cnt_bg.items() if K >= min_k):
        K, k = cnt_bg[t], cnt_c.get(t, 0)
        p = p_exact(k, M, K, n)
        assert math.isclose(p, float(hypergeom.sf(k - 1, M, K, n)), rel_tol=1e-9, abs_tol=1e-12)
        ogs = sorted(og for og in cand if t in d[og])
        e = n * K / M
        rows.append(dict(type=kind, rule=rule, term=t, name=label(kind, t), k_cand=k, n_cand=n, K_bg=K, M_bg=M,
                         expected=e, fold=k / e if e else float("nan"), rich_factor=k / K, p=p,
                         global_map=int(kind == "KEGG_Pathway" and t.startswith("map01")),
                         cand_ogs=",".join(ogs) or "-", cand_genes=",".join(gene[og] for og in ogs) or "-"))
    for r, q in zip(rows, bh([r["p"] for r in rows])):
        r["q_bh"] = q
    rows.sort(key=lambda r: (r["p"], r["term"]))
    fam = dict(M_bg=M, n_cand=n, terms_seen=len(cnt_bg), terms_tested=len(rows),
               min_p=min((r["p"] for r in rows), default=None), min_q=min((r["q_bh"] for r in rows), default=None),
               n_p_lt_0_05=sum(r["p"] < 0.05 for r in rows), n_q_lt_0_05=sum(r["q_bh"] < 0.05 for r in rows),
               n_q_lt_0_10=sum(r["q_bh"] < 0.10 for r in rows))
    return rows, fam


results, families, families_mink5 = [], {}, {}
for kind, rule in assign:
    if kind == "KEGG_ko":  # a KO is close to a single gene identity; kept for the per-OG table only
        continue
    rows, families[kind + "|" + rule] = enrich(kind, rule, MIN_K)
    results += rows
    families_mink5[kind + "|" + rule] = enrich(kind, rule, 5)[1]  # sensitivity: minimum term size 5
RES_COLS = ["type", "rule", "term", "name", "k_cand", "n_cand", "K_bg", "M_bg", "expected", "fold", "rich_factor",
            "p", "q_bh", "global_map", "cand_ogs", "cand_genes"]
write_tsv("enrich_results.tsv", [{c: fmt(r[c]) for c in RES_COLS} for r in results], RES_COLS)

# ---- 2. Per-OG annotation summary for the 99 analysed OGs ----
units = collections.defaultdict(set)
for r in A5:
    units[r["og"]].update(r["units"].split(";"))
og_rows = []
for og in U:
    fs, s = mem.get(og, []), sites.get(og, collections.Counter())
    assert len({f[0] for f in fs}) == len(fs)  # at most one annotated member per species and OG

    def maj(kind):
        return sorted(assign[(kind, "majority")].get(og, set()))

    nm, best, tot = top_value(og, "Preferred_name")
    og_rows.append(dict(
        og=og, candidate=int(og in Cset), sites_parallel=s["parallel"], sites_convergent=s["convergent"],
        sites_mixed=s["mixed"], units_with_candidates=";".join(sorted(units[og])) or "-", n_annot_members=len(fs),
        preferred_name=nm, name_support=f"{best}/{tot}", description=top_value(og, "Description")[0],
        kegg_ko_majority=",".join(maj("KEGG_ko")) or "-",
        kegg_map_majority="; ".join(f"{t} {KNAME.get(t, '')}" for t in maj("KEGG_Pathway")) or "-",
        pfam_majority=",".join(maj("PFAMs")) or "-", n_go_majority=len(maj("GOs"))))
assert {t for c in sites.values() for t in c} <= {"parallel", "convergent", "mixed"}
assert sum(r["sites_parallel"] + r["sites_convergent"] + r["sites_mixed"] for r in og_rows) == len(A5) == 127
write_tsv("enrich_og_terms.tsv", og_rows, list(og_rows[0]))

# ---- 3. Pathways named in the submitted Fig. 2c, Results P67 and Discussion P85 ----
# Flags: Fig. 2c parallel panel, Fig. 2c 'Postive' panel, P67 recurrent-set list, P67 branch-site list, P85 list.
# Fig. 2c 'Autophagy' is not species-qualified, so both KEGG autophagy maps are listed; P85 'transporter
# pathways' is matched to ABC transporters (approximation).
FLAGS = ["fig2c_parallel", "fig2c_positive", "p67_recurrent", "p67_branch_site", "p85"]
NAMED = [("map04020", 1, 0, 1, 0, 1), ("map04024", 1, 0, 1, 0, 1), ("map04152", 1, 0, 1, 0, 1),
         ("map04151", 1, 1, 1, 0, 1), ("map04141", 1, 1, 1, 1, 1), ("map04140", 1, 1, 1, 1, 1),
         ("map04138", 1, 1, 1, 1, 1), ("map04022", 1, 0, 0, 0, 0), ("map04015", 1, 0, 0, 0, 0),
         ("map04071", 1, 0, 0, 0, 0), ("map00330", 1, 0, 0, 0, 0), ("map00250", 1, 1, 0, 0, 0),
         ("map00010", 1, 0, 0, 0, 0), ("map00190", 1, 1, 0, 1, 1), ("map04210", 0, 1, 0, 1, 1),
         ("map04066", 0, 1, 0, 1, 0), ("map04010", 0, 1, 0, 1, 0), ("map00590", 0, 1, 0, 0, 0),
         ("map00071", 0, 1, 0, 0, 0), ("map02010", 0, 1, 0, 1, 1), ("map00770", 0, 1, 0, 1, 0),
         ("map00480", 0, 1, 0, 0, 0)]
hit = {(r["rule"], r["term"]): r for r in results if r["type"] == "KEGG_Pathway"}
named_rows = []
for mp, *flags in NAMED:
    assert mp in KNAME, mp
    row = dict(map_id=mp, name=KNAME[mp], **dict(zip(FLAGS, flags)))
    for rule in ("majority", "any"):
        d, fam = assign[("KEGG_Pathway", rule)], families["KEGG_Pathway|" + rule]
        ogs = sorted(og for og in d if mp in d[og])
        h = hit.get((rule, mp))
        row[rule + "_k_of_n"] = f"{sum(og in Cset for og in ogs)}/{fam['n_cand']}"
        row[rule + "_K_of_M"] = f"{len(ogs)}/{fam['M_bg']}"
        row[rule + "_p"] = fmt(h["p"]) if h else f"not tested (K<{MIN_K})"
        row[rule + "_q"] = fmt(h["q_bh"]) if h else "-"
        row[rule + "_ogs"] = ",".join(f"{og}:{gene[og]}{'*' if og in Cset else ''}" for og in ogs) or "-"
    named_rows.append(row)
write_tsv("enrich_named_pathways.tsv", named_rows, list(named_rows[0]))
