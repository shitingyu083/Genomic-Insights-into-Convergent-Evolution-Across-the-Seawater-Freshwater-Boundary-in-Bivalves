"""R1 revision, R2 M4 (comment 12) / R1 comment 8 / Fig. 5 / P45-P46, P76-P77: species-level amino-acid composition
inputs for PGLS (aa_pgls.R). Run from _R1_work. Inputs are the server outputs of aa_comp_v2.py (server_data/r1rev/aa_v2).

Datasets (one row per species x dataset):
  ALL_ANY    all protein models of the proteome used in OrthoFinder (42 bivalves)
  OGREP_ANY  longest protein per orthogroup per species (42 bivalves)
  ALL_EXT / ALL_MEM / ALL_CYT / OGREP_EXT  eggNOG-mapper GO cellular-component classes (EXT > MEM > CYT precedence);
             Sco and Mvi excluded (about 41-45 % of their eggNOG rows cannot be assigned to one protein)
  SC8_OGC / SC36_OGC  one-to-one orthologues: members of orthogroups that are single-copy in >= 8 (all) or >= 36
             bivalves; per protein CLR (pseudocount 0.5 per residue type), centred on the orthogroup mean over the
             bivalves present, then averaged per species (within-orthologue composition deviation; 42 bivalves)
Columns: clr_<aa>; lenc = ln(mean protein length) for species-level sets, or the orthogroup-centred ln(length) for the
OGC sets; gcp = ln(GARP) - ln(FYMINK) amino-acid GC proxy (Gly/Ala/Arg/Pro vs Phe/Tyr/Met/Ile/Asn/Lys; Singer & Hickey
2000), and gcp_<aa> = the same index with <aa> left out of its group (identical to gcp for the other ten residues).
"""
import csv, gzip, math
from collections import defaultdict

AA = "ACDEFGHIKLMNPQRSTVWY"
GARP, FYMINK = "GARP", "FYMINK"
SRC = "server_data/r1rev/aa_v2"
OUT = "analysis/aa"
SPLIT = {"Sco", "Mvi"}
SETS = [("ALL", "ANY"), ("OGREP", "ANY"), ("ALL", "EXT"), ("ALL", "MEM"), ("ALL", "CYT"), ("OGREP", "EXT")]


def gc_proxies(v):
    """v: dict aa -> positive value (frequency or count + pseudocount). Returns gcp and leave-one-out variants."""
    g, f = sum(v[a] for a in GARP), sum(v[a] for a in FYMINK)
    out = {"gcp": math.log(g) - math.log(f)}
    for a in AA:
        if a in GARP:
            out["gcp_" + a] = math.log(g - v[a]) - math.log(f)
        elif a in FYMINK:
            out["gcp_" + a] = math.log(g) - math.log(f - v[a])
        else:
            out["gcp_" + a] = out["gcp"]
    return out


def clr(v):
    lg = {a: math.log(v[a]) for a in AA}
    m = sum(lg.values()) / 20
    return {a: lg[a] - m for a in AA}


sal = {r["id"]: r["salinity3"] for r in csv.DictReader(open("analysis/pgls_span/traits_S1_parsed.tsv"), delimiter="\t")}
rows = []
for r in csv.DictReader(open(SRC + "/aa_species_sets.tsv"), delimiter="\t"):
    if (r["set"], r["class"]) not in SETS or r["species"] == "Lgi":
        continue
    if r["class"] != "ANY" and r["species"] in SPLIT:
        continue
    c = {a: int(r[a]) for a in AA}
    tot = sum(c.values())
    fr = {a: c[a] / tot for a in AA}
    o = {"species": r["species"], "salinity3": sal[r["species"]], "dataset": r["set"] + "_" + r["class"],
         "n_prot": int(r["n_prot"]), "n_res20": tot, "lenc": math.log(int(r["n_res"]) / int(r["n_prot"]))}
    o.update({"f_" + a: fr[a] for a in AA})
    o.update({"clr_" + a: x for a, x in clr(fr).items()})
    o.update(gc_proxies(fr))
    rows.append(o)

# one-to-one orthologues: orthogroup-centred composition (numpy; Lgi excluded before centring)
import numpy as np

IX = {a: i for i, a in enumerate(AA)}
gi, fi = [IX[a] for a in GARP], [IX[a] for a in FYMINK]
ogs, sps, lens, cnt = [], [], [], []
with gzip.open(SRC + "/aa_og_sc.tsv.gz", "rt") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["species"] == "Lgi":
            continue
        ogs.append(r["og"]); sps.append(r["species"]); lens.append(int(r["len"])); cnt.append([int(r[a]) for a in AA])
assert len(set(zip(ogs, sps))) == len(ogs), "duplicate orthogroup x species row"
V = np.asarray(cnt, dtype=float) + 0.5
LV = np.log(V)
g, fy = V[:, gi].sum(1), V[:, fi].sum(1)
loo = np.empty_like(V)
for a, i in IX.items():
    loo[:, i] = (np.log(g - V[:, i]) - np.log(fy)) if a in GARP else \
        (np.log(g) - np.log(fy - V[:, i])) if a in FYMINK else (np.log(g) - np.log(fy))
X = np.column_stack([LV - LV.mean(1, keepdims=True), np.log(lens), np.log(g) - np.log(fy), loo])
names = ["clr_" + a for a in AA] + ["lenc", "gcp"] + ["gcp_" + a for a in AA]
ou, oi = np.unique(ogs, return_inverse=True)
su, si = np.unique(sps, return_inverse=True)
nsc = np.bincount(oi)
assert nsc.min() >= 8, nsc.min()
S = np.zeros((len(ou), X.shape[1]))
np.add.at(S, oi, X)
Xc = X - (S / nsc[:, None])[oi]
for lab, keep in [("SC8_OGC", nsc[oi] >= 8), ("SC36_OGC", nsc[oi] >= 36)]:
    T = np.zeros((len(su), X.shape[1]))
    np.add.at(T, si[keep], Xc[keep])
    k = np.bincount(si[keep], minlength=len(su))
    for j, sp in enumerate(su):
        o = {"species": sp, "salinity3": sal[sp], "dataset": lab, "n_prot": int(k[j]),
             "n_res20": int((V[keep & (si == j)] - 0.5).sum())}
        o.update({nm: T[j, c] / k[j] for c, nm in enumerate(names)})
        rows.append(o)
    print(lab, "orthogroups", int((nsc >= (8 if lab == "SC8_OGC" else 36)).sum()), "proteins", int(keep.sum()),
          "per-species min/max", int(k.min()), int(k.max()))

cols = ["species", "salinity3", "dataset", "n_prot", "n_res20", "lenc", "gcp"] + ["gcp_" + a for a in AA] + \
       ["clr_" + a for a in AA] + ["f_" + a for a in AA]
with open(OUT + "/aa_species_clr.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for o in rows:
        w.writerow([("%.8g" % o[c]) if isinstance(o.get(c), float) else o.get(c, "NA") for c in cols])
n = defaultdict(int)
for o in rows:
    n[o["dataset"]] += 1
print({k: v for k, v in n.items()})
