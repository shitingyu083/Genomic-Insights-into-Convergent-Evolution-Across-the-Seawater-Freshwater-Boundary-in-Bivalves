import os
"""R1 revision, R2 M3 (Fig. 3b): are the gene-tree tips distinct loci?

Collapses the protein models (tips) of the three Fig. 3b orthogroups to annotation genes and to sequence-level loci.
  L1  protein models (= OrthoFinder copy number = tips in the submitted gene trees)
  L2  L1 joined when two models share an annotation gene ID (MAKER -R*, -T*, -mRNA*, .M*, .T*, Ensembl GENE:,
      GENE=, EVM TU), form one RefSeq "isoform X#" series, overlap on the same scaffold and strand, or are identical
  L3  L2 joined when pid_aligned >= PID, cov_shorter >= 0.90 and shared 15-mers >= 0.5 (single linkage; within species)
  full-length L3 locus = longest member >= 50% of the reference length (median over bivalves of the per-species
      longest member)
Inputs (all real server output): aa_v2/fig3b/fig3b_members.tsv, aa_v2/fig3b/fig3b_pairs.tsv (aa_comp_v2.py,
  md5 6c1958f4), fig3b_em/fig3b_members_emapper.tsv (fig3b_em.py), analysis/pgls_span/traits_S1_parsed.tsv.
"""
import collections, csv, re, statistics, sys

W = (os.environ.get("IZ_ROOT", ".") + "/_R1_work")
SD = W + "/server_data/r1rev"
OUT = W + "/analysis/fig3b"
FAM = {"OG0003585": "Reelin-like (RELN)", "OG0004390": "Delphilin-like (GRID2IP)", "OG0005024": "PMS2"}
# submitted Fig. 3b panel titles; panel-to-orthogroup assignment from the printed tip labels (fig3b_panel_tips.py)
SUB = {"OG0005024": "left: ATG8/LC3/GABARAP", "OG0003585": "middle: Beclin 1", "OG0004390": "right: Aquaporin/AQP"}
REFSEQ = {"Can", "Cgi", "Cvi", "Dpo", "Mar", "Mca", "Mme", "Oed", "Rph", "Yba"}
GENE_RE = {
    "Cni": r"^(CNI_\d+)-R[A-Z]+", "Awo": r"^(AWO\w+?G\d+)\.M\d+", "Lfo": r"^(.+?)\.T\d+$", "Pfu": r"^(.+?)\.T\d+$",
    "Mco": r"EVM(LG\d+\.G\d+)", "Mve": r"EVM\.TU\.([A-Z0-9]+\.\d+)", "Pst": r"^(PST_\d+)-T\d+",
    "Sgl": r"^(SGL\d+)-MRNA\d+", "Sso": r"GENE:(ENSAIOG\d+)", "Ude": r"GENE=(UDEG\d+)", "Sbr": r"GENE=(EVM\d+)",
    "Cbi": r"^(.+?)LOCUS=", "Cho": r"^(.+?)LOCUS=", "Sol": r"^(.+?)LOCUS="}
COORD_RE = {"Sso": r"PEPXBSPISOLI1\.1:([^:]+):(\d+):(\d+):(-?1)GENE", "Sol": r"LOCUS=(.+?):(\d+):(\d+):([+-])",
            "Cho": r"LOCUS=(.+?):(\d+):(\d+):([+-])", "Cbi": r"LOCUS=(.+?):(\d+):(\d+):([+-])",
            "Mve": r"(CM\d+\.\d+):(\d+)-(\d+)\(([+-])\)", "Mco": r"(LG\d+):(\d+)-(\d+)\(([+-])\)"}


def rd(p):
    with open(p, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


class UF:
    def __init__(s, xs):
        s.p = {x: x for x in xs}

    def f(s, x):
        while s.p[x] != x:
            s.p[x] = s.p[s.p[x]]
            x = s.p[x]
        return x

    def u(s, a, b):
        s.p[s.f(a)] = s.f(b)

    def groups(s):
        g = collections.defaultdict(list)
        for x in s.p:
            g[s.f(x)].append(x)
        return list(g.values())


def gene_id(sp, name):
    m = re.search(GENE_RE[sp], name) if sp in GENE_RE else None
    return m.group(1) if m else name


def refseq(name):
    m = re.match(r"^(XP_\d+)\.\d+(.*?)(?:ISOFORMX(\d+))?\[", name)
    return (m.group(2), m.group(3)) if m else (None, None)


def coord(sp, name):
    m = re.search(COORD_RE[sp], name) if sp in COORD_RE else None
    return (m.group(1), min(int(m.group(2)), int(m.group(3))), max(int(m.group(2)), int(m.group(3))), m.group(4)) if m else None


def model_num(sp, name):
    """(prefix, number) of the model / accession ID, used only to flag consecutive IDs."""
    s = re.match(r"^(XP_\d+)", name).group(1) if sp in REFSEQ and name.startswith("XP_") else gene_id(sp, name)
    s = re.sub(r"\.\d+$", "", s) if sp not in ("Lfo", "Pfu", "Mve", "Mco", "Air", "Apu", "Sco", "Mvi", "Csi", "Tna") else s
    m = re.match(r"^(.*?)(\d+)$", s)
    return (m.group(1), int(m.group(2))) if m else None


def main(pid_thr=0.95, write=True):
    mem = rd(SD + "/aa_v2/fig3b/fig3b_members.tsv")
    pairs = rd(SD + "/aa_v2/fig3b/fig3b_pairs.tsv")
    em = {r["iid"]: r for r in rd(SD + "/fig3b_em/fig3b_members_emapper.tsv")}
    sal = {r["id"]: r["salinity3"] for r in rd(W + "/analysis/pgls_span/traits_S1_parsed.tsv")}
    sal["Lgi"] = "OUTGROUP"
    for r in mem:
        r["len"] = int(r["len"])
    by = collections.defaultdict(list)
    for r in mem:
        by[(r["og"], r["species"])].append(r)
    ref = {og: statistics.median(max(x["len"] for x in v) for (o, sp), v in by.items() if o == og and sp != "Lgi")
           for og in FAM}
    pw = {}
    for p in pairs:
        if p["comparison"].startswith("within_"):
            pw[(p["og"], p["iid_a"], p["iid_b"])] = p
    rows, mrows = [], []
    for (og, sp), v in sorted(by.items()):
        ids = [x["iid"] for x in v]
        u2 = UF(ids)
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                a, b = v[i], v[j]
                same = gene_id(sp, a["name"]) == gene_id(sp, b["name"]) or a["md5"] == b["md5"]
                ca, cb = coord(sp, a["name"]), coord(sp, b["name"])
                if ca and cb and ca[0] == cb[0] and ca[3] == cb[3] and ca[1] <= cb[2] and cb[1] <= ca[2]:
                    same = True
                if same:
                    u2.u(a["iid"], b["iid"])
        if sp in REFSEQ:  # one RefSeq "isoform X#" series = one gene, unless an X# repeats among distinct sequences
            ser = collections.defaultdict(list)
            for x in v:
                st, iso = refseq(x["name"])
                if iso:
                    ser[st].append(x)
            for st, xs in ser.items():
                iso_md5 = collections.defaultdict(set)
                for x in xs:
                    iso_md5[refseq(x["name"])[1]].add(x["md5"])
                if all(len(s) == 1 for s in iso_md5.values()):
                    for x in xs[1:]:
                        u2.u(xs[0]["iid"], x["iid"])
        u3 = UF(ids)
        for g in u2.groups():
            for x in g[1:]:
                u3.u(g[0], x)
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                a, b = v[i]["iid"], v[j]["iid"]
                p = pw.get((og, a, b)) or pw.get((og, b, a))
                if p is None:
                    sys.exit("missing within-species pair %s %s %s" % (og, a, b))
                if (float(p["pid_aligned"]) >= pid_thr and float(p["cov_shorter"]) >= 0.90
                        and float(p["shared15mer_shorter"]) >= 0.5):
                    u3.u(a, b)
        g2, g3 = u2.groups(), u3.groups()
        ln = {x["iid"]: x["len"] for x in v}
        full = [g for g in g3 if max(ln[i] for i in g) >= 0.5 * ref[og]]
        # flags (not merged): adjacent non-overlapping models on one scaffold/strand (<100 kb) or consecutive IDs
        cid = {x["iid"]: next(k for k, g in enumerate(g3) if x["iid"] in g) for x in v}
        adj = set()
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                a, b = v[i], v[j]
                if cid[a["iid"]] == cid[b["iid"]]:
                    continue
                ca, cb = coord(sp, a["name"]), coord(sp, b["name"])
                na, nb = model_num(sp, a["name"]), model_num(sp, b["name"])
                if ca and cb:
                    if ca[0] == cb[0] and ca[3] == cb[3] and max(ca[1], cb[1]) - min(ca[2], cb[2]) < 1e5:
                        adj.add(tuple(sorted((cid[a["iid"]], cid[b["iid"]]))))
                elif na and nb and na[0] == nb[0]:
                    step = 10 if all(model_num(sp, x["name"]) and model_num(sp, x["name"])[1] % 10 == 0 for x in v) else 1
                    if abs(na[1] - nb[1]) <= 2 * step:
                        adj.add(tuple(sorted((cid[a["iid"]], cid[b["iid"]]))))
        notes = []
        if sp == "Ude":
            nm = collections.Counter(re.search(r"NAME=(LOC\d+)", x["name"]).group(1) for x in v if "NAME=LOC" in x["name"])
            notes += ["%s on %d models (RefSeq-derived name)" % (k, c) for k, c in nm.items()]
        rows.append(dict(og=og, family=FAM[og], submitted_title=SUB[og], species=sp, salinity=sal[sp], L1=len(v),
                         L2=len(g2), L3=len(g3), L3_full=len(full), L3_fragment=len(g3) - len(full),
                         adjacent_or_consecutive_cluster_pairs=len(adj), ref_len=ref[og],
                         max_len=max(ln.values()), notes="; ".join(notes)))
        for x in v:
            e = em.get(x["iid"], {})
            mrows.append(dict(og=og, family=FAM[og], species=sp, salinity=sal[sp], iid=x["iid"], name=x["name"][:90],
                              len=x["len"], md5=x["md5"][:10], gene_id=gene_id(sp, x["name"])[:40],
                              L2=next(k for k, g in enumerate(g2) if x["iid"] in g) + 1, L3=cid[x["iid"]] + 1,
                              fragment=int(x["len"] < 0.5 * ref[og]), emapper=e.get("Preferred_name", ""),
                              pfam=e.get("PFAMs", "")))
    if write:
        for fn, rr in (("fig3b_loci.tsv", rows), ("fig3b_loci_members.tsv", mrows)):
            with open(OUT + "/" + fn, "w", encoding="utf8", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rr[0]), delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rr)
    return rows, mrows, pairs, ref


if __name__ == "__main__":
    rows, mrows, pairs, ref = main()
    print("reference length (median of per-species longest, bivalves):", ref)
    for og in FAM:
        r = [x for x in rows if x["og"] == og]
        print(og, FAM[og], "species", len(r), "L1/L2/L3/L3_full totals", [sum(x[k] for x in r) for k in ("L1", "L2", "L3", "L3_full")])
