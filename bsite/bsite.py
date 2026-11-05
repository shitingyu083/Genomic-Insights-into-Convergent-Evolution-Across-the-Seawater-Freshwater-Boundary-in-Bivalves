#!/usr/bin/env python3
"""R1 revision (R2-M1, Fig 2a): genome-wide branch-site positive-selection and parallel/convergent-substitution
screens on the bivalves with coding sequences, following the three-evidence-line logic of Yuan et al. (2020,
Mol. Ecol. 29:3954-3969): positively selected genes (PSGs), parallel/convergent substituted genes (PCGs) and
expanded gene families are identified separately on the same orthologue universe and then intersected.

Read-only on inputs; writes only <out_dir> under ~/R1_revision_20261024/bsite/.
Taxa: freshwater (foreground) Pst (Unionida), Dpo + Dro (Dreissenidae), Cfl (Corbicula), Lfo (Limnoperna);
comparators Cgi, Csi, Tna, Sco, Sbr (always) and Mye (when single-copy and CDS-mapped).
Universe: bivalve-node hierarchical orthogroups (OrthoFinder N1 HOGs) with 1-3 CDS-mapped copies in each of the
10 core taxa (Mye used when 1-3 mapped copies); multi-copy taxa reduced to the copy with the highest mean identity to
the single-copy taxa.
Per OG: MAFFT v7.490 (--auto) protein alignment, codons threaded on the alignment; a sequence missing > 30% of the
alignment core (columns occupied in >= 50% of taxa) is dropped (OG dropped if it is a core taxon); codon columns with a
gap, stop or ambiguous base in any retained taxon removed; >= 150 clean nt.
Tree: OrthoFinder rooted species tree pruned to the taxa present, unrooted; foreground (#1) = the four freshwater-unit
stem branches (Pst, Cfl, Lfo terminal branches; Dreissenidae stem).
Branch-site model A (model=2, NSsites=2, CodonFreq=2): alternative fix_omega=0 (omega 1.5) vs null fix_omega=1.
ASR: aaml (JTT, empirical+F, gamma 4) on the translated clean alignment, marginal reconstruction (RateAncestor=1).
A site is parallel when >= 2 freshwater-unit branches change from the same ancestral residue to the same derived
residue, convergent when they change from different ancestral residues to the same derived residue.
usage: python3 bsite.py <out_dir> [nproc] [--count]
"""
import os, re, subprocess, sys, time
from multiprocessing import Pool

HOME = os.path.expanduser("~")
R = os.environ.get("R0_DATA", HOME + "/###qutongjinhua_0710cxcf") + ""
RES = R + "/ortho_clean/orthofinder_output/Results_Dec31"
TREE = R + "/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/SpeciesTree_rooted_node_labels.nwk"
JONES = R + "/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/jones.dat"
MAP = os.environ.get("R1_WORK", HOME + "/R1_revision_20261024") + "/bsite/out/cds_map"
FW = ["Pst", "Dpo", "Dro", "Cfl", "Lfo"]
CORE = FW + ["Cgi", "Csi", "Tna", "Sco", "Sbr"]
OPT = ["Mye"]
UNITS = {"Unionida": ["Pst"], "Dreissenidae": ["Dpo", "Dro"], "Corbicula": ["Cfl"], "Limnoperna": ["Lfo"]}
B = "TCAG"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB = {a + b + c: CODE[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
BS = """seqfile = codon.phy
treefile = fg.tree
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = 2
NSsites = 2
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = %d
omega = %s
cleandata = 1
fix_blength = 0
method = 0
"""
AAML = """seqfile = aa.phy
treefile = plain.tree
outfile = out.txt
noisy = 0
verbose = 0
runmode = 0
seqtype = 2
clock = 0
aaDist = 0
model = 3
aaRatefile = jones.dat
fix_alpha = 0
alpha = 0.5
ncatG = 4
RateAncestor = 1
cleandata = 0
method = 1
"""


# ---------------- newick ----------------
class Node:
    def __init__(self, name=None):
        self.name, self.ch, self.fg = name, [], False


def parse(s):
    s = s.strip().rstrip(";")
    pos = [0]

    def rd():
        n = Node()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                n.ch.append(rd())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                pos[0] += 1  # ')'
                break
        st = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",)":
            pos[0] += 1
        lab = s[st:pos[0]].split(":")[0]
        if not n.ch:
            n.name = lab.replace(".pep", "")
        return n
    return rd()


def prune(n, keep):
    if not n.ch:
        return n if n.name in keep else None
    ch = [c for c in (prune(c, keep) for c in n.ch) if c is not None]
    if not ch:
        return None
    if len(ch) == 1:
        return ch[0]
    n.ch = ch
    return n


def tips(n):
    return [n.name] if not n.ch else [t for c in n.ch for t in tips(c)]


def unroot(n):
    if len(n.ch) == 2:
        a, b = n.ch
        if a.ch and not a.fg:
            n.ch = a.ch + [b]
        elif b.ch and not b.fg:
            n.ch = [a] + b.ch
    return n


def mark(n):
    present = set(tips(n))
    for unit, members in UNITS.items():
        mem = sorted(m for m in members if m in present)
        if not mem:
            continue
        target = None
        stack = [n]
        while stack:
            x = stack.pop()
            if sorted(tips(x)) == mem:
                target = x
            stack.extend(x.ch)
        if target is not None:
            target.fg = True


def nwk(n, fg=True):
    s = n.name if not n.ch else "(" + ",".join(nwk(c, fg) for c in n.ch) + ")"
    return s + (" #1" if fg and n.fg else "")


# ---------------- io ----------------
def load_map(sp):
    d = {}
    p = MAP + "/" + sp + ".tsv"
    for l in open(p):
        k, cid, s = l.rstrip("\n").split("\t")
        d[k] = s
    return d


def key(s):
    s = re.sub(r"[^\x20-\x7e]", "", s)
    return re.sub(r"[\s:,()]", "_", s).upper()


MAXCOPY = 3


def universe():
    """Hierarchical orthogroups at the bivalve node (N1): 1..MAXCOPY CDS-mapped copies in every core taxon;
    Mye used when 1..MAXCOPY mapped copies. genes[sp] = list of CDS (representative chosen later)."""
    cds = {sp: load_map(sp) for sp in CORE + OPT}
    f = open(RES + "/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", encoding="utf-8", errors="replace")
    hdr = [h.replace(".pep", "") for h in f.readline().rstrip("\n").split("\t")]
    col = {h: i for i, h in enumerate(hdr)}
    uni = []
    for l in f:
        x = l.rstrip("\n").split("\t")
        genes = {}
        ok = True
        for sp in CORE + OPT:
            v = x[col[sp]] if col[sp] < len(x) else ""
            g = [key(t) for t in v.split(", ") if t]
            m = [cds[sp][t] for t in g if t in cds[sp]]
            if 1 <= len(g) <= MAXCOPY and len(m) == len(g):
                genes[sp] = m
            elif sp in CORE:
                ok = False
                break
        if ok:
            uni.append((x[0].replace("N1.", ""), genes))
    return uni


def pick(genes, wd):
    """One representative per taxon: single-copy taxa kept; for multi-copy taxa, the copy with the highest mean
    identity (MAFFT alignment of all copies, non-gap positions) to the single-copy taxa's sequences."""
    if all(len(v) == 1 for v in genes.values()):
        return {sp: v[0] for sp, v in genes.items()}, 0
    with open(wd + "/cand.fa", "w") as f:
        for sp, v in genes.items():
            for i, s in enumerate(v):
                f.write(">%s|%d\n%s\n" % (sp, i, tr(s)))
    aln = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", "cand.fa"], cwd=wd,
                         capture_output=True, text=True, timeout=3600).stdout
    al, nm = {}, None
    for l in aln.split("\n"):
        if l.startswith(">"):
            nm = l[1:].strip()
            al[nm] = []
        elif nm:
            al[nm].append(l.strip().upper())
    al = {k: "".join(v) for k, v in al.items()}
    ref = [al["%s|0" % sp] for sp, v in genes.items() if len(v) == 1]

    def ident(a, b):
        n = m = 0
        for x, y in zip(a, b):
            if x != "-" and y != "-":
                n += 1
                m += x == y
        return m / n if n else 0.0
    out = {}
    for sp, v in genes.items():
        if len(v) == 1:
            out[sp] = v[0]
            continue
        sc = [sum(ident(al["%s|%d" % (sp, i)], r) for r in ref) / max(len(ref), 1) for i in range(len(v))]
        out[sp] = v[max(range(len(v)), key=lambda i: sc[i])]
    return out, sum(len(v) > 1 for v in genes.values())


def tr(s):
    return "".join(TAB.get(s[i:i + 3], "X") for i in range(0, len(s), 3))


def write_phy(path, d, order):
    with open(path, "w") as f:
        f.write(" %d %d\n" % (len(order), len(d[order[0]])))
        for sp in order:
            f.write("%s  %s\n" % (sp, d[sp]))


def lnl(path):
    for l in open(path):
        if l.startswith("lnL("):
            return float(l.split(")")[1].split(":")[1].split()[0])
    return None


def beb(path):
    txt = open(path).read()
    i = txt.find("Bayes Empirical Bayes (BEB)")
    if i < 0:
        return []
    j = txt.find("Positive sites for foreground lineages", i)
    out = []
    for l in txt[j:].split("\n")[1:]:
        m = re.match(r"^\s+(\d+)\s+(\S)\s+(\d\.\d+)", l)
        if m:
            out.append((int(m.group(1)), m.group(2), float(m.group(3))))
        elif l.strip().startswith("The grid") or (out and not l.strip()):
            break
    return out


def run(cmd, wd, timeout=7200):
    return subprocess.run(cmd, cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)


def rst_parse(wd, order):
    """Return parent map, node names (tip index -> species) and per-site best state + prob for internal nodes."""
    txt = open(wd + "/rst").read()
    m = re.search(r"tree with node labels for Rod Page's TreeView\n(.+);", txt)
    t = m.group(1)
    # tokens: "(1_Pst, 2_Dpo) 12 ,"
    stack, parent, cur = [], {}, []
    toks = re.findall(r"\(|\)\s*\d+|\d+_[A-Za-z]+|,", t)
    for tk in toks:
        if tk == "(":
            stack.append([])
        elif tk.startswith(")"):
            nid = int(tk[1:].strip())
            for c in stack.pop():
                parent[c] = nid
            if stack:
                stack[-1].append(nid)
            else:
                root = nid
        elif tk == ",":
            continue
        else:
            nid = int(tk.split("_")[0])
            stack[-1].append(nid)
    nt = len(order)
    i = txt.find("Prob of best state at each node, listed by site")
    sites = []
    for l in txt[i:].split("\n")[3:]:
        mm = re.match(r"^\s+(\d+)\s+\d+\s+(\S+):\s+(.*)$", l)
        if not mm:
            if sites:
                break
            continue
        states = re.findall(r"([A-Z\-\*])\(([\d\.]+)\)", mm.group(3))
        sites.append((mm.group(2), states))
    return parent, root, nt, sites


def one(item):
    og, genes, outdir = item
    wd = outdir + "/og/" + og
    os.makedirs(wd, exist_ok=True)
    res = {"og": og}
    try:
        genes, res["n_multicopy_taxa"] = pick(genes, wd)
        order = [sp for sp in CORE + OPT if sp in genes]
        with open(wd + "/prot.fa", "w") as f:
            for sp in order:
                f.write(">%s\n%s\n" % (sp, tr(genes[sp])))
        aln = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", "prot.fa"], cwd=wd,
                             capture_output=True, text=True, timeout=3600).stdout
        al, nm = {}, None
        for l in aln.split("\n"):
            if l.startswith(">"):
                nm = l[1:].strip()
                al[nm] = []
            elif nm:
                al[nm].append(l.strip().upper())
        al = {k: "".join(v) for k, v in al.items()}
        cod = {}
        for sp in order:
            s, k, c = genes[sp], 0, []
            for a in al[sp]:
                if a == "-":
                    c.append("---")
                else:
                    c.append(s[3 * k:3 * k + 3])
                    k += 1
            cod[sp] = c
        L0 = len(al[order[0]])
        # columns where at least half of the taxa have a residue define the comparable core of the alignment
        core_cols = [i for i in range(L0) if sum(al[sp][i] != "-" for sp in order) >= 0.5 * len(order)]
        L = len(core_cols)
        drop = [sp for sp in order if sum(al[sp][i] == "-" for i in core_cols) > 0.5 * L]
        order = [sp for sp in order if sp not in drop]
        res["dropped"] = ",".join(drop)
        # keep the OG if >= 3 of the 4 freshwater units and >= 4 comparators survive
        # (all-four-unit OGs are flagged by the 'units' column for a sensitivity subset)
        nunit = sum(any(x in order for x in m) for m in UNITS.values())
        if nunit < 3 or sum(sp not in FW for sp in order) < 4:
            res["status"] = "too_few_taxa:" + ",".join(drop)
            return res
        keepcol = [i for i in range(L0) if all(TAB.get(cod[sp][i], "*") != "*" for sp in order)]
        res["ncodon"] = len(keepcol)
        res["ntaxa"] = len(order)
        if 3 * len(keepcol) < 150:
            res["status"] = "short"
            return res
        cc = {sp: "".join(cod[sp][i] for i in keepcol) for sp in order}
        aa = {sp: "".join(TAB[cod[sp][i]] for i in keepcol) for sp in order}
        write_phy(wd + "/codon.phy", cc, order)
        write_phy(wd + "/aa.phy", aa, order)
        t = unroot(prune(parse(open(TREE).read()), set(order)))
        mark(t)
        open(wd + "/fg.tree", "w").write(nwk(t) + ";\n")
        open(wd + "/plain.tree", "w").write(nwk(t, False) + ";\n")
        # branch-site
        # alternative model from three omega starts (best lnL kept) to avoid local optima
        out = {}
        for tag, fix, om in (("alt", 0, "0.5"), ("alt2", 0, "1.5"), ("alt3", 0, "4"), ("null", 1, "1")):
            d = wd + "/" + tag
            os.makedirs(d, exist_ok=True)
            for f in ("codon.phy", "fg.tree"):
                os.symlink("../" + f, d + "/" + f) if not os.path.exists(d + "/" + f) else None
            open(d + "/bs.ctl", "w").write(BS % (fix, om))
            run(["codeml", "bs.ctl"], d)
            out[tag] = lnl(d + "/out.txt")
        best = max(("alt", "alt2", "alt3"), key=lambda k: out[k] if out[k] is not None else -1e18)
        res["lnL_alt"], res["lnL_null"] = out[best], out["null"]
        res["beb"] = beb(wd + "/" + best + "/out.txt")
        # ASR
        d = wd + "/asr"
        os.makedirs(d, exist_ok=True)
        for f in ("aa.phy", "plain.tree"):
            os.symlink("../" + f, d + "/" + f) if not os.path.exists(d + "/" + f) else None
        if not os.path.exists(d + "/jones.dat"):
            os.symlink(JONES, d + "/jones.dat")
        open(d + "/aaml.ctl", "w").write(AAML)
        run(["codeml", "aaml.ctl"], d)
        parent, root, nt, sites = rst_parse(d, order)
        idx = {sp: i + 1 for i, sp in enumerate(order)}
        # branch (node id) of each freshwater unit stem
        def mrca(members):
            anc = []
            for sp in members:
                a, x = [], idx[sp]
                while x in parent:
                    a.append(x)
                    x = parent[x]
                a.append(x)
                anc.append(a)
            common = [n for n in anc[0] if all(n in a for a in anc[1:])]
            return common[0]
        ub = {u: mrca([x for x in m if x in idx]) for u, m in UNITS.items() if any(x in idx for x in m)}
        res["units"] = len(ub)
        events = []
        for s, (data, states) in enumerate(sites, 1):
            def st(n):
                if n <= nt:
                    return data[n - 1], 1.0
                k = n - nt - 1
                return states[k][0], float(states[k][1])
            ch = {}
            for u, n in ub.items():
                a, pa = st(parent[n])
                b, pb = st(n)
                if a != b:
                    ch[u] = (a, b, min(pa, pb))
            us = sorted(ch)
            for i in range(len(us)):
                for j in range(i + 1, len(us)):
                    a1, b1, p1 = ch[us[i]]
                    a2, b2, p2 = ch[us[j]]
                    if b1 == b2:
                        comp = [data[idx[sp] - 1] for sp in order if sp not in FW]
                        events.append((s, us[i] + "-" + us[j], "parallel" if a1 == a2 else "convergent",
                                       a1 + ">" + b1 + "|" + a2 + ">" + b2, min(p1, p2), int(b1 not in comp)))
        res["events"] = events
        res["status"] = "ok"
    except Exception as e:  # noqa
        res["status"] = "error:" + type(e).__name__ + ":" + str(e)[:120]
    return res


def main():
    out = sys.argv[1]
    npr = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 48
    os.makedirs(out + "/og", exist_ok=True)
    uni = universe()
    with open(out + "/universe.tsv", "w") as f:
        f.write("og\tntaxa_cds\thas_Mye\n")
        for og, g in uni:
            f.write("%s\t%d\t%d\n" % (og, len(g), int("Mye" in g)))
    if len(sys.argv) > 3 and sys.argv[3].startswith("--test"):
        uni = uni[::max(1, len(uni) // int(sys.argv[3].split("=")[1]))][:int(sys.argv[3].split("=")[1])]
    print("universe", len(uni), flush=True)
    if "--count" in sys.argv:
        return
    done = set()
    rp = out + "/bsite_results.tsv"
    if os.path.exists(rp):
        done = {l.split("\t")[0] for l in open(rp)}
    new = not os.path.exists(rp)
    fo = open(rp, "a")
    fe = open(out + "/pc_events.tsv", "a")
    if new:
        fo.write("og\tstatus\tntaxa\tncodon\tlnL_alt\tlnL_null\tbeb\tn_events\tn_multicopy_taxa\tdropped\tunits\n")
        fe.write("og\tsite\tpair\ttype\tchange\tminpp\tabsent_in_comparators\n")
    todo = [(og, g, out) for og, g in uni if og not in done]
    t0 = time.time()
    with Pool(npr) as p:
        for k, r in enumerate(p.imap_unordered(one, todo, chunksize=1), 1):
            b = r.get("beb", [])
            fo.write("\t".join(str(x) for x in [r["og"], r.get("status"), r.get("ntaxa", ""), r.get("ncodon", ""),
                     r.get("lnL_alt", ""), r.get("lnL_null", ""),
                     ";".join("%d%s:%.3f" % x for x in b), len(r.get("events", [])),
                     r.get("n_multicopy_taxa", ""), r.get("dropped", ""), r.get("units", "")]) + "\n")
            for e in r.get("events", []):
                fe.write(r["og"] + "\t" + "\t".join(str(x) for x in e) + "\n")
            fo.flush()
            fe.flush()
            if k % 200 == 0:
                print(k, len(todo), round(time.time() - t0), flush=True)
    open(out + "/BSITE.DONE", "w").close()


if __name__ == "__main__":
    main()
