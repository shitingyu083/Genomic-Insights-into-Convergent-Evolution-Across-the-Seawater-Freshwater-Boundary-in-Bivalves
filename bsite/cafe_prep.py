"""R1 revision (Fig. 2a, third evidence line): gene-family input for CAFE5.

Counts genes (not protein models) per orthogroup and species from OrthoFinder Orthogroups.tsv:
  RefSeq species   isoforms collapsed: proteins sharing the product description once "ISOFORMX#" is removed
  GENE_RE species  proteins collapsed by annotation gene ID (same regexes as analysis/fig3b/fig3b_loci.py)
  other species    one protein = one gene
Filters (CAFE5 tutorial practice): present in >= 2 species including at least one non-Lgi bivalve;
families with < 100 genes in every species go to counts_small.tsv (lambda estimation), the rest to counts_large.tsv.
Output: cafe/counts_small.tsv, cafe/counts_large.tsv, cafe/collapse_summary.tsv
"""
import os, re, collections

HOME = os.path.expanduser("~")
RES = os.environ.get("R0_DATA", HOME + "/###qutongjinhua_0710cxcf") + "/ortho_clean/orthofinder_output/Results_Dec31"
OUT = os.environ.get("R1_WORK", HOME + "/R1_revision_20261024") + "/bsite/cafe"
os.makedirs(OUT, exist_ok=True)
REFSEQ = {"Can", "Cgi", "Cvi", "Dpo", "Mar", "Mca", "Mme", "Oed", "Rph", "Yba"}
GENE_RE = {
    "Cni": r"^(CNI_\d+)-R[A-Z]+", "Awo": r"^(AWO\w+?G\d+)\.M\d+", "Lfo": r"^(.+?)\.T\d+$", "Pfu": r"^(.+?)\.T\d+$",
    "Mco": r"EVM(LG\d+\.G\d+)", "Mve": r"EVM\.TU\.([A-Z0-9]+\.\d+)", "Pst": r"^(PST_\d+)-T\d+",
    "Sgl": r"^(SGL\d+)-MRNA\d+", "Sso": r"GENE:(ENSAIOG\d+)", "Ude": r"GENE=(UDEG\d+)", "Sbr": r"GENE=(EVM\d+)",
    "Cbi": r"^(.+?)LOCUS=", "Cho": r"^(.+?)LOCUS=", "Sol": r"^(.+?)LOCUS="}


def unit(sp, name):
    if sp in REFSEQ:
        m = re.match(r"^XP_\d+\.\d+(.*?)ISOFORMX\d+(\[.*)?$", name)
        return ("iso", m.group(1)) if m else name
    if sp in GENE_RE:
        m = re.search(GENE_RE[sp], name)
        return m.group(1) if m else name
    return name


f = open(RES + "/Orthogroups/Orthogroups.tsv")
hdr = f.readline().rstrip("\n").split("\t")
spp = [h.replace(".pep", "") for h in hdr[1:]]
prot = collections.Counter()
gene = collections.Counter()
small = open(OUT + "/counts_small.tsv", "w")
large = open(OUT + "/counts_large.tsv", "w")
for o in (small, large):
    o.write("Desc\tFamily ID\t" + "\t".join(spp) + "\n")
ns = nl = 0
for line in f:
    p = line.rstrip("\n").split("\t")
    og, cells = p[0], p[1:] + [""] * (len(spp) - len(p) + 1)
    cnt = []
    for sp, c in zip(spp, cells):
        names = [x for x in c.split(", ") if x]
        g = len({unit(sp, x) for x in names})
        prot[sp] += len(names)
        gene[sp] += g
        cnt.append(g)
    present = sum(1 for x in cnt if x > 0)
    if present < 2 or all(x == 0 for sp, x in zip(spp, cnt) if sp != "Lgi"):
        continue
    o = small if max(cnt) < 100 else large
    o.write("(null)\t" + og + "\t" + "\t".join(map(str, cnt)) + "\n")
    if o is small:
        ns += 1
    else:
        nl += 1
with open(OUT + "/collapse_summary.tsv", "w") as o:
    o.write("species\tproteins_in_OGs\tgenes_in_OGs\n")
    for sp in spp:
        o.write("%s\t%d\t%d\n" % (sp, prot[sp], gene[sp]))
print("small", ns, "large", nl)
