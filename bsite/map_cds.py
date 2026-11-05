#!/usr/bin/env python3
"""R1 revision (R2-M1/Fig 2a): map OrthoFinder proteins to coding sequences by exact translation.
Read-only on inputs; writes only <out_dir> under ~/R1_revision_20261024/.
For each species with a CDS file on the server, every CDS is translated (standard code, frame 1, terminal stop
removed). An OrthoFinder protein is mapped to a CDS when the translation is identical to the protein sequence
(protein 'X' and terminal '*' tolerated). Output: <out>/cds_map/<sp>.tsv (protein_key, cds_id, cds_seq) and
<out>/cds_map_summary.tsv (proteins, mapped, fraction).
usage: python3 map_cds.py <out_dir>
"""
import gzip, os, re, sys

HOME = os.path.expanduser("~")
R = os.environ.get("R0_DATA", HOME + "/###qutongjinhua_0710cxcf") + ""
WD = R + "/ortho_clean/orthofinder_output/Results_Dec31/WorkingDirectory"
CDS = {  # species -> candidate CDS files (all tried, union of matches)
    "Cfl": ["species_cds/cfl.cds.fa", "sty0721/Cfl/Cfl.cds.fasta"],
    "Dpo": ["species_cds/Dpo.cds.fa", "sty0721/Dpo/Dpo_cds.fasta"],
    "Dro": ["species_cds/Dro.cds.fa"],
    "Lfo": ["species_cds/lfo.cds.fa", "sty0721/Lfo/lfortunei_genome_cds.fa"],
    "Pst": ["species_cds/Pstr.cds.fa", "sty0721/Pst/Potamilus_streckersoni.genomic_cds.fa"],
    "Cgi": ["species_cds/Cgi.cds.fa", "sty0721/Cgi/Cgi_cds.fasta"],
    "Csi": ["species_cds/Csi.cds.fa"],
    "Tna": ["species_cds/Tna.cds.fa"],
    "Sco": ["species_cds/sco.cds.fa", "sty0721/Sco/Sinonovacula_constricta.genomic_cds.fa", "sty0721/Sco/Sco_proteins_cds.fasta"],
    "Sbr": ["sty0721/Sbr/Sbr_assembly_cds.fasta"],
    "Mme": ["species_cds/Mme.cds.fa"],
    "Mye": ["sty0721/Mye/Mye_cds.fa"],
    "Mco": ["sty0721/Mco/Mco.fna.TBtools.fa", "sty0721/Mco/Mco-exon.fna"],
}
B = "TCAG"
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TAB = {a + b + c: CODE[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def key(s):
    s = re.sub(r"[^\x20-\x7e]", "", s)
    return re.sub(r"[\s:,()]", "_", s).upper()


def tr(s):
    return "".join(TAB.get(s[i:i + 3], "X") for i in range(0, len(s) - len(s) % 3, 3))


def main():
    out = sys.argv[1]
    os.makedirs(out + "/cds_map", exist_ok=True)
    spid = {}
    for l in open(WD + "/SpeciesIDs.txt"):
        i, n = l.rstrip("\n").split(": ", 1)
        spid[n.replace(".pep.fa", "").replace(".fa", "")] = i
    seqid = {}
    with open(WD + "/SequenceIDs.txt", encoding="utf-8", errors="replace") as f:
        for l in f:
            l = l.rstrip("\n")
            if ": " in l:
                k, h = l.split(": ", 1)
                seqid[k] = h
    summ = open(out + "/cds_map_summary.tsv", "w")
    summ.write("species\tn_protein\tn_mapped\tfraction\tfiles\n")
    for sp, files in CDS.items():
        idx = {}
        used = []
        for rel in files:
            p = os.path.join(HOME, rel)
            if not os.path.exists(p):
                continue
            used.append(rel)
            for h, s in fasta(p):
                s = s.upper().replace("U", "T")
                if len(s) % 3 or len(s) < 150:
                    continue
                t = tr(s)
                if t.endswith("*"):
                    t, s = t[:-1], s[:-3]
                if "*" in t:
                    continue
                idx.setdefault(t, (h.split()[0], s))
        n = m = 0
        fo = open(out + "/cds_map/" + sp + ".tsv", "w")
        for h, s in fasta(WD + "/Species" + spid[sp] + ".fa"):
            n += 1
            s = s.upper().rstrip("*").replace("*", "X")
            hit = idx.get(s)
            if hit is None:
                continue
            m += 1
            fo.write("%s\t%s\t%s\n" % (key(seqid[h.split()[0]]), hit[0], hit[1]))
        fo.close()
        summ.write("%s\t%d\t%d\t%.4f\t%s\n" % (sp, n, m, m / max(n, 1), ";".join(used)))
        summ.flush()
    open(out + "/MAP.DONE", "w").close()


if __name__ == "__main__":
    main()
