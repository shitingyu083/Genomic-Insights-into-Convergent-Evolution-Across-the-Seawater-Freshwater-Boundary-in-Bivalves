"""Write one FASTA per family for the rebuilt Fig. 3b gene trees (R1 revision, 2026-10-30).

Inputs: OrthoFinder SequenceIDs.txt (iid -> header), aa_v2/pfam_seqs.fa.gz (iid -> sequence),
fig3b_tips.tsv (one representative protein, the longest, per species x family x L2 gene), and the
original proteomes <rdir>/<sp>.pep.fa (the OrthoFinder inputs before whitespace removal and upper-casing).
Tip labels use the species code and the first token of the species' own protein header, never a transferred name.
"""
import csv
import gzip
import re
import sys

seqids, pfam, tips, rdir = sys.argv[1:5]
want = {r["iid"]: r for r in csv.DictReader(open(tips, encoding="utf-8"), delimiter="\t")}

hdr = {}
with open(seqids, encoding="utf-8", errors="replace") as f:
    for line in f:
        iid, _, h = line.rstrip("\n").partition(": ")
        if iid in want:
            hdr[iid] = h

seq, cur = {}, None
with gzip.open(pfam, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            iid = line[1:].split()[0]
            cur = iid if iid in want else None
            if cur:
                seq[cur] = []
        elif cur:
            seq[cur].append(line.strip())


def squash(x):
    return re.sub(r"\s", "", x).upper()


orig = {}
for sp in sorted({r["species"] for r in want.values()}):
    by_key, by_seq = {}, {}
    h, buf = None, []
    with open(rdir + "/" + sp + ".pep.fa", encoding="utf-8", errors="replace") as f:
        lines = f.read().split("\n")
    for line in lines + [">"]:
        if line.startswith(">"):
            if h is not None:
                by_key.setdefault(squash(h), h)
                by_seq.setdefault("".join(buf).upper().rstrip("*"), h)
            h, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    orig[sp] = (by_key, by_seq)

how = {}
for iid, r in want.items():
    by_key, by_seq = orig[r["species"]]
    k = squash(hdr.get(iid, ""))
    s = "".join(seq.get(iid, [])).upper().rstrip("*")
    if k in by_key:
        hdr[iid], how[iid] = by_key[k], "header"
    elif s in by_seq:
        hdr[iid], how[iid] = by_seq[s], "sequence"
    else:
        how[iid] = "orthofinder_header"
print("original header matched by", {v: sum(1 for x in how.values() if x == v) for v in sorted(set(how.values()))})

missing = [i for i in want if i not in seq or i not in hdr]
print("tips", len(want), "with sequence", len(seq), "with header", len(hdr), "missing", len(missing))
if missing:
    print("MISSING", " ".join(missing[:20]))
    sys.exit(1)

AA = set("ACDEFGHIKLMNPQRSTVWY")
used, fams = set(), {}
with open("fig3b_tip_labels.tsv", "w", encoding="utf-8", newline="") as m:
    w = csv.writer(m, delimiter="\t", lineterminator="\n")
    w.writerow(["family", "label", "species", "habitat", "iid", "header_first_token", "len_aa", "n_members_merged",
                "rep_full", "og", "n_nonstandard_to_X", "header_match"])
    for iid, r in sorted(want.items(), key=lambda kv: (kv[1]["family"], kv[1]["species"], kv[0])):
        tok = hdr[iid].split()[0] if hdr[iid].split() else iid
        lab = r["species"] + "_" + re.sub(r"[^A-Za-z0-9_.-]", "_", tok)
        if lab in used:
            lab = lab + "_" + iid
        used.add(lab)
        s = "".join(seq[iid]).upper().rstrip("*")
        n_x = sum(1 for c in s if c not in AA)
        s = "".join(c if c in AA else "X" for c in s)
        fams.setdefault(r["family"], []).append((lab, s))
        w.writerow([r["family"], lab, r["species"], r["habitat"], iid, tok, len(s), r["n_members_merged"],
                    r["rep_full"], r["og"], n_x, how[iid]])

for fam, items in fams.items():
    with open(fam + ".fa", "w", encoding="utf-8", newline="\n") as f:
        for lab, s in items:
            f.write(">" + lab + "\n" + s + "\n")
    print(fam, len(items))
