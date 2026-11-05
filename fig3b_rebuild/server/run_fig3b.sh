#!/bin/bash
# Rebuilt Fig. 3b gene trees on the Pfam-defined ATG8, BECN1 and AQP families (R1 revision, 2026-10-30).
# One tip per species x family x L2 gene (longest protein); MAFFT L-INS-i; IQ-TREE 2.4.0 ModelFinder + 1,000 UFBoot.
set -u
T=${R1_WORK:-$HOME/R1_revision_20261024}/fig3b_trees
R="${R0_DATA:-$HOME/###qutongjinhua_0710cxcf}"
WD="$R/ortho_clean/orthofinder_output/Results_Dec31/WorkingDirectory"
IQ=~/mambaforge/bin/iqtree2
cd "$T" || exit 1
{ echo "start $(date -Is)"; mafft --version 2>&1; $IQ --version | head -1; } > versions.txt 2>&1
python3 extract.py "$WD/SequenceIDs.txt" ${R1_WORK:-$HOME/R1_revision_20261024}/aa_v2/pfam_seqs.fa.gz fig3b_tips.tsv "$R" > extract.log 2>&1 || exit 1
for fam in ATG8 BECN1 AQP; do
  th=8; [ "$fam" = AQP ] && th=16
  ( mafft --localpair --maxiterate 1000 --thread "$th" "$fam.fa" > "$fam.aln.fa" 2> "$fam.mafft.log"
    $IQ -s "$fam.aln.fa" -m MFP -B 1000 -T "$th" --seed 20261030 --prefix "$fam" > "$fam.iq.stdout" 2>&1
    echo "$fam $? $(date -Is)" >> iq_exit.txt ) &
done
wait
echo "end $(date -Is)" >> versions.txt
touch FIG3B.DONE
