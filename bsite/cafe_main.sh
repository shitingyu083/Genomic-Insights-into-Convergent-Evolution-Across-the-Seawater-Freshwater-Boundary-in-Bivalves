#!/bin/bash
# R1 final CAFE5 run (bivalve-rooted). Lgi (gastropod outgroup) pruned so that the CAFE5 root filter keeps families
# present at the bivalve MRCA, matching the N1 (bivalve-node) universe of the branch-site analysis.
# Stage 1: families with copy-number range <= 20 -> estimate lambda + global error model (-e).
# Stage 2: remaining families (range > 20, and >= 100-copy families) -> lambda and error model fixed from stage 1.
set -e
cd ${R1_WORK:-$HOME/R1_revision_20261024}/bsite/cafe
E=${R1_WORK:-$HOME/R1_revision_20261024}/tools/cafe_env/bin
$E/Rscript -e 'library(ape); t <- read.tree("ultra.tre"); t <- drop.tip(t, "Lgi"); stopifnot(is.ultrametric(t, tol=1e-6)); write.tree(t, "ultra_biv.tre", digits=8); cat(Ntip(t), max(node.depth.edgelength(t)), "\n")'
mkdir -p main
for f in counts_small counts_large; do
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Lgi") c=i} {s=""; for(i=1;i<=NF;i++) if(i!=c) s=s (s==""?"":"\t") $i; print s}' $f.tsv > main/$f.biv.tsv
done
awk -F'\t' 'NR==1{print > "main/s1.tsv"; print > "main/s2.tsv"; next}
  {mx=0;mn=1e9;for(i=3;i<=NF;i++){if($i>mx)mx=$i;if($i<mn)mn=$i}
   if(mx-mn<=20) print > "main/s1.tsv"; else print > "main/s2.tsv"}' main/counts_small.biv.tsv
tail -n +2 main/counts_large.biv.tsv >> main/s2.tsv
wc -l main/s1.tsv main/s2.tsv
$E/cafe5 -i main/s1.tsv -t ultra_biv.tre -e -c 24 -o main/s1 > main/s1.log 2>&1
L=$(awk '/^Lambda:/{print $2}' main/s1/Base_results.txt)
echo "stage1 lambda $L"
$E/cafe5 -i main/s2.tsv -t ultra_biv.tre -l $L -emain/s1/Base_error_model.txt -c 24 -o main/s2 > main/s2.log 2>&1
touch main/CAFE.DONE
