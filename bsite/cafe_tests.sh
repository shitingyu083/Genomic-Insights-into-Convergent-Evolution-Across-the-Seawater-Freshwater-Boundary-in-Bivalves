#!/bin/bash
# CAFE5 diagnostics: the -e run on counts_small gave -lnL = inf at every lambda. Test whether the inf comes from
# the error model or from families with large copy-number ranges.
cd ${R1_WORK:-$HOME/R1_revision_20261024}/bsite/cafe
E=${R1_WORK:-$HOME/R1_revision_20261024}/tools/cafe_env/bin
mkdir -p tests
for d in 10 20 40; do
  awk -F'\t' -v D=$d 'NR==1{print;next}{mx=0;mn=1e9;for(i=3;i<=NF;i++){if($i>mx)mx=$i;if($i<mn)mn=$i} if(mx-mn<=D)print}' \
    counts_small.tsv > tests/d$d.tsv
done
wc -l tests/d*.tsv > tests/sizes.txt
$E/cafe5 -i tests/d20.tsv -t ultra.tre -c 8 -o tests/d20_noe > tests/d20_noe.log 2>&1 &
$E/cafe5 -i tests/d20.tsv -t ultra.tre -e -c 8 -o tests/d20_e > tests/d20_e.log 2>&1 &
$E/cafe5 -i counts_small.tsv -t ultra.tre -c 8 -o tests/all_noe > tests/all_noe.log 2>&1 &
$E/cafe5 -i tests/d40.tsv -t ultra.tre -e -c 8 -o tests/d40_e > tests/d40_e.log 2>&1 &
wait
touch tests/TESTS.DONE
