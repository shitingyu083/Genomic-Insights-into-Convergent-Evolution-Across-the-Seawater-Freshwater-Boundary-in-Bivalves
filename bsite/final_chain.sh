#!/bin/bash
# Waits for branch-site run1 and CAFE5 two-stage main run, then builds the final three-line evidence tables.
cd ${R1_WORK:-$HOME/R1_revision_20261024}/bsite
while [ ! -f run1/BSITE.DONE ] || [ ! -f cafe/main/CAFE.DONE ]; do sleep 120; done
set -e
mkdir -p final_3l pcctrl
python3 cafe_fw.py cafe/main/s1 Base cafe/main/s1_fw
python3 cafe_fw.py cafe/main/s2 Base cafe/main/s2_fw
(cat cafe/main/s1_fw/fw_expansion.tsv; tail -n +2 cafe/main/s2_fw/fw_expansion.tsv) > cafe/main/fw_expansion.tsv
for k in 2 3 4; do python3 pc_control.py run1 pcctrl/final_K${k}_pp0.8.tsv $k 0.8 2>pcctrl/final_K${k}.err | tee pcctrl/final_K${k}.txt; done
python3 pc_control.py run1 pcctrl/final_K3_pp0.95.tsv 3 0.95 2>pcctrl/final_K3_95.err | tee pcctrl/final_K3_pp0.95.txt
python3 pc_strata.py run1 > pcctrl/final_strata.txt
python3 three_lines.py run1 kegg cafe/main/fw_expansion.tsv kegg_pathway_list.tsv final_3l pcctrl/final_K3_pp0.8.tsv > final_3l/stdout.txt 2>&1
touch final_3l/FINAL.DONE
