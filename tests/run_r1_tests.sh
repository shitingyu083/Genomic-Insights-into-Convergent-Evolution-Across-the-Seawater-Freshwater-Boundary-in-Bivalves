#!/bin/bash
# R1 reproduction test (Linux). Runs every downstream script of R1_reanalysis/ on a scratch copy of the
# repository (outputs overwrite the shipped files in the copy; compare_outputs.py then diffs against the pristine copy).
# Usage: bash run_r1_tests.sh <scratch_repo_dir>
set -u
REPO=$(cd "$1" && pwd)
export IZ_ROOT="$REPO/R1_reanalysis"
W="$IZ_ROOT/_R1_work"; A="$W/analysis"
PY=${PY:-python3}
RS=${RS:-Rscript}
export MPLBACKEND=Agg PYTHONIOENCODING=utf-8 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-4} MKL_NUM_THREADS=${MKL_NUM_THREADS:-4}
ONLY=${ONLY:-}   # optional space-separated step ids; others are skipped and status.tsv is appended to
LOG="$REPO/_test_logs"; mkdir -p "$LOG"; [ -z "$ONLY" ] && : > "$LOG/status.tsv"
# compressed copies -> plain names expected by scripts (DATA_README.md)
find "$IZ_ROOT" \( -name "Orthogroups.tsv.gz" -o -name "deseq2_all_vs_NS.tsv.gz" -o -name "vst.tsv.gz" \) -exec gunzip -kf {} \;

run() {  # run <id> <cwd> <cmd...>
  local id=$1 d=$2; shift 2
  if [ -n "$ONLY" ] && [[ " $ONLY " != *" $id "* ]]; then return 0; fi
  local t0=$(date +%s)
  (cd "$d" && timeout 3600 "$@") > "$LOG/$id.log" 2>&1
  local rc=$?
  printf "%s\t%s\t%ss\t%s\n" "$id" "$rc" "$(( $(date +%s) - t0 ))" "$*" >> "$LOG/status.tsv"
  echo "[$id] rc=$rc"
}

# Fsite screen (Table S6) and unit criterion
run fsite_scan     "$A/r2m1_fsite" $PY scan_reproduce.py
run fsite_unit     "$A/r2m1_fsite" $PY unit_criterion.py
run fsite_s6steps  "$A/r2m1_fsite" $PY s6_steps.py
# convergence null model (v4, reported) and stem-specific check
mkdir -p "$REPO/_test_out/A_convergence_v4"
R1_LP="$W/server_data/asr_v4" R1_OUT="$REPO/_test_out/A_convergence_v4" run conv_A "$W/reanalysis" $PY A_convergence.py
R1_OUT="$REPO/_test_out/A_convergence_v4" run conv_A8 "$W/reanalysis" $PY A8_stem_specific.py
# SSLS statistics, enrichment
run ssls_stats     "$A/ssls" $PY ssls_stats.py
run enrich         "$A/enrich" $PY enrich_candidates.py
# Table S1 traits -> PGLS of assembly span
run s1_parse       "$A/pgls_span" $PY parse_S1_R1.py
T="$REPO/_test_out/pgls_span_nov1"; mkdir -p "$T"
cp "$W/server_data/x/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/SpeciesTree_rooted_node_labels.nwk" "$T/"
cp "$A/iqtree_rerun/out/fixed_topology/iqtree_species_tree_codes.nwk" "$A/pgls_span/traits_S1_parsed.tsv" "$T/"
R1_WORK="$REPO/_test_out" run pgls_span_ref  "$T" $RS "$A/pgls_span/pgls_span_nov1.R"
if [ -z "$ONLY" ] || [[ " $ONLY " == *" pgls_span_ref "* ]]; then
  cp "$T/pgls_span_results.tsv" "$T/ref_pgls_span_results.tsv"; cp "$T/pgls_span_results_rawMb.tsv" "$T/ref_pgls_span_results_rawMb.tsv"
fi
R1_WORK="$REPO/_test_out" run pgls_span_iq   "$T" $RS "$A/pgls_span/pgls_span_nov1_iqtree.R"
run span_lambda    "$A/pgls_span" $RS span_lambda_curve.R
run fig4a_check    "$A/pgls_span" $RS fig4a_p_check.R
# IQ-TREE rerun support and topology
run node_support   "$A/iqtree_rerun" $RS node_support.R
# gene families (Fig. 3a, Table S4)
run fig3b_loci     "$A/fig3b" $PY fig3b_loci.py
run gf20_counts    "$A/genefam20" $PY genefam20_counts.py
R1_TREE="$A/iqtree_rerun/out/fixed_topology/iqtree_species_tree_codes.nwk" R1_OUTW="$A/genefam20/out" \
  run gf20_pgls    "$A/genefam20" $RS genefam20_pgls.R
run gf20_summary   "$A/genefam20" $PY summarise_pgls20.py
R1_TREE="$A/iqtree_rerun/out/fixed_topology/iqtree_species_tree_codes.nwk" R1_OUTW="$REPO/_test_out/busco_tmp" \
  run busco_sens   "$A/genefam20/busco_rerun_sens" $RS busco_rerun_sens.R
# amino-acid composition (Fig. 5, Table S5)
run aa_prep        "$W" $PY analysis/aa/aa_prep.py
run aa_pgls        "$A/aa" $RS aa_pgls.R
run aa_gc_check    "$A/aa" $RS aa_gc_check.R
run aa_gc_all20    "$A/aa" $RS aa_gc_all20.R
run aa_tables      "$A/aa" $PY aa_report_tables.py
# supplementary tables
for t in S4 S6 S7 S8; do run table_$t "$A/tables" $PY build_$t.py; done
run annotate_ogs   "$A" $PY annotate_ogs.py
echo "done"; cat "$LOG/status.tsv"
