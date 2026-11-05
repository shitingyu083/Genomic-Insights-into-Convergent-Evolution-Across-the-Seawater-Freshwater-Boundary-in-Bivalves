# Changelog

## R1 revision, 2026-11-05d (data and execution test)

- Added the genomic data used by the R1 reanalysis under `R1_reanalysis/` (`DATA_README.md`, `DATA_MANIFEST.csv`). Salinity-assay measurements remain templates.
- Added the Fig. 3b gene-tree files, the Fig. 3b eggNOG-mapper table and the IQ-TREE support files that the shipped scripts read but that were missing.
- Path strings in six Python scripts changed from Windows backslashes to forward slashes; `figS1_R1_abc.py` reads its bold font from `IZ_BOLD_FONT`; `annotate_ogs.py` creates its output folder.
- `R/05`, `R/07` and `R/08` skip tests on all-NA template columns instead of stopping, so `run_all.R` completes on the templates.
- `analysis/ssls/ssls_stats.py` reads its result files in sorted order, which makes the permutation P values independent of the operating system.
- Added `tests/run_r1_tests.sh`, `tests/compare_outputs.py` and `tests/compare_xlsx.py`; results in `QA_SUMMARY.txt`.
- The salinity-assay methods state that animals were randomly assigned to tanks with identical conditions and treated as independent experimental units.

## R1 revision, 2026-11-05 (INZ-2026-07-492)

- Fsite screen (`R/02_fsite_branch_parallel_ssls.R`) rewritten to the R1 definition: derived residue relative to the FW-MRCA state, shared by at least three of the seven freshwater target taxa; filter 1 = origin on a freshwater-unit stem; filter 2 = absent from all non-freshwater taxa. Sharing by two or more freshwater units is a descriptor. The R0 unit-recurrence and seawater-alternative thresholds were removed from `config/analysis_config.csv`.
- Salinity assays (`R/07_salinity_autophagy.R`): Kruskal-Wallis and Wilcoxon tests replaced by two-sided Welch t tests with per-panel adjustment (Fig. 6a Bonferroni over three comparisons; Fig. 6d Benjamini-Hochberg over two comparisons), configured in `salinity_comparisons.csv`.
- New templates: `fsite_origin_branches.csv`, `salinity_comparisons.csv`; `ancestral_states.csv` now holds FW-MRCA states; the assay template records a single sampling time point.
- Added `R1_reanalysis/` with the scripts of the R1 reanalyses (species tree, ancestral reconstruction, Fsite scan, null-model tests, branch-site tests, CAFE5, SSLS, Fig. 3b gene trees, Table S1 BUSCO rerun, PGLS of assembly span, gene families and amino-acid composition, Fig. 6f DESeq2, figure and Supplementary Table builders). Absolute paths were replaced by environment variables (see README).
- Removed the duplicate folder `nature_code_submission_R_workflow/` (an older copy of the root workflow).
- Updated `METHOD_CODE_MAP.csv`, `CODE_AVAILABILITY_TEXT.txt`, `QA_SUMMARY.txt`, `external/external_software_commands.md`, `CODE_FILE_MANIFEST.csv`; added `LICENSE`.

## Original submission, 2026-05-17

- R workflow (`run_all.R`, `R/00`-`R/09`), templates, configuration and documentation.
