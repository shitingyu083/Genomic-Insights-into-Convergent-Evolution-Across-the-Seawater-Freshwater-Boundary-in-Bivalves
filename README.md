# Code Submission R Workflow

This directory contains a reproducible R workflow for the manuscript analyses associated with freshwater colonization in bivalves. It is designed for submission as Supplementary Software or for archiving with a DOI.

## Scope

The workflow recomputes the manuscript's R-side analyses from standardized input tables:

- genome sampling and quality-control summaries;
- Fsite screen (derived residues shared by at least three freshwater target taxa, filtered by stem origin and absence from non-freshwater taxa);
- branch-site positive-selection signal summaries from PAML outputs;
- candidate parallel-substitution signal summaries;
- site-wise log-likelihood support summaries;
- pathway enrichment with Benjamini-Hochberg correction;
- gene-family copy-number association tables;
- genome architecture association tables;
- TE-family, local TE-density and Kimura-divergence proxy summaries;
- amino-acid centred log-ratio composition analyses;
- salinity-transfer assay summaries and tests;
- public transcriptome heatmap and DESeq2 summary tables;
- ATG7 FoldX ddG summary tables;
- Source Data index for figures.

External programs used upstream (BUSCO, OrthoMCL, IQ-TREE, PAML, RAxML, RepeatModeler2, RepeatMasker, LTR_retriever, DESeq2, ColabFold and FoldX) are documented in `external/external_software_commands.md`. Their parsed outputs are the standardized inputs consumed by this R workflow.

## Directory Layout

```text
./
  run_all.R
  R/
  config/
  data/
    templates/
  external/
  outputs/
  R1_reanalysis/
```

`METHOD_CODE_MAP.csv` provides a manuscript-component-to-code mapping table for reviewers.

Place manuscript input tables in `data/`. If a table is absent, `run_all.R` falls back to the corresponding file in `data/templates/`, which is useful for checking the workflow structure but is not intended to reproduce the manuscript results. For journal submission, the final archived package should include the real input tables in `data/`, not only the templates.

For manuscript reproduction, replace every template with the archived analysis table generated from the final data freeze. The templates are schema examples only. The package is considered complete for reproduction only when all real input tables listed below are present in `data/` and the upstream command logs are archived.

## Required Input Tables

Use these filenames in `data/`:

- `species_metadata.csv`
- `species_tree.nwk`
- `genome_qc_metrics.csv`
- `ancestral_states.csv`
- `fsite_origin_branches.csv`
- `extant_site_states.csv`
- `site_annotations.csv`
- `branch_site_results.csv`
- `parallel_substitution_results.csv`
- `ssls_results.csv`
- `gene_pathway_annotations.csv`
- `gene_family_copy_number.csv`
- `genome_architecture_traits.csv`
- `te_family_abundance.csv`
- `local_te_density.csv`
- `te_kimura_bins.csv`
- `amino_acid_counts.csv`
- `salinity_assay_measurements.csv`
- `salinity_comparisons.csv`
- `public_transcriptome_expression.csv`
- `public_transcriptome_deseq2.csv`
- `atg7_foldx_results.csv`

Column definitions are shown in `data/templates/`.

## Data in this repository

The genomic data used by the R1 reanalysis are included under `R1_reanalysis/` (about 176 MB, no file above 50 MB).
`DATA_README.md` describes what is included and what is not, and `DATA_MANIFEST.csv` lists each file with its SHA-256 and the scripts that read it.
The salinity-transfer assay measurements are not included; `data/templates/` shows their schema.

`tests/run_r1_tests.sh` re-runs the downstream R1 scripts on a scratch copy of this repository, and `tests/compare_outputs.py` compares the regenerated outputs with the shipped ones (see `QA_SUMMARY.txt`):

    cp -r <repo> /tmp/repo_run && bash tests/run_r1_tests.sh /tmp/repo_run
    python3 tests/compare_outputs.py <repo> /tmp/repo_run > compare.tsv
    python3 tests/compare_xlsx.py <repo> /tmp/repo_run > compare_xlsx.tsv

## R Dependencies

The workflow uses base R for table handling and requires these R packages for the phylogenetically informed analyses reported in the manuscript:

- `ape`
- `nlme`
- `phytools`

If these packages are missing, the scripts stop for analyses that require phylogenetic correction rather than silently replacing them with non-phylogenetic tests.

## Running

From the repository root:

```bash
Rscript run_all.R
```

Main outputs are written to:

```text
outputs/tables/
outputs/sessionInfo.txt
```

## Important Interpretation Notes

The script keeps the manuscript's conservative language:

- `Fsite` is an ancestral-state-guided candidate-site filter, not an independent proof of adaptation.
- `branch-site positive-selection signals` are PAML screening signals conditional on the tested model and alignment quality.
- `candidate parallel-substitution signals` are Fsite records that pass both filters; sharing by two or more freshwater units is a descriptor.
- TE Kimura divergence is used as a relative age proxy, not as an absolute insertion-time estimate.
- Salinity-transfer analyses in *Mercenaria mercenaria* test acute low-salinity responsiveness, not evolutionary freshwater adaptation.

## Audit Documents

- `METHOD_CODE_MAP.csv`: maps Methods/Results components to scripts, inputs, outputs and upstream tools.
- `POLICY_AND_METHOD_REFERENCES.md`: summarizes the policy and method basis for this package.
- `QA_SUMMARY.txt`: records static checks and local execution status.
- `CHANGELOG_R1.md`: changes made for the R1 revision.

## Reproducibility

The workflow writes `outputs/sessionInfo.txt` at the end of `run_all.R`. For archival submission, include:

- this whole directory;
- all input tables in `data/`;
- upstream command logs and parameter files from external software;
- final Source Data files used for figures.

## R1 revision (INZ-2026-07-492, November 2026)

This release matches the revised manuscript. Changes to the R workflow:

- `R/02_fsite_branch_parallel_ssls.R` implements the R1 Fsite definition. A record is one derived residue at one alignment column. The candidate set contains residues that differ from the reconstructed state at the most recent common ancestor of the freshwater taxa (FW-MRCA) and are shared by at least three of the seven freshwater target taxa (`Awo;Cfl;Dpo;Dro;Lfo;Pst;Ude`; see `config/analysis_config.csv`). Filter 1 keeps residues that arose on a freshwater-unit stem (`fsite_origin_branches.csv`, written by `R1_reanalysis/_R1_work/analysis/r2m1_fsite/scan_reproduce.py`). Filter 2 keeps residues absent from all non-freshwater taxa, Lottia included. Sharing by two or more freshwater units is reported as a descriptor (`descriptor_two_or_more_units`) and is not a filter.
- `R/07_salinity_autophagy.R` uses two-sided Welch t tests. The comparisons and the adjustment are read from `salinity_comparisons.csv`: Fig. 6a (Atg7 qPCR) uses Bonferroni adjustment across its three pairwise comparisons and Fig. 6d (FAA) uses Benjamini-Hochberg adjustment across its two comparisons. Fig. 6c (flow cytometry) is summarised without a test.
- New inputs: `fsite_origin_branches.csv` and `salinity_comparisons.csv` (templates in `data/templates/`). `ancestral_states.csv` now holds the FW-MRCA state (`node = FW_MRCA`).

### R1 reanalysis scripts

`R1_reanalysis/_R1_work/` holds the Python, R and shell scripts used for the R1 reanalyses, in their original directory layout so that relative imports work. `R1_reanalysis/SCRIPT_INDEX.csv` lists each script with the first line of its header comment, which names the figure, table or check it serves. `METHOD_CODE_MAP.csv` links manuscript components to these scripts. Set `IZ_ROOT` to the absolute path of `R1_reanalysis/` (or run the scripts from that directory, the default) so that the scripts find their inputs under `_R1_work/`.

Absolute paths from the authors' machines were replaced by environment variables with the original layout as default:

| Variable | Meaning | Default |
|---|---|---|
| `IZ_ROOT` | directory that contains `_R1_work/` (here `R1_reanalysis/`) | `.` |
| `R1_WORK` | server working directory for the R1 runs | `$HOME/R1_revision_20261024` |
| `R0_DATA` | server directory with the original-submission alignments, trees and codeml outputs | `$HOME/###qutongjinhua_0710cxcf` |
| `CONDA_BASE` | conda installation used by the shell scripts | `$HOME/anaconda3` |

Python dependencies: numpy, pandas, scipy, matplotlib, openpyxl, pymupdf. R dependencies: ape, nlme, phytools, and DESeq2 for Fig. 6f.

`R1_reanalysis/_R1_work/fig6f_de/` re-quantifies the public *Mercenaria mercenaria* RNA-seq data (Fig. 6f) with salmon v1.4.0 against RefSeq GCF_014805675.1 and fits DESeq2 (`run2.sh`, `de.R`, `samples.tsv`, `mapping.txt`, `sessionInfo.txt`). The gene count table, the DESeq2 results and the VST matrix are included (`gene_counts.tsv`, `deseq2_all_vs_NS.tsv.gz`, `vst.tsv.gz`); the salmon `quant/` folders are regenerated from PRJNA975821 by `run2.sh`.

Not included: genome assemblies, alignments, trees, codeml outputs and assay measurements (inputs are described by the templates and in the manuscript's Data Availability statement); scripts that only audited the original-submission figure files.
