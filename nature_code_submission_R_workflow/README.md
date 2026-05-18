# Nature Code Submission R Workflow

This directory contains a reproducible R workflow for the manuscript analyses associated with freshwater colonization in bivalves. It is designed for submission as Supplementary Software or for archiving with a DOI.

## Scope

The workflow recomputes the manuscript's R-side analyses from standardized input tables:

- genome sampling and quality-control summaries;
- Fsite recurrence summaries across four phylogenetically collapsed freshwater units;
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
nature_code_submission_R_workflow/
  run_all.R
  R/
  config/
  data/
    templates/
  external/
  outputs/
```

Place manuscript input tables in `data/`. If a table is absent, `run_all.R` falls back to the corresponding file in `data/templates/`, which is useful for checking the workflow structure but is not intended to reproduce the manuscript results.

## Required Input Tables

Use these filenames in `data/`:

- `species_metadata.csv`
- `genome_qc_metrics.csv`
- `ancestral_states.csv`
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
- `public_transcriptome_expression.csv`
- `public_transcriptome_deseq2.csv`
- `atg7_foldx_results.csv`

Column definitions are shown in `data/templates/`.

## Running

From inside `nature_code_submission_R_workflow/`:

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
- `candidate parallel-substitution signals` are recurrent candidate states across informative collapsed freshwater units.
- TE Kimura divergence is used as a relative age proxy, not as an absolute insertion-time estimate.
- Salinity-transfer analyses in *Mercenaria mercenaria* test acute low-salinity responsiveness, not evolutionary freshwater adaptation.

## Reproducibility

The workflow writes `outputs/sessionInfo.txt` at the end of `run_all.R`. For archival submission, include:

- this whole directory;
- all input tables in `data/`;
- upstream command logs and parameter files from external software;
- final Source Data files used for figures.
