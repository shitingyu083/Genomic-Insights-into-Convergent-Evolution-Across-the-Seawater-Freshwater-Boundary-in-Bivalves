# Validation Report

Date: 2026-05-17

## Static Checks Completed

- Workflow directory created with a clean R modular structure.
- Input templates are provided for all manuscript analysis modules.
- No development-marker strings were detected in the workflow scripts or documentation.
- No Unicode replacement characters were detected.
- Workflow documentation uses ASCII-only text after normalization.
- The code distinguishes upstream external software outputs from R-side statistical summaries.

## Execution Status

The workflow could not be executed in the current local environment because `Rscript` was not available on PATH and no local `Rscript.exe` installation was found under the searched standard Windows locations.

Run the workflow on a machine with R installed:

```bash
cd nature_code_submission_R_workflow
Rscript run_all.R
```

or on Windows PowerShell:

```powershell
cd nature_code_submission_R_workflow
.\run_workflow_windows.ps1
```

## Expected Outputs

After a successful run, the workflow writes:

- `outputs/tables/source_data_index.csv`
- `outputs/sessionInfo.txt`
- source tables for Figs. 1-6 and Supplementary Fig. 1 where corresponding inputs are supplied
- analysis tables for Fsite, branch-site positive-selection signals, candidate parallel-substitution signals, enrichment, gene-family copy number, genome architecture, TE profiles, amino-acid CLR composition, salinity assays, public transcriptome summaries and ATG7 FoldX summaries
