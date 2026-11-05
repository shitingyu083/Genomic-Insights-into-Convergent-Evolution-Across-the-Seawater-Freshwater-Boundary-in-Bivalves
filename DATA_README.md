# Data included in this repository

All files are under `R1_reanalysis/`. Set `IZ_ROOT` to that folder before running the R1 scripts:

```bash
export IZ_ROOT="$(pwd)/R1_reanalysis"
```

`DATA_MANIFEST.csv` lists every data file with its size, SHA-256, origin and the scripts that read it.
Every file is below GitHub's 50 MB warning size, so the repository can be pushed without Git LFS.

## What is included

| Folder (under `R1_reanalysis/`) | Content |
|---|---|
| `_R1_work/server_data/x/.../Fsite_method/` | OrthoFinder single-copy alignments and PAML ancestral states used by the Fsite screen |
| `_R1_work/server_data/x/ortho_clean/.../Results_Dec31/` | OrthoFinder species tree, `Orthogroups.tsv.gz`, `Orthogroups.GeneCount.tsv` |
| `_R1_work/server_data/asr_v4/` | R1 ancestral-state reconstructions (input to `A_convergence.py`) |
| `_R1_work/server_data/r1rev/` | R1 server outputs: amino-acid composition, GC check, 20-family copy numbers, Fig. 3b gene trees, annotations |
| `_R1_work/server_data/headers.tsv.gz` | Original FASTA headers used to annotate orthogroups |
| `_R1_work/reanalysis/A_convergence_v4/` | Reported convergence reanalysis (observed vs expected sharing, site tests) |
| `_R1_work/analysis/*/` | Inputs and outputs of each R1 analysis (PGLS, branch-site, SSLS, enrichment, gene families, IQ-TREE rerun, tables) |
| `_R1_work/fig6f_de/` | Public *M. mercenaria* RNA-seq (PRJNA975821) quantified with salmon on GCF_014805675.1: gene counts, DESeq2 results, VST matrix |
| `SUBMIT_R1_revised_20261024/` | Supplementary tables S1, S2, S4 to S9 and the figure source-data tables as published |

Compressed files: `Orthogroups.tsv.gz`, `deseq2_all_vs_NS.tsv.gz` and `vst.tsv.gz` are stored compressed to keep each file small.
Scripts that open the plain-text name need them decompressed first:

```bash
find R1_reanalysis -name "Orthogroups.tsv.gz" -o -name "deseq2_all_vs_NS.tsv.gz" -o -name "vst.tsv.gz" | xargs gunzip -k
```

Absolute paths in the original files were replaced by placeholders: `$R1_WORK` (R1 server working folder), `$R0_DATA` (original-submission data folder), `$HOME` and `<IZ_ROOT>`.

## What is not included, and why

| Not included | Reason / where to get it |
|---|---|
| Genome assemblies, proteomes and BUSCO inputs | Public; the accession or repository of every retained genome is in Table S1 (`SUBMIT_R1_revised_20261024/Supplementary_Table_S1_*`) |
| OrthoFinder `WorkingDirectory` sequences | Regenerated from the public proteomes with the commands in `external/external_software_commands.md` |
| RefSeq annotation `gff.gz` for GCF_014805675.1 | Downloaded from NCBI by `fig6f_de/run2.sh` |
| Salmon `quant/` folders | Regenerated from PRJNA975821 by `fig6f_de/run2.sh`; `de.R` starts from them, and the shipped `gene_counts.tsv` is its output |
| Salinity-transfer assay measurements (Fig. 6) | Not deposited in this repository; `data/templates/` gives the table schema used by `run_all.R` |
| Simulated or placeholder tables (e.g. `brackish_as_marine_simulated.csv`, pseudo-motif controls) | Development checks only; no manuscript number depends on them |
| TE annotation and FoldX output files | Large intermediate files; regenerated with the commands in `external/external_software_commands.md` |
| Run logs, backups and internal QA files | Not needed to reproduce results |
