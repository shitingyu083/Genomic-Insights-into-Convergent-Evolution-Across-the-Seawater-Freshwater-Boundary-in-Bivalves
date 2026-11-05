# External Software Used Before the R Workflow

The R workflow in this directory recomputes the statistical summaries, false-discovery-rate correction, candidate-set overlaps, Source Data tables and figure-input tables. Several upstream operations in the manuscript use established command-line software and should be run before `run_all.R`.

Replace paths, thread counts and database locations with the archived values used for the manuscript. The resulting tabular outputs should be copied into `data/` using the filenames listed in `README.md`.

## BUSCO

```bash
busco -i genome.fa -l metazoa_odb10 -m genome -o BUSCO_species_genome --cpu 16
busco -i proteins.fa -l metazoa_odb10 -m proteins -o BUSCO_species_proteins --cpu 16
```

Required R input derived from BUSCO:

- `data/genome_qc_metrics.csv`

## BLASTP and OrthoMCL

```bash
makeblastdb -in all_proteins.fa -dbtype prot -out all_proteins
blastp -query all_proteins.fa -db all_proteins -evalue 1e-7 -seg yes -outfmt 6 -num_threads 32 -out all_vs_all.blastp.tsv
orthomclFilterFasta compliantFasta 10 20
orthomclBlastParser all_vs_all.blastp.tsv compliantFasta > similarSequences.txt
orthomclMclToGroups OG 1 mclOutput > Orthogroups.tsv
```

Required R input derived from orthogroups:

- `data/gene_family_copy_number.csv`
- site-level alignment/state tables used in `data/extant_site_states.csv`

## IQ-TREE and MCMCTree

```bash
iqtree2 -s concatenated_single_copy_alignment.phy -p partitions.nex -m MFP -B 1000 -alrt 1000 -T AUTO
mcmctree mcmctree.ctl
```

Required R input or manuscript source:

- `data/species_metadata.csv`
- species tree files archived with Supplementary Software or Source Data

## PAML codeml ancestral reconstruction and branch-site tests

```bash
codeml ancestral_reconstruction.ctl
codeml branch_site_null.ctl
codeml branch_site_alt.ctl
```

Required R input derived from PAML:

- `data/ancestral_states.csv`
- `data/branch_site_results.csv`

The R workflow reports both the conservative chi-square(1) P value and the 50:50 mixture P value for the modified branch-site test; the manuscript uses the conservative chi-square(1) screen unless otherwise stated.

## RAxML site-wise log-likelihood support

```bash
raxmlHPC-PTHREADS -f G -m PROTGAMMAJTT -s alignment.phy -z candidate_topologies.tre -n gene_sitewise -T 16
```

Required R input derived from RAxML:

- `data/ssls_results.csv`

## RepeatModeler2, RepeatMasker and LTR_retriever

```bash
RepeatModeler -database species_db -LTRStruct -pa 16
RepeatMasker -pa 16 -s -lib species_repeat_library.fa genome.fa
LTR_retriever -genome genome.fa -inharvest candidates.scn
```

Required R input derived from TE annotation:

- `data/te_family_abundance.csv`
- `data/local_te_density.csv`
- `data/te_kimura_bins.csv`

Kimura divergence from consensus is treated as a relative age proxy in the R workflow, not as an absolute insertion time.

## Public RNA-seq Processing

```bash
fastp -i reads_R1.fastq.gz -I reads_R2.fastq.gz -o clean_R1.fastq.gz -O clean_R2.fastq.gz
hisat2 -x genome_index -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz -S sample.sam
featureCounts -a annotation.gtf -o counts.txt *.bam
```

DESeq2 outputs should be exported as:

- `data/public_transcriptome_expression.csv`
- `data/public_transcriptome_deseq2.csv`

## AlphaFold/ColabFold and FoldX

```bash
colabfold_batch ATG7.fasta colabfold_out/
foldx --command=RepairPDB --pdb=ATG7_model.pdb
foldx --command=BuildModel --pdb=ATG7_model_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=5
```

Required R input:

- `data/atg7_foldx_results.csv`

## R1 revision: upstream programs and versions

The R1 reanalyses used the programs below. Exact command lines are in the shell and Python scripts under `R1_reanalysis/_R1_work/`.

- Species tree: IQ-TREE v2.4.0 on the concatenated single-copy supermatrix (`analysis/iqtree_rerun/server/run_iqtree.sh`, `run_iqtree_phase2.sh`); the OrthoFinder species tree is the reference topology for ancestral reconstruction (`analysis/iqtree_rerun/compare_topology.R`, `node_support.R`).
- Ancestral reconstruction and Fsite records: PAML v4.9j codeml (AAML, JTT+F) on 99 orthogroups (`remote/build_asr.py`, `server_data/asr_v4/run_all.sh`, `analysis/paml_parse.py`); Fsite scan reproduced by `analysis/r2m1_fsite/scan_reproduce.py`; Table S6 step counts by `analysis/r2m1_fsite/s6_steps.py` and `analysis/tables/build_S6.py`.
- Null-model tests of recurrent substitutions (JTT-fgene, JTT-fsite): `analysis/conv_core.py`, `reanalysis/A_convergence.py`, `reanalysis/A8_stem_specific.py`, `analysis/tables/build_S7.py`.
- Branch-site tests: PAML codeml branch-site model A (`analysis/bsite/bsite.py`, `analysis/bsite/final_chain.sh`).
- Gene-family expansion: CAFE5 on an ultrametric OrthoFinder tree (`analysis/bsite/ultra_tree.R`, `cafe_prep.py`, `cafe_main.sh`, `cafe_tests.sh`, `cafe_fw.py`).
- SSLS: RAxML v8.2.12 `-f G` (`analysis/ssls/results/final/ssls.sh`, `ssls_gene.py`, `ssls_pool.py`, `analysis/ssls/ssls_stats.py`).
- Fig. 3b gene trees: MAFFT L-INS-i and IQ-TREE 2.4.0 with ModelFinder and 1,000 ultrafast bootstraps (`analysis/fig3b_rebuild/server/extract.py`, `run_fig3b.sh`, `analysis/fig3b/fig3b_trees_R1.py`).
- Genome completeness for Table S1: BUSCO 5.7.1, protein mode, mollusca_odb10 (`analysis/tableS1/server/run_busco.sh`, `busco_collapse.py`).
- Functional annotation: orthogroup names from the original FASTA headers (`analysis/annotate_ogs.py`); eggNOG-mapper/KEGG annotation (`analysis/bsite/kegg_map.py`, `analysis/enrich/enrich_candidates.py`).
- PGLS (assembly span, gene families, amino-acid composition): R packages ape and nlme (`analysis/pgls_span/`, `analysis/genefam20/`, `analysis/aa/`).

## R1 revision: public RNA-seq (Fig. 6f)

The R1 analysis replaces the fastp/HISAT2/featureCounts template above with transcript quantification:

```bash
salmon index -t GCF_014805675.1_rna.fna -i salmon_idx
salmon quant -i salmon_idx -l A -1 reads_1.fastq -2 reads_2.fastq --validateMappings -o quant/<sample>
Rscript de.R
```

`R1_reanalysis/_R1_work/fig6f_de/run2.sh` downloads the RefSeq GCF_014805675.1 transcripts and annotation, fetches the 21 SRA runs listed in `samples.tsv`, and runs salmon v1.4.0. `de.R` sums transcript counts to genes, keeps genes with at least 10 counts in at least three samples, fits DESeq2 with design `~ trt` (reference NS) and writes the contrasts AH1, AH5, CH, AL1, AL5 and CL against NS.
