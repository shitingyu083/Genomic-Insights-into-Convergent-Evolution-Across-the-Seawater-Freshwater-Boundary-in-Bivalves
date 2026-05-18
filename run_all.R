source(file.path("R", "00_utils.R"))
source(file.path("R", "01_sampling_qc.R"))
source(file.path("R", "02_fsite_branch_parallel_ssls.R"))
source(file.path("R", "03_enrichment.R"))
source(file.path("R", "04_gene_family_copy_number.R"))
source(file.path("R", "05_genome_architecture_te.R"))
source(file.path("R", "06_amino_acid_clr.R"))
source(file.path("R", "07_salinity_autophagy.R"))
source(file.path("R", "08_structure_foldx.R"))
source(file.path("R", "09_source_data_index.R"))

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "data")
template_dir <- file.path(data_dir, "templates")
outdir <- file.path(root, "outputs", "tables")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

config <- read_config(file.path(root, "config", "analysis_config.csv"))
set.seed(as.integer(config[["random_seed"]]))

path_or_template <- function(filename) {
  primary <- file.path(data_dir, filename)
  if (file.exists(primary)) return(primary)
  file.path(template_dir, filename)
}

paths <- list(
  species_metadata = path_or_template("species_metadata.csv"),
  species_tree = path_or_template("species_tree.nwk"),
  genome_qc = path_or_template("genome_qc_metrics.csv"),
  ancestral_states = path_or_template("ancestral_states.csv"),
  extant_site_states = path_or_template("extant_site_states.csv"),
  site_annotations = path_or_template("site_annotations.csv"),
  branch_site_results = path_or_template("branch_site_results.csv"),
  parallel_substitution_results = path_or_template("parallel_substitution_results.csv"),
  ssls_results = path_or_template("ssls_results.csv"),
  gene_pathway_annotations = path_or_template("gene_pathway_annotations.csv"),
  gene_family_copy_number = path_or_template("gene_family_copy_number.csv"),
  genome_architecture_traits = path_or_template("genome_architecture_traits.csv"),
  te_family_abundance = path_or_template("te_family_abundance.csv"),
  local_te_density = path_or_template("local_te_density.csv"),
  te_kimura_bins = path_or_template("te_kimura_bins.csv"),
  amino_acid_counts = path_or_template("amino_acid_counts.csv"),
  salinity_assay_measurements = path_or_template("salinity_assay_measurements.csv"),
  public_transcriptome_expression = path_or_template("public_transcriptome_expression.csv"),
  public_transcriptome_deseq2 = path_or_template("public_transcriptome_deseq2.csv"),
  atg7_foldx_results = path_or_template("atg7_foldx_results.csv")
)

message_step("1/9 sampling and genome QC")
run_sampling_qc(paths, outdir)

message_step("2/9 Fsite, branch-site positive-selection signals, candidate parallel substitutions and SSLS")
run_site_level_workflow(paths, outdir, config)

message_step("3/9 pathway enrichment")
run_enrichment(paths, outdir, q_threshold = as_num(config[["enrichment_q_threshold"]]))

message_step("4/9 gene-family copy number")
run_gene_family_copy_number(paths, outdir, config)

message_step("5/9 genome architecture")
run_genome_architecture(paths, outdir)

message_step("6/9 TE summaries and Kimura divergence proxy")
run_te_summaries(paths, outdir)

message_step("7/9 amino-acid CLR composition")
run_amino_acid_clr(paths, outdir)

message_step("8/9 salinity assays and public transcriptome")
run_salinity_assays(paths, outdir)

message_step("9/9 ATG7 FoldX/structure summaries and output index")
run_structure_foldx(paths, outdir)
run_source_data_index(outdir)

session <- utils::capture.output(utils::sessionInfo())
writeLines(session, file.path(root, "outputs", "sessionInfo.txt"))

message_step("complete. Tables written to ", outdir)
