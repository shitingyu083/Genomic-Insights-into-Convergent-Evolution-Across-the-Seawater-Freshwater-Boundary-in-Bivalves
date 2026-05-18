source(file.path("R", "00_utils.R"))

fit_trait_model <- function(dat, trait) {
  pgls <- fit_pgls_habitat(dat, trait, attr(dat, "species_tree"), species_col = "species_id")
  data.frame(
    trait = trait,
    estimate = pgls$estimate_freshwater_minus_seawater_group,
    standard_error = pgls$standard_error,
    statistic = pgls$statistic,
    p = pgls$p,
    model = pgls$model,
    lambda = pgls$lambda
  )
}

run_genome_architecture <- function(paths, outdir) {
  metadata <- standardize_species_metadata(read_table_auto(paths$species_metadata))
  tree <- read_species_tree(paths$species_tree, required = TRUE)
  traits <- read_table_auto(paths$genome_architecture_traits)
  require_cols(traits, c("species_id"), "genome_architecture_traits")
  x <- merge(traits, metadata[, c("species_id", "habitat_binary")], by = "species_id", all.x = TRUE)
  x <- subset(x, habitat_binary %in% c("freshwater", "seawater_group"))
  trait_cols <- setdiff(names(traits), c("species_id", "busco_complete", "busco_duplicated",
                                         "scaffold_n50", "annotated_gene_number"))
  trait_cols <- trait_cols[vapply(traits[trait_cols], function(z) any(!is.na(as_num(z))), logical(1))]
  if (!length(trait_cols)) {
    empty <- data.frame(
      trait = character(), estimate = numeric(), statistic = numeric(),
      p = numeric(), model = character(), q = numeric()
    )
    write_table(empty, file.path(outdir, "table_genome_architecture_associations.csv"))
    write_table(x, file.path(outdir, "source_fig4a_genome_architecture.csv"))
    return(invisible(empty))
  }
  for (cc in trait_cols) x[[cc]] <- as_num(x[[cc]])
  attr(x, "species_tree") <- tree
  res <- do.call(rbind, lapply(trait_cols, function(tr) fit_trait_model(x, tr)))
  res$q <- bh(res$p)
  write_table(res, file.path(outdir, "table_genome_architecture_associations.csv"))
  write_table(x, file.path(outdir, "source_fig4a_genome_architecture.csv"))
  invisible(res)
}

run_te_summaries <- function(paths, outdir) {
  fam <- read_table_auto(paths$te_family_abundance, required = FALSE)
  if (!is.null(fam)) {
    require_cols(fam, c("species_id", "te_family", "te_bp", "assembly_span_bp"), "te_family_abundance")
    fam$te_bp <- as_num(fam$te_bp)
    fam$assembly_span_bp <- as_num(fam$assembly_span_bp)
    fam$te_fraction <- fam$te_bp / fam$assembly_span_bp
    fam$te_fraction_scaled_by_family <- ave(fam$te_fraction, fam$te_family, FUN = scale_within_group)
    write_table(fam, file.path(outdir, "source_fig4b_te_family_profiles.csv"))
  }

  local <- read_table_auto(paths$local_te_density, required = FALSE)
  if (!is.null(local)) {
    require_cols(local, c("species_id", "gene", "region", "te_family", "te_bp", "region_length_bp"), "local_te_density")
    local$te_bp <- as_num(local$te_bp)
    local$region_length_bp <- as_num(local$region_length_bp)
    local$te_density <- local$te_bp / local$region_length_bp
    local$te_density_scaled <- ave(local$te_density, paste(local$gene, local$region, local$te_family), FUN = scale_within_group)
    write_table(local, file.path(outdir, "source_fig4c_local_te_density.csv"))
  }

  kim <- read_table_auto(paths$te_kimura_bins, required = FALSE)
  if (!is.null(kim)) {
    require_cols(kim, c("species_id", "kimura_divergence_percent", "te_bp", "assembly_span_bp"), "te_kimura_bins")
    kim$kimura_divergence_percent <- as_num(kim$kimura_divergence_percent)
    kim$te_bp <- as_num(kim$te_bp)
    kim$assembly_span_bp <- as_num(kim$assembly_span_bp)
    kim$assembly_fraction <- kim$te_bp / kim$assembly_span_bp
    kim$kimura_bin <- cut(kim$kimura_divergence_percent, breaks = seq(0, 60, by = 2),
                          include.lowest = TRUE, right = FALSE)
    summary <- aggregate(assembly_fraction ~ species_id + kimura_bin, kim, sum, na.rm = TRUE)
    write_table(summary, file.path(outdir, "source_fig4d_te_kimura_relative_age_proxy.csv"))
  }
  invisible(TRUE)
}
