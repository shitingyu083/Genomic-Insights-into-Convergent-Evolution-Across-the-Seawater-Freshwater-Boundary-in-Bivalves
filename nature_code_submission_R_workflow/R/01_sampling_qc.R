source(file.path("R", "00_utils.R"))

run_sampling_qc <- function(paths, outdir) {
  metadata <- standardize_species_metadata(read_table_auto(paths$species_metadata))
  qc <- read_table_auto(paths$genome_qc, required = FALSE)

  retained <- subset(metadata, include_42_bivalves & !is_outgroup)
  habitat_counts <- as.data.frame(table(retained$habitat_binary), stringsAsFactors = FALSE)
  names(habitat_counts) <- c("habitat_binary", "n_species")

  unit_counts <- aggregate(
    species_id ~ freshwater_unit,
    subset(retained, habitat_binary == "freshwater" & !is.na(freshwater_unit) & freshwater_unit != ""),
    function(z) length(unique(z))
  )
  names(unit_counts)[2] <- "n_species"

  summary <- data.frame(
    metric = c(
      "retained_bivalve_species",
      "freshwater_taxa",
      "seawater_group_taxa",
      "phylogenetically_collapsed_freshwater_units"
    ),
    value = c(
      nrow(retained),
      sum(retained$habitat_binary == "freshwater"),
      sum(retained$habitat_binary == "seawater_group"),
      length(unique(unit_counts$freshwater_unit))
    )
  )

  write_table(summary, file.path(outdir, "table_sampling_summary.csv"))
  write_table(habitat_counts, file.path(outdir, "table_sampling_habitat_counts.csv"))
  write_table(unit_counts, file.path(outdir, "table_freshwater_collapsed_units.csv"))

  if (!is.null(qc)) {
    require_cols(qc, c("species_id"), "genome_qc")
    qc_join <- merge(retained, qc, by = "species_id", all.x = TRUE)
    write_table(qc_join, file.path(outdir, "source_fig1_genome_qc.csv"))
  }

  invisible(list(metadata = metadata, summary = summary))
}

