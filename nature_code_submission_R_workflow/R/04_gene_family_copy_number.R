source(file.path("R", "00_utils.R"))

run_gene_family_copy_number <- function(paths, outdir, config) {
  metadata <- standardize_species_metadata(read_table_auto(paths$species_metadata))
  copy <- read_table_auto(paths$gene_family_copy_number)
  require_cols(copy, c("gene_family", "species_id", "copy_number"), "gene_family_copy_number")
  copy$copy_number <- as_num(copy$copy_number)
  x <- merge(copy, metadata[, c("species_id", "habitat_binary")], by = "species_id", all.x = TRUE)
  x <- subset(x, habitat_binary %in% c("freshwater", "seawater_group"))

  families <- unique(x$gene_family)
  res <- lapply(families, function(fam) {
    y <- x[x$gene_family == fam, ]
    occupancy <- mean(!is.na(y$copy_number))
    fw_mean <- mean(y$copy_number[y$habitat_binary == "freshwater"], na.rm = TRUE)
    sw_mean <- mean(y$copy_number[y$habitat_binary == "seawater_group"], na.rm = TRUE)
    wt <- safe_wilcox(y$copy_number, y$habitat_binary)
    data.frame(
      gene_family = fam,
      n_species = sum(!is.na(y$copy_number)),
      occupancy = occupancy,
      freshwater_mean = fw_mean,
      seawater_group_mean = sw_mean,
      delta_freshwater_minus_seawater_group = fw_mean - sw_mean,
      wilcoxon_statistic = wt[["statistic"]],
      p_wilcoxon = wt[["p"]]
    )
  })
  res <- do.call(rbind, res)
  res$q_wilcoxon <- bh(res$p_wilcoxon)
  res$candidate_copy_number_difference <- res$occupancy >= as_num(config[["gene_family_occupancy_min"]]) &
    res$q_wilcoxon < as_num(config[["gene_family_q_threshold"]]) &
    res$delta_freshwater_minus_seawater_group > 0
  res <- res[order(res$q_wilcoxon, -res$delta_freshwater_minus_seawater_group), ]
  write_table(res, file.path(outdir, "table_gene_family_copy_number_associations.csv"))

  x_scaled <- x
  x_scaled$copy_number_scaled_by_family <- ave(x_scaled$copy_number, x_scaled$gene_family, FUN = scale_within_group)
  write_table(x_scaled, file.path(outdir, "source_fig3a_gene_family_heatmap.csv"))
  invisible(res)
}

