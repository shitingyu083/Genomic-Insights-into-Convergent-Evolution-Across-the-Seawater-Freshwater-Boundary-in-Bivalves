source(file.path("R", "00_utils.R"))

clr_transform <- function(counts) {
  x <- as_num(counts)
  if (any(is.na(x))) return(rep(NA_real_, length(x)))
  if (any(x == 0)) {
    min_nonzero <- min(x[x > 0], na.rm = TRUE)
    if (!is.finite(min_nonzero)) min_nonzero <- 1
    x <- x + 0.5 * min_nonzero
  }
  freq <- x / sum(x)
  log(freq) - mean(log(freq))
}

run_amino_acid_clr <- function(paths, outdir) {
  metadata <- standardize_species_metadata(read_table_auto(paths$species_metadata))
  aa <- read_table_auto(paths$amino_acid_counts)
  require_cols(aa, c("species_id", "protein_set", "amino_acid", "count"), "amino_acid_counts")
  aa$count <- as_num(aa$count)
  aa <- merge(aa, metadata[, c("species_id", "habitat_binary")], by = "species_id", all.x = TRUE)
  aa <- subset(aa, habitat_binary %in% c("freshwater", "seawater_group"))

  split_key <- paste(aa$species_id, aa$protein_set, sep = "::")
  aa$clr <- NA_real_
  for (key in unique(split_key)) {
    idx <- which(split_key == key)
    aa$clr[idx] <- clr_transform(aa$count[idx])
  }

  tests <- do.call(rbind, lapply(split(aa, list(aa$protein_set, aa$amino_acid), drop = TRUE), function(y) {
    wt <- safe_wilcox(y$clr, y$habitat_binary)
    data.frame(
      protein_set = y$protein_set[1],
      amino_acid = y$amino_acid[1],
      freshwater_mean_clr = mean(y$clr[y$habitat_binary == "freshwater"], na.rm = TRUE),
      seawater_group_mean_clr = mean(y$clr[y$habitat_binary == "seawater_group"], na.rm = TRUE),
      p = wt[["p"]],
      statistic = wt[["statistic"]]
    )
  }))
  tests$q <- ave(tests$p, tests$protein_set, FUN = bh)
  tests$significant_q_0_05 <- tests$q < 0.05
  write_table(aa, file.path(outdir, "source_fig5_amino_acid_clr_long.csv"))
  write_table(tests, file.path(outdir, "table_amino_acid_clr_associations.csv"))
  invisible(tests)
}

