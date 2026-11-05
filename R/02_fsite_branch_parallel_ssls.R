source(file.path("R", "00_utils.R"))

AMINO_ACIDS <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

empty_fsite_table <- function() {
  data.frame(
    gene = character(), position = character(), fw_mrca_state = character(), derived_state = character(),
    n_target_taxa = integer(), target_taxa = character(), n_units_sharing = integer(),
    units_sharing = character(), filter1_unit_stem_origin = logical(),
    filter2_absent_nonfreshwater = logical(), descriptor_two_or_more_units = logical(),
    fsite_pass = logical()
  )
}

# R1 Fsite screen (Methods, Fsite screening). One record is one derived freshwater residue at one column.
# Candidate set: the residue differs from the reconstructed FW-MRCA state and is shared by at least
# `min_target_taxa` freshwater target taxa; no lineage-number, probability or FDR threshold is applied.
# Filter 1: the derived residue arose on a freshwater-unit stem (from `origin`, written by
#   R1_reanalysis/_R1_work/analysis/r2m1_fsite/scan_reproduce.py from the codeml node states).
# Filter 2: the derived residue is absent from all non-freshwater taxa, the outgroup included.
# Sharing by two or more freshwater units is recorded as a descriptor, not used as a filter.
screen_fsite <- function(extant, metadata, ancestral, origin = NULL, target_taxa, min_target_taxa = 3) {
  require_cols(extant, c("gene", "position", "species_id", "amino_acid"), "extant_site_states")
  require_cols(ancestral, c("gene", "position", "node", "ancestral_state"), "ancestral_states")
  metadata <- standardize_species_metadata(metadata)
  extant$position <- as.character(extant$position)
  extant <- merge(extant, metadata[, c("species_id", "habitat_binary", "freshwater_unit")],
                  by = "species_id", all.x = TRUE)
  extant$key <- paste(extant$gene, extant$position, sep = "::")
  mrca <- ancestral[ancestral$node == "FW_MRCA", , drop = FALSE]
  mrca_state <- stats::setNames(as.character(mrca$ancestral_state),
                                paste(mrca$gene, as.character(mrca$position), sep = "::"))
  is_residue <- function(a) !is.na(a) & a %in% AMINO_ACIDS

  res <- lapply(split(extant, extant$key), function(z) {
    anc <- unname(mrca_state[z$key[1]])
    if (!length(anc) || is.na(anc)) return(NULL)
    tg <- z[z$species_id %in% target_taxa & is_residue(z$amino_acid), , drop = FALSE]
    derived <- tg[tg$amino_acid != anc, , drop = FALSE]
    if (!nrow(derived)) return(NULL)
    tab <- table(derived$amino_acid)
    keep <- names(tab)[tab >= min_target_taxa]
    if (!length(keep)) return(NULL)
    nonfw <- z$amino_acid[!is.na(z$habitat_binary) & z$habitat_binary != "freshwater" &
                            is_residue(z$amino_acid)]
    do.call(rbind, lapply(keep, function(s) {
      sp <- derived$species_id[derived$amino_acid == s]
      units <- unique(derived$freshwater_unit[derived$amino_acid == s])
      units <- units[!is.na(units) & nzchar(units)]
      data.frame(
        gene = z$gene[1], position = z$position[1], fw_mrca_state = anc, derived_state = s,
        n_target_taxa = length(sp), target_taxa = paste(sort(sp), collapse = ";"),
        n_units_sharing = length(units), units_sharing = paste(sort(units), collapse = ";"),
        filter2_absent_nonfreshwater = !(s %in% nonfw)
      )
    }))
  })
  out <- do.call(rbind, res)
  if (is.null(out) || !nrow(out)) return(empty_fsite_table())

  if (!is.null(origin)) {
    require_cols(origin, c("gene", "position", "derived_state", "origin_on_unit_stem"), "fsite_origin_branches")
    origin$position <- as.character(origin$position)
    origin$filter1_unit_stem_origin <- as.logical(origin$origin_on_unit_stem)
    out <- merge(out, origin[, c("gene", "position", "derived_state", "filter1_unit_stem_origin")],
                 by = c("gene", "position", "derived_state"), all.x = TRUE)
  } else {
    out$filter1_unit_stem_origin <- NA
  }
  out$descriptor_two_or_more_units <- out$n_units_sharing >= 2
  out$fsite_pass <- out$filter1_unit_stem_origin & out$filter2_absent_nonfreshwater
  rownames(out) <- NULL
  out[, names(empty_fsite_table())]
}

summarize_branch_site <- function(branch_site, q_threshold = 0.05, beb_threshold = 0.95) {
  require_cols(branch_site, c("gene", "position", "lnL_null", "lnL_alt", "beb_posterior"), "branch_site_results")
  branch_site$lnL_null <- as_num(branch_site$lnL_null)
  branch_site$lnL_alt <- as_num(branch_site$lnL_alt)
  branch_site$beb_posterior <- as_num(branch_site$beb_posterior)
  branch_site$lrt <- pmax(0, 2 * (branch_site$lnL_alt - branch_site$lnL_null))
  branch_site$p_chisq1_conservative <- stats::pchisq(branch_site$lrt, df = 1, lower.tail = FALSE)
  branch_site$p_mixture_50_50 <- 0.5 * stats::pchisq(branch_site$lrt, df = 1, lower.tail = FALSE)
  branch_site$q_chisq1_conservative <- bh(branch_site$p_chisq1_conservative)
  branch_site$q_mixture_50_50 <- bh(branch_site$p_mixture_50_50)
  branch_site$branch_site_signal <- branch_site$q_chisq1_conservative < q_threshold &
    branch_site$beb_posterior >= beb_threshold
  branch_site
}

# Candidate parallel signal = Fsite record that passes both filters; with an external parallel-substitution
# table, the record must also appear there.
summarize_parallel <- function(fsite, extant_parallel = NULL) {
  out <- fsite
  out$candidate_parallel_signal <- out$fsite_pass %in% TRUE
  if (!is.null(extant_parallel) && nrow(out)) {
    require_cols(extant_parallel, c("gene", "position"), "parallel_substitution_results")
    ext_key <- paste(extant_parallel$gene, as.character(extant_parallel$position), sep = "::")
    out$in_external_parallel_table <- paste(out$gene, out$position, sep = "::") %in% ext_key
    out$candidate_parallel_signal <- out$candidate_parallel_signal & out$in_external_parallel_table
  }
  out
}

summarize_ssls <- function(ssls) {
  require_cols(ssls, c("gene", "position", "lnL_species_tree", "lnL_freshwater_constrained", "site_class"), "ssls_results")
  ssls$lnL_species_tree <- as_num(ssls$lnL_species_tree)
  ssls$lnL_freshwater_constrained <- as_num(ssls$lnL_freshwater_constrained)
  ssls$delta_ssls <- ssls$lnL_species_tree - ssls$lnL_freshwater_constrained
  classes <- unique(ssls$site_class)
  if (length(classes) > 1) {
    kw <- stats::kruskal.test(delta_ssls ~ site_class, data = ssls)
    attr(ssls, "class_test") <- data.frame(test = "Kruskal-Wallis", statistic = unname(kw$statistic), p = kw$p.value)
  }
  ssls
}

run_site_level_workflow <- function(paths, outdir, config) {
  metadata <- standardize_species_metadata(read_table_auto(paths$species_metadata))
  extant <- read_table_auto(paths$extant_site_states)
  ancestral <- read_table_auto(paths$ancestral_states)
  origin <- read_table_auto(paths$fsite_origin_branches, required = FALSE)
  fsite <- screen_fsite(
    extant, metadata, ancestral, origin = origin,
    target_taxa = strsplit(config[["fsite_target_taxa"]], ";")[[1]],
    min_target_taxa = as.integer(config[["fsite_min_target_taxa"]])
  )
  write_table(fsite, file.path(outdir, "table_fsite_screen.csv"))

  branch <- read_table_auto(paths$branch_site_results, required = FALSE)
  if (!is.null(branch)) {
    branch_sum <- summarize_branch_site(
      branch,
      q_threshold = as_num(config[["branch_site_q_threshold"]]),
      beb_threshold = as_num(config[["branch_site_beb_threshold"]])
    )
    write_table(branch_sum, file.path(outdir, "table_branch_site_positive_selection_signals.csv"))
  } else {
    branch_sum <- NULL
  }

  parallel_external <- read_table_auto(paths$parallel_substitution_results, required = FALSE)
  parallel_sum <- summarize_parallel(fsite, parallel_external)
  write_table(parallel_sum, file.path(outdir, "table_candidate_parallel_substitution_signals.csv"))

  # Fig. 2a overlap is counted per alignment position.
  par_pos <- if (nrow(parallel_sum)) {
    aggregate(candidate_parallel_signal ~ gene + position, parallel_sum, any)
  } else {
    data.frame(gene = character(), position = character(), candidate_parallel_signal = logical())
  }
  br_pos <- if (is.null(branch_sum) || !nrow(branch_sum)) {
    data.frame(gene = character(), position = character(), branch_site_signal = logical())
  } else {
    b <- branch_sum[, c("gene", "position", "branch_site_signal")]
    b$position <- as.character(b$position)
    b$branch_site_signal <- b$branch_site_signal %in% TRUE
    aggregate(branch_site_signal ~ gene + position, b, any)
  }
  overlap <- merge(par_pos, br_pos, by = c("gene", "position"), all = TRUE)
  overlap$candidate_parallel_signal[is.na(overlap$candidate_parallel_signal)] <- FALSE
  overlap$branch_site_signal[is.na(overlap$branch_site_signal)] <- FALSE
  overlap$both_filters <- overlap$candidate_parallel_signal & overlap$branch_site_signal
  write_table(overlap, file.path(outdir, "source_fig2a_site_filter_overlap.csv"))

  ssls <- read_table_auto(paths$ssls_results, required = FALSE)
  if (!is.null(ssls)) {
    ssls_sum <- summarize_ssls(ssls)
    write_table(ssls_sum, file.path(outdir, "source_fig2b_ssls.csv"))
    class_test <- attr(ssls_sum, "class_test")
    if (!is.null(class_test)) write_table(class_test, file.path(outdir, "table_ssls_class_test.csv"))
  }

  invisible(list(fsite = fsite, branch = branch_sum, parallel = parallel_sum, overlap = overlap))
}
