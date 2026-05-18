source(file.path("R", "00_utils.R"))

infer_unit_state <- function(states, metadata) {
  x <- merge(states, metadata[, c("species_id", "habitat_binary", "freshwater_unit")],
             by = "species_id", all.x = TRUE)
  fw <- subset(x, habitat_binary == "freshwater" & !is.na(freshwater_unit) & freshwater_unit != "")
  if (!nrow(fw)) return(data.frame())
  split_fw <- split(fw, fw$freshwater_unit)
  out <- lapply(names(split_fw), function(unit) {
    y <- split_fw[[unit]]
    aa <- unique(y$amino_acid[!is.na(y$amino_acid) & y$amino_acid != "" & y$amino_acid != "-"])
    state <- if (length(aa) == 1) aa else NA_character_
    data.frame(freshwater_unit = unit, unit_state = state, informative = !is.na(state))
  })
  do.call(rbind, out)
}

screen_fsite <- function(extant, metadata, min_units = 2) {
  require_cols(extant, c("gene", "position", "species_id", "amino_acid"), "extant_site_states")
  metadata <- standardize_species_metadata(metadata)
  extant <- merge(extant, metadata[, c("species_id", "habitat_binary", "freshwater_unit")],
                  by = "species_id", all.x = TRUE)
  extant$key <- paste(extant$gene, extant$position, sep = "::")
  keys <- unique(extant$key)

  res <- lapply(keys, function(k) {
    z <- extant[extant$key == k, ]
    units <- infer_unit_state(z[, c("gene", "position", "species_id", "amino_acid")], metadata)
    informative <- subset(units, informative)
    if (!nrow(informative)) {
      return(data.frame(
        gene = z$gene[1], position = z$position[1], candidate_state = NA,
        n_informative_units = 0, n_recurrent_units = 0,
        recurrent_units = "", fsite_pass = FALSE
      ))
    }
    tab <- sort(table(informative$unit_state), decreasing = TRUE)
    candidate_state <- names(tab)[1]
    n_rec <- as.integer(tab[1])
    rec_units <- informative$freshwater_unit[informative$unit_state == candidate_state]
    data.frame(
      gene = z$gene[1],
      position = z$position[1],
      candidate_state = candidate_state,
      n_informative_units = nrow(informative),
      n_recurrent_units = n_rec,
      recurrent_units = paste(rec_units, collapse = ";"),
      fsite_pass = n_rec >= min_units
    )
  })
  do.call(rbind, res)
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

summarize_parallel <- function(fsite, extant_parallel = NULL, min_units = 2) {
  if (is.null(extant_parallel)) {
    out <- fsite
    out$candidate_parallel_signal <- out$fsite_pass & out$n_recurrent_units >= min_units
    return(out[, c("gene", "position", "candidate_state", "n_informative_units",
                   "n_recurrent_units", "recurrent_units", "candidate_parallel_signal")])
  }
  require_cols(extant_parallel, c("gene", "position"), "parallel_substitution_results")
  fsite$key <- paste(fsite$gene, fsite$position, sep = "::")
  extant_parallel$key <- paste(extant_parallel$gene, extant_parallel$position, sep = "::")
  out <- merge(fsite, extant_parallel, by = "key", all.x = TRUE, suffixes = c("", "_external"))
  out$candidate_parallel_signal <- out$fsite_pass & !is.na(out$gene_external)
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
  fsite <- screen_fsite(extant, metadata, min_units = as.integer(config[["freshwater_recurrence_min_units"]]))
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
  parallel_sum <- summarize_parallel(fsite, parallel_external,
                                     min_units = as.integer(config[["freshwater_recurrence_min_units"]]))
  write_table(parallel_sum, file.path(outdir, "table_candidate_parallel_substitution_signals.csv"))

  overlap <- merge(
    parallel_sum[, c("gene", "position", "candidate_parallel_signal")],
    if (is.null(branch_sum)) data.frame(gene = character(), position = integer(), branch_site_signal = logical())
    else branch_sum[, c("gene", "position", "branch_site_signal")],
    by = c("gene", "position"), all = TRUE
  )
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

