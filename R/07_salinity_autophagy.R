source(file.path("R", "00_utils.R"))

# Two-sided Welch t test for one pair of treatments (Methods, salinity assays).
welch_pair <- function(y, t1, t2) {
  a <- y$value[y$treatment_psu == t1 & !is.na(y$value)]
  b <- y$value[y$treatment_psu == t2 & !is.na(y$value)]
  base <- data.frame(
    assay = y$assay[1], treatment_1 = t1, treatment_2 = t2, test = "Welch two-sample t test (two-sided)",
    n_1 = length(a), n_2 = length(b),
    mean_1 = if (length(a)) mean(a) else NA, mean_2 = if (length(b)) mean(b) else NA,
    t = NA, df = NA, p = NA
  )
  if (length(a) < 2 || length(b) < 2) return(base)
  tt <- stats::t.test(a, b, var.equal = FALSE)
  base$t <- unname(tt$statistic)
  base$df <- unname(tt$parameter)
  base$p <- tt$p.value
  base
}

# Comparisons and multiplicity adjustment follow `salinity_comparisons.csv`
# (assay, treatment_1, treatment_2, adjust_method in bonferroni / BH / none).
# Without that table, every pair within an assay is tested with Bonferroni adjustment.
run_salinity_assays <- function(paths, outdir) {
  assay <- read_table_auto(paths$salinity_assay_measurements, required = FALSE)
  if (!is.null(assay)) {
    require_cols(assay, c("assay", "individual_id", "treatment_psu", "time_point", "value"), "salinity_assay_measurements")
    assay$value <- as_num(assay$value)
    assay$treatment_psu <- as_num(assay$treatment_psu)
    obs <- assay[!is.na(assay$value), , drop = FALSE]
    write_table(assay, file.path(outdir, "source_fig6a_c_d_salinity_assay_raw.csv"))
    if (nrow(obs)) {
      summary <- aggregate(value ~ assay + treatment_psu + time_point, obs, function(z) {
        c(n = length(z), mean = mean(z), sd = if (length(z) > 1) stats::sd(z) else NA)
      })
      summary <- do.call(data.frame, summary)
      names(summary) <- sub("value\\.", "", names(summary))
      write_table(summary, file.path(outdir, "table_salinity_assay_summary.csv"))
    }

    comps <- read_table_auto(paths$salinity_comparisons, required = FALSE)
    if (is.null(comps)) {
      comps <- do.call(rbind, lapply(split(obs, obs$assay), function(y) {
        tr <- sort(unique(y$treatment_psu), decreasing = TRUE)
        if (length(tr) < 2) return(NULL)
        pp <- utils::combn(tr, 2)
        data.frame(assay = y$assay[1], treatment_1 = pp[1, ], treatment_2 = pp[2, ], adjust_method = "bonferroni")
      }))
    }
    if (!is.null(comps) && nrow(comps)) {
      require_cols(comps, c("assay", "treatment_1", "treatment_2", "adjust_method"), "salinity_comparisons")
      tests <- do.call(rbind, lapply(seq_len(nrow(comps)), function(i) {
        y <- assay[assay$assay == comps$assay[i], , drop = FALSE]
        if (!nrow(y)) y <- data.frame(assay = comps$assay[i], treatment_psu = NA, value = NA)
        r <- welch_pair(y, as_num(comps$treatment_1[i]), as_num(comps$treatment_2[i]))
        r$adjust_method <- comps$adjust_method[i]
        r
      }))
      grp <- paste(tests$assay, tests$adjust_method)
      tests$q <- NA_real_
      for (g in unique(grp)) {
        idx <- which(grp == g)
        m <- tests$adjust_method[idx[1]]
        tests$q[idx] <- if (identical(m, "none")) tests$p[idx] else stats::p.adjust(tests$p[idx], method = m)
      }
      write_table(tests, file.path(outdir, "table_salinity_assay_pairwise_tests.csv"))
    }
  }

  expr <- read_table_auto(paths$public_transcriptome_expression, required = FALSE)
  de <- read_table_auto(paths$public_transcriptome_deseq2, required = FALSE)
  if (!is.null(expr)) {
    require_cols(expr, c("gene_symbol", "treatment", "replicate", "normalized_expression"), "public_transcriptome_expression")
    expr$normalized_expression <- as_num(expr$normalized_expression)
    expr <- expr[!is.na(expr$normalized_expression), , drop = FALSE]
    if (!nrow(expr)) {
      message("[transcriptome] public_transcriptome_expression has no values; Fig. 6f tables skipped")
      expr <- NULL
    }
  }
  if (!is.null(expr)) {
    expr$log2_expression <- log2(expr$normalized_expression + 1)
    treatment_class <- ifelse(expr$treatment %in% c("AL1", "AL5", "CL"), "low_salinity_10psu", "normal_high_salinity_30_45psu")
    expr$treatment_class <- treatment_class
    means <- aggregate(log2_expression ~ gene_symbol + treatment + treatment_class, expr, mean, na.rm = TRUE)
    means$z_by_gene <- ave(means$log2_expression, means$gene_symbol, FUN = scale_within_group)
    class_means <- aggregate(log2_expression ~ gene_symbol + treatment_class, expr, mean, na.rm = TRUE)
    wide <- reshape(class_means, idvar = "gene_symbol", timevar = "treatment_class", direction = "wide")
    low_col <- "log2_expression.low_salinity_10psu"
    high_col <- "log2_expression.normal_high_salinity_30_45psu"
    if (all(c(low_col, high_col) %in% names(wide))) {
      wide$delta_low_minus_normal_high <- wide[[low_col]] - wide[[high_col]]
    }
    if (!is.null(de)) {
      require_cols(de, c("gene_symbol", "contrast", "log2FoldChange", "pvalue", "padj"), "public_transcriptome_deseq2")
      de$log2FoldChange <- as_num(de$log2FoldChange)
      de$pvalue <- as_num(de$pvalue)
      de$padj <- as_num(de$padj)
      write_table(de, file.path(outdir, "table_public_transcriptome_deseq2_autophagy_genes.csv"))
    }
    write_table(means, file.path(outdir, "source_fig6f_public_transcriptome_heatmap.csv"))
    write_table(wide, file.path(outdir, "table_public_transcriptome_low_salinity_delta.csv"))
  }
  invisible(TRUE)
}
