source(file.path("R", "00_utils.R"))

run_salinity_assays <- function(paths, outdir) {
  assay <- read_table_auto(paths$salinity_assay_measurements, required = FALSE)
  if (!is.null(assay)) {
    require_cols(assay, c("assay", "individual_id", "treatment_psu", "time_point", "value"), "salinity_assay_measurements")
    assay$value <- as_num(assay$value)
    assay$treatment_psu <- as_num(assay$treatment_psu)
    summary <- aggregate(value ~ assay + treatment_psu + time_point, assay, function(z) {
      c(n = sum(!is.na(z)), mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
    })
    summary <- do.call(data.frame, summary)
    names(summary) <- sub("value\\.", "", names(summary))
    tests <- do.call(rbind, lapply(split(assay, assay$assay), function(y) {
      if (length(unique(y$treatment_psu[!is.na(y$value)])) < 2) {
        return(data.frame(assay = y$assay[1], test = "Kruskal-Wallis", statistic = NA, p = NA))
      }
      kt <- stats::kruskal.test(value ~ factor(treatment_psu), data = y)
      data.frame(assay = y$assay[1], test = "Kruskal-Wallis", statistic = unname(kt$statistic), p = kt$p.value)
    }))
    tests$q <- bh(tests$p)
    pairwise <- do.call(rbind, lapply(split(assay, assay$assay), function(y) {
      tr <- sort(unique(y$treatment_psu[!is.na(y$value)]))
      if (length(tr) < 2) {
        return(data.frame(
          assay = y$assay[1], treatment_1 = NA, treatment_2 = NA,
          test = "Wilcoxon rank-sum", statistic = NA, p = NA
        ))
      }
      pairs <- utils::combn(tr, 2, simplify = FALSE)
      do.call(rbind, lapply(pairs, function(pp) {
        z <- y[y$treatment_psu %in% pp, ]
        wt <- safe_wilcox(z$value, factor(z$treatment_psu))
        data.frame(
          assay = y$assay[1],
          treatment_1 = pp[1],
          treatment_2 = pp[2],
          test = "Wilcoxon rank-sum",
          statistic = wt[["statistic"]],
          p = wt[["p"]]
        )
      }))
    }))
    pairwise$q_within_assay <- ave(pairwise$p, pairwise$assay, FUN = bh)
    write_table(assay, file.path(outdir, "source_fig6a_c_d_salinity_assay_raw.csv"))
    write_table(summary, file.path(outdir, "table_salinity_assay_summary.csv"))
    write_table(tests, file.path(outdir, "table_salinity_assay_tests.csv"))
    write_table(pairwise, file.path(outdir, "table_salinity_assay_pairwise_tests.csv"))
  }

  expr <- read_table_auto(paths$public_transcriptome_expression, required = FALSE)
  de <- read_table_auto(paths$public_transcriptome_deseq2, required = FALSE)
  if (!is.null(expr)) {
    require_cols(expr, c("gene_symbol", "treatment", "replicate", "normalized_expression"), "public_transcriptome_expression")
    expr$normalized_expression <- as_num(expr$normalized_expression)
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
