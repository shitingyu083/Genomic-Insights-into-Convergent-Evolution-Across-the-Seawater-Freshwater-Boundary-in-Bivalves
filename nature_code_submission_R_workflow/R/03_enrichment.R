source(file.path("R", "00_utils.R"))

hypergeom_enrichment <- function(query_genes, background_genes, gene_pathways) {
  require_cols(gene_pathways, c("gene", "pathway_id", "pathway_name"), "gene_pathway_annotations")
  query_genes <- unique(query_genes[!is.na(query_genes) & nzchar(query_genes)])
  background_genes <- unique(background_genes[!is.na(background_genes) & nzchar(background_genes)])
  gp <- unique(gene_pathways[gene_pathways$gene %in% background_genes, ])
  if (!length(query_genes) || !nrow(gp)) {
    return(data.frame(
      pathway_id = character(), pathway_name = character(), functional_category = character(),
      overlap = integer(), pathway_size = integer(), query_size = integer(),
      background_size = integer(), p = numeric(), q = numeric()
    ))
  }
  pathways <- split(gp$gene, gp$pathway_id)
  out <- lapply(names(pathways), function(pid) {
    pathway_genes <- unique(pathways[[pid]])
    k <- sum(query_genes %in% pathway_genes)
    m <- length(pathway_genes)
    n <- length(setdiff(background_genes, pathway_genes))
    qn <- length(query_genes)
    p <- stats::phyper(k - 1, m, n, qn, lower.tail = FALSE)
    name <- gp$pathway_name[match(pid, gp$pathway_id)]
    category <- if ("functional_category" %in% names(gp)) collapse_semicolon(gp$functional_category[gp$pathway_id == pid]) else ""
    data.frame(pathway_id = pid, pathway_name = name, functional_category = category,
               overlap = k, pathway_size = m, query_size = qn, background_size = length(background_genes), p = p)
  })
  res <- do.call(rbind, out)
  res$q <- bh(res$p)
  res[order(res$q, res$p), ]
}

run_enrichment <- function(paths, outdir, q_threshold = 0.05) {
  annotations <- read_table_auto(paths$gene_pathway_annotations)
  background <- unique(annotations$gene)
  parallel <- read_table_auto(file.path(outdir, "table_candidate_parallel_substitution_signals.csv"), required = FALSE)
  branch <- read_table_auto(file.path(outdir, "table_branch_site_positive_selection_signals.csv"), required = FALSE)

  if (!is.null(parallel)) {
    qgenes <- unique(parallel$gene[parallel$candidate_parallel_signal %in% TRUE])
    e <- hypergeom_enrichment(qgenes, background, annotations)
    e$significant <- e$q < q_threshold
    write_table(e, file.path(outdir, "table_enrichment_candidate_parallel_substitution_genes.csv"))
  }
  if (!is.null(branch)) {
    qgenes <- unique(branch$gene[branch$branch_site_signal %in% TRUE])
    e <- hypergeom_enrichment(qgenes, background, annotations)
    e$significant <- e$q < q_threshold
    write_table(e, file.path(outdir, "table_enrichment_branch_site_positive_selection_signal_genes.csv"))
  }
  invisible(TRUE)
}
