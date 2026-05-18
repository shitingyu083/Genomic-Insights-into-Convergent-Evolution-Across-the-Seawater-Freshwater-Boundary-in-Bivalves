source(file.path("R", "00_utils.R"))

run_structure_foldx <- function(paths, outdir) {
  foldx <- read_table_auto(paths$atg7_foldx_results, required = FALSE)
  if (is.null(foldx)) return(invisible(NULL))
  require_cols(foldx, c("gene", "site", "substitution", "run_id", "ddg_kcal_mol"), "atg7_foldx_results")
  foldx$ddg_kcal_mol <- as_num(foldx$ddg_kcal_mol)
  summary <- aggregate(ddg_kcal_mol ~ gene + site + substitution, foldx, function(z) {
    c(n = sum(!is.na(z)), mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
  })
  summary <- do.call(data.frame, summary)
  names(summary) <- sub("ddg_kcal_mol\\.", "", names(summary))
  write_table(foldx, file.path(outdir, "source_fig6e_foldx_runs.csv"))
  write_table(summary, file.path(outdir, "table_atg7_foldx_ddg_summary.csv"))
  invisible(summary)
}

