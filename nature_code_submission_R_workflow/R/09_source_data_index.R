source(file.path("R", "00_utils.R"))

run_source_data_index <- function(outdir) {
  files <- list.files(outdir, pattern = "\\.(csv|tsv)$", recursive = FALSE, full.names = FALSE)
  index <- data.frame(
    file = files,
    role = ifelse(grepl("^source_fig", files), "figure_source_data", "analysis_table"),
    stringsAsFactors = FALSE
  )
  write_table(index, file.path(outdir, "source_data_index.csv"))
  invisible(index)
}

