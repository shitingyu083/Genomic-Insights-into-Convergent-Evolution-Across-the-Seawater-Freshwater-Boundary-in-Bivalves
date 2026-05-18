options(stringsAsFactors = FALSE)

quiet_library <- function(pkg, required = TRUE) {
  ok <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
  if (!ok && required) {
    stop(
      "Required R package not installed: ", pkg,
      ". Install it before running this module.",
      call. = FALSE
    )
  }
  ok
}

read_table_auto <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Input file not found: ", path, call. = FALSE)
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    utils::read.csv(path, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, check.names = FALSE)
  } else if (ext %in% c("rds")) {
    readRDS(path)
  } else {
    stop("Unsupported table extension for ", path, call. = FALSE)
  }
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

require_cols <- function(x, cols, name = deparse(substitute(x))) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

standardize_species_metadata <- function(metadata) {
  require_cols(
    metadata,
    c("species_id", "species_name", "habitat", "habitat_binary",
      "freshwater_unit", "is_outgroup", "include_42_bivalves", "include_site_level"),
    "species_metadata"
  )
  parse_bool <- function(x) {
    if (is.logical(x)) return(x)
    y <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(y))
    out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[y %in% c("false", "f", "0", "no", "n", "")] <- FALSE
    out
  }
  metadata$is_outgroup <- parse_bool(metadata$is_outgroup)
  metadata$include_42_bivalves <- parse_bool(metadata$include_42_bivalves)
  metadata$include_site_level <- parse_bool(metadata$include_site_level)
  metadata
}

read_config <- function(path) {
  cfg <- read_table_auto(path)
  require_cols(cfg, c("key", "value"), "analysis_config")
  values <- stats::setNames(cfg$value, cfg$key)
  values
}

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

bh <- function(p) {
  stats::p.adjust(p, method = "BH")
}

safe_lm_p <- function(formula, data) {
  fit <- stats::lm(formula, data = data)
  sm <- summary(fit)
  coef_table <- sm$coefficients
  list(fit = fit, coef_table = coef_table)
}

safe_wilcox <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]
  g <- as.factor(g[ok])
  if (length(levels(g)) != 2 || any(table(g) < 1)) {
    return(c(statistic = NA_real_, p = NA_real_))
  }
  wt <- suppressWarnings(stats::wilcox.test(x ~ g, exact = FALSE))
  c(statistic = unname(wt$statistic), p = wt$p.value)
}

scale_within_group <- function(x) {
  if (all(is.na(x))) return(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

collapse_semicolon <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  paste(x, collapse = ";")
}

message_step <- function(...) {
  message("[workflow] ", paste0(..., collapse = ""))
}
