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

read_species_tree <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Species tree not found: ", path, call. = FALSE)
    return(NULL)
  }
  quiet_library("ape", required = TRUE)
  ape::read.tree(path)
}

drop_unmatched_tree_tips <- function(tree, data, species_col = "species_id") {
  if (is.null(tree)) return(list(tree = NULL, data = data))
  keep <- intersect(tree$tip.label, as.character(data[[species_col]]))
  if (length(keep) < 3) {
    stop("Fewer than three species overlap between tree and data.", call. = FALSE)
  }
  drop <- setdiff(tree$tip.label, keep)
  tree2 <- if (length(drop)) ape::drop.tip(tree, drop) else tree
  data2 <- data[match(tree2$tip.label, as.character(data[[species_col]])), , drop = FALSE]
  list(tree = tree2, data = data2)
}

fit_pgls_habitat <- function(data, response_col, tree, species_col = "species_id") {
  quiet_library("ape", required = TRUE)
  quiet_library("nlme", required = TRUE)
  require_cols(data, c(species_col, "habitat_binary", response_col), "pgls_data")
  x <- data[!is.na(data[[response_col]]) &
              data$habitat_binary %in% c("seawater_group", "freshwater"), , drop = FALSE]
  if (nrow(x) < 3 || length(unique(x$habitat_binary)) != 2) {
    return(data.frame(
      estimate_freshwater_minus_seawater_group = NA_real_,
      standard_error = NA_real_,
      statistic = NA_real_,
      p = NA_real_,
      model = "PGLS_not_estimated_insufficient_data",
      lambda = NA_real_
    ))
  }
  x$habitat_factor <- factor(x$habitat_binary, levels = c("seawater_group", "freshwater"))
  matched <- tryCatch(
    drop_unmatched_tree_tips(tree, x, species_col = species_col),
    error = function(e) NULL
  )
  if (is.null(matched)) {
    return(data.frame(
      estimate_freshwater_minus_seawater_group = NA_real_,
      standard_error = NA_real_,
      statistic = NA_real_,
      p = NA_real_,
      model = "PGLS_not_estimated_tree_data_mismatch",
      lambda = NA_real_
    ))
  }
  x <- matched$data
  tree <- matched$tree
  if (length(unique(x$habitat_factor)) != 2) {
    return(data.frame(
      estimate_freshwater_minus_seawater_group = NA_real_,
      standard_error = NA_real_,
      statistic = NA_real_,
      p = NA_real_,
      model = "PGLS_not_estimated_single_habitat_after_tree_matching",
      lambda = NA_real_
    ))
  }
  corr <- ape::corPagel(value = 1, phy = tree, form = stats::as.formula(paste0("~", species_col)), fixed = FALSE)
  form <- stats::as.formula(paste(response_col, "~ habitat_factor"))
  fit <- tryCatch(
    nlme::gls(form, data = x, correlation = corr, method = "ML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(
      estimate_freshwater_minus_seawater_group = NA_real_,
      standard_error = NA_real_,
      statistic = NA_real_,
      p = NA_real_,
      model = "PGLS_not_estimated_model_fit_failed",
      lambda = NA_real_
    ))
  }
  tab <- summary(fit)$tTable
  coef_name <- rownames(tab)[grepl("freshwater", rownames(tab), ignore.case = TRUE)][1]
  data.frame(
    estimate_freshwater_minus_seawater_group = unname(tab[coef_name, "Value"]),
    standard_error = unname(tab[coef_name, "Std.Error"]),
    statistic = unname(tab[coef_name, "t-value"]),
    p = unname(tab[coef_name, "p-value"]),
    model = "PGLS_Pagel_lambda_ML",
    lambda = unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  )
}
