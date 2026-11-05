# R2-M3(ii): which simple model reproduces the submitted Fig. 4a assembly-span P = 0.014 (q = 0.031)?
# Read-only check on the authors' archived trait table (step3_pgls/input/genomic_traits.csv) and Table S1 spans.
suppressMessages({library(ape); library(nlme)})
W <- file.path(Sys.getenv("IZ_ROOT", "."), "_R1_work")
IN <- file.path(W, "server_data/x/orthofinder1229/workflow/convergence_analysis_workflow/step3_pgls/input")
g <- read.csv(file.path(IN, "genomic_traits.csv"), stringsAsFactors = FALSE)
g$id <- sub("\\.pep$", "", g$species)
g <- g[g$habitat != "outgroup", ]
g$FW <- as.integer(g$habitat == "freshwater")
s1 <- read.delim(file.path(W, "analysis/pgls_span/traits_S1_parsed.tsv"), stringsAsFactors = FALSE)
m <- merge(g[, c("id", "FW", "genome_size")], s1[, c("id", "span_Mb")], by = "id", all = TRUE)
cat("authors genome_size vs Table S1 span_Mb (n both):", sum(complete.cases(m)),
    " max |diff| Mb:", max(abs(m$genome_size - m$span_Mb), na.rm = TRUE), "\n")
print(m[!complete.cases(m) | abs(m$genome_size - m$span_Mb) > 1, ], row.names = FALSE)
out <- list()
add <- function(name, p, n) out[[length(out) + 1]] <<- data.frame(test = name, n = n, P = p)
for (src in c("authors_genomic_traits", "TableS1_span")) {
  d <- if (src == "authors_genomic_traits") data.frame(id = g$id, FW = g$FW, y = g$genome_size) else
    data.frame(id = s1$id, FW = as.integer(s1$salinity3 == "FRESHWATER"), y = s1$span_Mb)
  d <- d[!is.na(d$y), ]
  n <- nrow(d)
  add(paste(src, "Wilcoxon raw"), wilcox.test(d$y[d$FW == 1], d$y[d$FW == 0])$p.value, n)
  add(paste(src, "Welch t raw"), t.test(d$y[d$FW == 1], d$y[d$FW == 0])$p.value, n)
  add(paste(src, "Welch t log10"), t.test(log10(d$y[d$FW == 1]), log10(d$y[d$FW == 0]))$p.value, n)
  add(paste(src, "Student t raw"), t.test(d$y[d$FW == 1], d$y[d$FW == 0], var.equal = TRUE)$p.value, n)
  add(paste(src, "Student t log10"), t.test(log10(d$y[d$FW == 1]), log10(d$y[d$FW == 0]), var.equal = TRUE)$p.value, n)
  for (tf in c("authors_species_tree", "reference_tree")) {
    tr <- if (tf == "authors_species_tree") read.tree(file.path(IN, "species_tree.nwk")) else
      read.tree(file.path(W, "analysis/iqtree_rerun/out/fixed_topology/iqtree_species_tree_codes.nwk"))
    tr$tip.label <- sub("\\.pep$", "", tr$tip.label); tr$node.label <- NULL
    tr <- keep.tip(tr, intersect(tr$tip.label, d$id))
    dd <- d[match(tr$tip.label, d$id), ]; dd <- dd[!is.na(dd$id), ]
    if (is.null(tr$edge.length)) next
    for (y in c("raw", "log10")) {
      dd$yy <- if (y == "raw") dd$y else log10(dd$y)
      for (lam in c(0, 0.442, 1)) {
        f <- try(gls(yy ~ FW, data = dd, correlation = corPagel(lam, tr, form = ~id, fixed = TRUE), method = "REML"), silent = TRUE)
        if (!inherits(f, "try-error")) add(sprintf("%s %s PGLS %s lambda=%.3f", src, tf, y, lam), summary(f)$tTable["FW", "p-value"], nrow(dd))
      }
    }
  }
}
res <- do.call(rbind, out)
res$near_0.014 <- abs(res$P - 0.014) < 0.002
print(res, digits = 4, row.names = FALSE)
write.table(res, file.path(W, "analysis/pgls_span/fig4a_p_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
