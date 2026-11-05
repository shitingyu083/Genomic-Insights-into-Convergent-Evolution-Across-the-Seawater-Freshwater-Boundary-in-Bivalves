# R1 revision, R2 comment 12: extension of aa_gc_check.R (4) from Gly/Ala/Pro/Gln to all 20 residues, because the
# 960-fit table shows sign-consistent FW estimates for Asp (48/48 negative), Ala (48/48 positive), Leu (45/48) and
# Tyr (47/48 negative). Same inputs, tree, subsets and fitting rule as aa_gc_check.R; descriptive only (n = 9-12).
# Output: aa_gc_all20.tsv (set, aa, Spearman rho/P of GC4 vs CLR, FW-only and FW + GC4 PGLS).
suppressMessages({library(ape); library(nlme)})
W <- file.path(Sys.getenv("IZ_ROOT", "."), "_R1_work")
TREE <- file.path(W, "server_data/x/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/SpeciesTree_rooted_node_labels.nwk")
OUT <- file.path(W, "analysis/aa")
tr0 <- read.tree(TREE)
tr0$tip.label <- sub("[.]pep$", "", tr0$tip.label)
tr0$node.label <- NULL
tr <- drop.tip(tr0, "Lgi")
g <- read.delim(file.path(W, "server_data/r1rev/aa_gc/aa_gc_species.tsv"), stringsAsFactors = FALSE)
a0 <- read.delim(file.path(OUT, "aa_species_clr.tsv"), stringsAsFactors = FALSE)
a0 <- a0[a0$dataset == "ALL_ANY", ]
d <- merge(g, a0, by = "species")
stopifnot(nrow(d) == 12, all(d$species %in% tr$tip.label))
d$id <- d$species
d$FW <- as.integer(d$salinity3 == "FRESHWATER")
AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

fit1 <- function(dd, fml, tree) {  # identical to aa_gc_check.R / aa_pgls.R
  vars <- all.vars(fml)
  dd <- dd[complete.cases(dd[, vars]), ]
  tt <- keep.tip(tree, dd$id)
  dd <- dd[match(tt$tip.label, dd$id), ]
  ll <- function(lam) {
    m <- try(gls(fml, data = dd, correlation = corPagel(lam, tt, form = ~id, fixed = TRUE), method = "ML"), silent = TRUE)
    if (inherits(m, "try-error")) -Inf else as.numeric(logLik(m))
  }
  o <- optimize(ll, c(0, 1), maximum = TRUE, tol = 1e-5)
  lam <- if (ll(0) > o$objective) 0 else if (ll(1) > o$objective) 1 else o$maximum
  m <- try(gls(fml, data = dd, correlation = corPagel(lam, tt, form = ~id, fixed = TRUE), method = "REML"), silent = TRUE)
  if (inherits(m, "try-error")) return(c(lambda = NA, beta = NA, lo = NA, hi = NA, P = NA))
  tab <- summary(m)$tTable
  ci <- tab["FW", "Value"] + c(-1, 1) * qt(0.975, nrow(dd) - nrow(tab)) * tab["FW", "Std.Error"]
  c(lambda = lam, beta = tab["FW", "Value"], lo = ci[1], hi = ci[2], P = tab["FW", "p-value"])
}

SETS <- list(S10 = d$species[d$matches_proteome_count == 1], S9 = setdiff(d$species[d$matches_proteome_count == 1], "Lfo"),
             S12 = d$species)
res <- list()
for (s in names(SETS)) {
  dd <- d[d$species %in% SETS[[s]], ]
  for (a in AA) {
    dd$y <- dd[[paste0("clr_", a)]]
    k <- cor.test(dd$GC4, dd$y, method = "spearman", exact = TRUE)
    f0 <- fit1(dd, y ~ FW, tr); f1 <- fit1(dd, y ~ FW + GC4, tr)
    res[[length(res) + 1]] <- data.frame(set = s, aa = a, n = nrow(dd), n_FW = sum(dd$FW),
      rho_GC4 = unname(k$estimate), P_rho = k$p.value,
      beta_FW = f0["beta"], lo_FW = f0[3], hi_FW = f0[4], P_FW = f0["P"],
      beta_FW_adjGC4 = f1["beta"], lo_adj = f1[3], hi_adj = f1[4], P_adj = f1["P"], row.names = NULL)
  }
}
res <- do.call(rbind, res)
write.table(res, file.path(OUT, "aa_gc_all20.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(res[res$set != "S9", c("set", "aa", "rho_GC4", "P_rho", "beta_FW", "P_FW", "beta_FW_adjGC4", "lo_adj", "hi_adj", "P_adj")],
      digits = 3, row.names = FALSE)
