# R1 revision, R2 M4 (comment 12) / R1 comment 8 / Fig. 5 / P76-P77: PGLS of species-level amino-acid composition (CLR,
# aa_species_clr.tsv from aa_prep.py) on the OrthoFinder rooted species tree; same tree and fitting rule as
# analysis/genefam/genefam_pgls.R (gls + Pagel's lambda by bounded ML on [0, 1], REML coefficients at lambda-hat).
# beta_FW = phylogenetically adjusted difference in CLR (freshwater minus comparator group).
# Adjusted models add lenc (ln mean protein length, or the orthogroup-centred ln length for the one-to-one orthologue
# sets) and gcx = the amino-acid GC proxy ln(GARP) - ln(FYMINK) with the response residue left out of its group.
# BH across the 20 amino acids within each dataset x model.
suppressMessages({library(ape); library(nlme)})
W <- file.path(Sys.getenv("IZ_ROOT", "."), "_R1_work")
TREE <- file.path(W, "server_data/x/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/SpeciesTree_rooted_node_labels.nwk")
OUT <- file.path(W, "analysis/aa")
tr0 <- read.tree(TREE)
tr0$tip.label <- sub("[.]pep$", "", tr0$tip.label)
tr0$node.label <- NULL
tr <- drop.tip(tr0, "Lgi")
tra <- read.delim(file.path(W, "analysis/pgls_span/traits_S1_parsed.tsv"), stringsAsFactors = FALSE)
d0 <- read.delim(file.path(OUT, "aa_species_clr.tsv"), stringsAsFactors = FALSE)
d0 <- merge(d0, tra[, c("id", "salinity3", "BUSCO_C")], by.x = "species", by.y = "id", suffixes = c("", ".S1"))
stopifnot(all(d0$salinity3 == d0$salinity3.S1), all(d0$species %in% tr$tip.label))
d0$id <- d0$species
d0$FW <- as.integer(d0$salinity3 == "FRESHWATER")

fit1 <- function(dd, fml, tree) {
  vars <- all.vars(fml)
  dd <- dd[complete.cases(dd[, vars]), ]
  tt <- keep.tip(tree, dd$id)
  dd <- dd[match(tt$tip.label, dd$id), ]
  na <- data.frame(n = nrow(dd), n_FW = sum(dd$FW), lambda = NA, LR_lambda0 = NA, P_lambda0 = NA, beta_FW = NA, se = NA,
                   ci_lo = NA, ci_hi = NA, t = NA, P = NA, other_terms = "")
  ll <- function(lam) {
    m <- try(gls(fml, data = dd, correlation = corPagel(lam, tt, form = ~id, fixed = TRUE), method = "ML"), silent = TRUE)
    if (inherits(m, "try-error")) -Inf else as.numeric(logLik(m))
  }
  o <- optimize(ll, c(0, 1), maximum = TRUE, tol = 1e-5)
  lam <- if (ll(0) > o$objective) 0 else if (ll(1) > o$objective) 1 else o$maximum
  m <- try(gls(fml, data = dd, correlation = corPagel(lam, tt, form = ~id, fixed = TRUE), method = "REML"), silent = TRUE)
  if (inherits(m, "try-error")) return(within(na, other_terms <- "gls failed"))
  tab <- summary(m)$tTable
  ci <- tab["FW", "Value"] + c(-1, 1) * qt(0.975, nrow(dd) - nrow(tab)) * tab["FW", "Std.Error"]
  lr0 <- 2 * (ll(lam) - ll(0))
  ot <- if (nrow(tab) > 2) paste(sprintf("%s=%.4g(P=%.3g)", rownames(tab)[-(1:2)], tab[-(1:2), "Value"],
                                         tab[-(1:2), "p-value"]), collapse = ";") else ""
  data.frame(n = nrow(dd), n_FW = sum(dd$FW), lambda = lam, LR_lambda0 = lr0, P_lambda0 = pchisq(lr0, 1, lower.tail = FALSE),
             beta_FW = tab["FW", "Value"], se = tab["FW", "Std.Error"], ci_lo = ci[1], ci_hi = ci[2],
             t = tab["FW", "t-value"], P = tab["FW", "p-value"], other_terms = ot)
}

AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
DS <- c("ALL_ANY", "OGREP_ANY", "SC8_OGC", "SC36_OGC", "ALL_EXT", "OGREP_EXT", "ALL_MEM", "ALL_CYT")
U <- y ~ FW
A <- y ~ FW + lenc + gcx
res <- list(); vif <- list()
for (ds in DS) for (a in AA) {
  d <- d0[d0$dataset == ds, ]
  d$y <- d[[paste0("clr_", a)]]
  d$gcx <- d[[paste0("gcp_", a)]]
  S <- d$salinity3
  hq <- d$BUSCO_C >= 85
  desc <- data.frame(mean_FW = mean(d$y[d$FW == 1]), mean_comp = mean(d$y[d$FW == 0]),
                     mean_euh = mean(d$y[S == "MARINE"]), mean_brack = mean(d$y[S == "BRACKISH"]))
  mods <- list(
    list("1 main: FW vs brackish+euhaline", d, U),
    list("2 comparator euhaline only", d[S != "BRACKISH", ], U),
    list("3 BUSCO C < 85% excluded", d[hq, ], U),
    list("4 adjusted (ln length, GC proxy)", d, A),
    list("5 adjusted, euhaline only", d[S != "BRACKISH", ], A),
    list("6 adjusted, BUSCO C < 85% excluded", d[hq, ], A))
  for (x in mods) res[[length(res) + 1]] <- cbind(dataset = ds, aa = a, model = x[[1]], desc, fit1(x[[2]], x[[3]], tr))
  if (a == "A") {  # VIF of the adjusted-model predictors (full GC proxy)
    X <- data.frame(FW = d$FW, lenc = d$lenc, gcp = d$gcp)
    for (v in names(X)) vif[[length(vif) + 1]] <- data.frame(dataset = ds, n = nrow(X), term = v,
      VIF = 1 / (1 - summary(lm(X[[v]] ~ ., data = X[, setdiff(names(X), v), drop = FALSE]))$r.squared))
  }
}
res <- do.call(rbind, res)
res$diff_means <- res$mean_FW - res$mean_comp
res$q_BH <- ave(res$P, res$dataset, res$model, FUN = function(p) p.adjust(p, "BH"))
vif <- do.call(rbind, vif)
write.table(res, file.path(OUT, "aa_pgls_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(vif, file.path(OUT, "aa_vif.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("fits", nrow(res), "failed", sum(is.na(res$P)), "\n")
f <- res[res$aa %in% c("G", "A", "P", "Q"), ]
f$model <- substr(f$model, 1, 2)
print(format(f[, c("dataset", "aa", "model", "n", "lambda", "diff_means", "beta_FW", "ci_lo", "ci_hi", "P", "q_BH")],
             digits = 3), row.names = FALSE)
cat("\nq < 0.05 by dataset x model:\n")
s <- res[!is.na(res$q_BH) & res$q_BH < 0.05, ]
print(table(s$dataset, substr(s$model, 1, 2)))
print(s[, c("dataset", "aa", "model", "n", "beta_FW", "P", "q_BH")], row.names = FALSE, digits = 3)
print(vif, digits = 3)
