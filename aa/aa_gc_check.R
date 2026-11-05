# R1 revision, R2 comment 12 (nucleotide composition) / P76-P77 / Fig. 5: CDS-level check for the 12 species whose
# CDS FASTA is on the server (server_data/r1rev/aa_gc/aa_gc_species.tsv from remote/s27_aa_gc.py, run in
# ~/R1_revision_20261024/aa_gc; CDS with length %% 3 != 0 or internal stops skipped, terminal stop stripped).
# (1) translated-CDS amino-acid frequencies vs the OrthoFinder proteome used in aa_pgls.R (ALL_ANY f_*);
# (2) GC, GC3, GC4 freshwater vs comparators: exact Wilcoxon and PGLS (same tree and fitting rule as aa_pgls.R);
# (3) Spearman correlation of GC3/GC4 with the proteome GC proxy ln(GARP) - ln(FYMINK) used as covariate in aa_pgls.R
#     and with CLR(G/A/P/Q);
# (4) CLR(G/A/P/Q) ~ FW and ~ FW + GC4 within the same subsets (descriptive; n = 9-12).
# Sets: S10 = the 10 species whose CDS header count equals the proteome count (5 FW: Cfl Dpo Dro Lfo Pst;
# comparators Cgi Csi Tna Sco brackish, Sbr euhaline); S9 = S10 without Lfo (20,802 of 60,717 CDS with internal stops);
# S12 = S10 + Mme, Mye (euhaline; CDS files hold fewer records than the proteome).
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

fit1 <- function(dd, fml, tree) {  # identical to aa_pgls.R
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
  if (inherits(m, "try-error")) return(data.frame(lambda = NA, beta_FW = NA, ci_lo = NA, ci_hi = NA, P = NA, other = "gls failed"))
  tab <- summary(m)$tTable
  ci <- tab["FW", "Value"] + c(-1, 1) * qt(0.975, nrow(dd) - nrow(tab)) * tab["FW", "Std.Error"]
  ot <- if (nrow(tab) > 2) paste(sprintf("%s=%.4g(P=%.3g)", rownames(tab)[-(1:2)], tab[-(1:2), "Value"],
                                         tab[-(1:2), "p-value"]), collapse = ";") else ""
  data.frame(lambda = lam, beta_FW = tab["FW", "Value"], ci_lo = ci[1], ci_hi = ci[2], P = tab["FW", "p-value"], other = ot)
}

# (1) translated CDS vs proteome
trc <- as.matrix(d[, paste0("tr_", AA)]); ftr <- trc / rowSums(trc)
fpr <- as.matrix(d[, paste0("f_", AA)])
dif <- abs(ftr - fpr)
t1 <- data.frame(species = d$species, salinity3 = d$salinity3, header_count_matches = d$matches_proteome_count,
                 n_cds = d$n_seq, n_prot = d$n_prot, n_cds_used = d$n_used, used_over_prot = d$n_used / d$n_prot,
                 TVD = 0.5 * rowSums(dif), max_abs_diff = apply(dif, 1, max), max_aa = AA[apply(dif, 1, which.max)])
write.table(t1, file.path(OUT, "aa_gc_trcheck.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

SETS <- list(S10 = d$species[d$matches_proteome_count == 1], S9 = setdiff(d$species[d$matches_proteome_count == 1], "Lfo"),
             S12 = d$species)
# (2) GC metrics FW vs comparators; (4) CLR(G/A/P/Q) ~ FW (+ GC4)
t2 <- list()
for (s in names(SETS)) {
  dd <- d[d$species %in% SETS[[s]], ]
  for (v in c("GC", "GC3", "GC4", paste0("clr_", c("G", "A", "P", "Q")))) {
    dd$y <- dd[[v]]
    w <- wilcox.test(dd$y[dd$FW == 1], dd$y[dd$FW == 0], exact = TRUE)
    for (f in c("FW", "FW+GC4")) {
      if (f == "FW+GC4" && !grepl("^clr_", v)) next
      fml <- if (f == "FW") y ~ FW else y ~ FW + GC4
      p <- fit1(dd, fml, tr)
      t2[[length(t2) + 1]] <- data.frame(set = s, response = v, model = f, n = nrow(dd), n_FW = sum(dd$FW),
        mean_FW = mean(dd$y[dd$FW == 1]), mean_comp = mean(dd$y[dd$FW == 0]), wilcox_W = unname(w$statistic),
        wilcox_P = w$p.value, p)
    }
  }
}
t2 <- do.call(rbind, t2)
write.table(t2, file.path(OUT, "aa_gc_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# (3) proxy validation and GC-composition correlations
t3 <- list()
for (s in names(SETS)) {
  dd <- d[d$species %in% SETS[[s]], ]
  for (x in c("GC3", "GC4")) for (y in c("gcp", paste0("clr_", c("G", "A", "P", "Q")))) {
    k <- cor.test(dd[[x]], dd[[y]], method = "spearman", exact = TRUE)
    t3[[length(t3) + 1]] <- data.frame(set = s, x = x, y = y, n = nrow(dd), rho = unname(k$estimate), P = k$p.value)
  }
}
t3 <- do.call(rbind, t3)
write.table(t3, file.path(OUT, "aa_gc_proxy_cor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

options(width = 200)
cat("## (1) translated CDS vs proteome\n"); print(t1, digits = 4, row.names = FALSE)
cat("\n## (2)/(4) FW vs comparators\n")
print(transform(t2[, c("set", "response", "model", "n", "n_FW", "mean_FW", "mean_comp", "wilcox_P", "lambda", "beta_FW",
                       "ci_lo", "ci_hi", "P", "other")]), digits = 3, row.names = FALSE)
cat("\n## (3) Spearman\n"); print(t3, digits = 3, row.names = FALSE)
