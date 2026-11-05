# Audit point 21 (2026-11-04): sensitivity of the 20-family PGLS to the proteome-mode BUSCO rerun (Table S1,
# "BUSCO rerun C (%)"). Re-uses the set-up and fit1() of ../genefam20_pgls.R unchanged (evaluated up to the model loop),
# then excludes taxa with rerun BUSCO C < 85% (unadjusted and adjusted with rerun BUSCO C as covariate). BH across 20.
src <- readLines(file.path(Sys.getenv("IZ_ROOT", "."), "_R1_work/analysis/genefam20/genefam20_pgls.R"), encoding = "UTF-8")
eval(parse(text = src[1:(which(startsWith(src, "res <- list()"))[1] - 1)]))
s1 <- read.delim(file.path(Sys.getenv("IZ_ROOT", "."), "SUBMIT_R1_revised_20261024/Supplementary_Table_S1_retained_bivalve_genomes_R1.tsv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
rr <- setNames(as.numeric(s1[["BUSCO rerun C (%)"]]), s1$ID)
d0$busco_rr <- rr[d0$id]
stopifnot(!any(is.na(d0$busco_rr)))
Ar <- lc ~ FW + busco_rr + l_n50 + l_nprot
res <- list()
for (lv in LEVELS) for (fm in FAMS) {
  d <- d0[d0$family == fm, ]
  if (lv != "OG_L1") d <- d[d$pfam_route_complete == 1, ]
  d$cnt <- d[[lv]]; d$lc <- log1p(d$cnt)
  hr <- d$busco_rr >= 85
  mods <- list(list("1 main", d, U, tr),
               list("R4 rerun BUSCO C < 85% excluded", d[hr, ], U, tr),
               list("R6 adjusted (rerun BUSCO C, log10 N50, log10 n_prot)", d, Ar, tr),
               list("R8 adjusted, rerun BUSCO C < 85% excluded", d[hr, ], Ar, tr))
  for (x in mods) res[[length(res) + 1]] <- cbind(level = lv, family = fm, model = x[[1]], fit1(x[[2]], x[[3]], x[[4]]))
}
res <- do.call(rbind, res)
q <- rep(NA_real_, nrow(res))
for (g in unique(paste(res$level, res$model))) {
  i <- which(paste(res$level, res$model) == g); m <- i[res$family[i] %in% FAMS20]; q[m] <- p.adjust(res$P[m], "BH") }
res$q_BH <- q
write.table(res, "busco_rerun_sens_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("taxa rerun C<85:", paste(sort(unique(d0$id[d0$busco_rr < 85])), collapse = ","), "\n")
k <- c("family", "model", "n", "n_FW", "lambda", "ratio_FW", "ratio_lo", "ratio_hi", "P", "q_BH")
print(res[res$level == "L2" & res$family %in% c("AQP", "CTS", "PDI"), k], digits = 3, row.names = FALSE)
cat("L2 families q<0.05 per model:\n"); print(tapply(res$q_BH[res$level=="L2"] < 0.05, res$model[res$level=="L2"], sum, na.rm=TRUE))
