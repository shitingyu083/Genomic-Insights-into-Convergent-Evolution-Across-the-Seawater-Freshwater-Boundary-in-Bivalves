# R1 revision: ultrametric species tree for CAFE5 from the OrthoFinder rooted species tree (STAG/STRIDE).
# Penalised-likelihood dating (ape::chronos, correlated model, lambda = 1) with the root fixed at 100 relative
# time units (no fossil calibration; CAFE5 lambda is therefore per relative time unit).
library(ape)
a <- commandArgs(TRUE)
t <- read.tree(a[1])
t$tip.label <- sub("\\.pep$", "", t$tip.label)
t$node.label <- NULL
set.seed(1)
cal <- makeChronosCalib(t, node = "root", age.min = 100, age.max = 100)
u <- chronos(t, lambda = 1, model = "correlated", calibration = cal)
class(u) <- "phylo"
stopifnot(is.ultrametric(u, tol = 1e-6), is.binary(u))
write.tree(u, a[2], digits = 8)
cat("tips", Ntip(u), "root depth", max(node.depth.edgelength(u)), "\n")
