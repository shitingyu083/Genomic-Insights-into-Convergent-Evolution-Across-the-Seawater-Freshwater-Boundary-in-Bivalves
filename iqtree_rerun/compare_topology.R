# R1 phylogeny option A: compare the R0 IQ-TREE species tree with the OrthoFinder species tree
# that every archived analysis used (ASR/Fsite, PGLS, chronos).
# Usage: Rscript compare_topology.R <iqtree.treefile|contree> [tip_map.tsv: iqtree_label<TAB>code]
# Writes to iqtree_rerun/out/: topology_report.txt, split_table.tsv, iqtree_species_tree_codes.nwk
suppressPackageStartupMessages(library(ape))
a <- commandArgs(TRUE)
stopifnot(length(a) >= 1, file.exists(a[1]))
W <- file.path(Sys.getenv("IZ_ROOT", "."), "_R1_work")
REF <- file.path(W, "server_data/x/orthofinder1229/workflow/convergence_analysis_workflow/step2_convergence/Fsite_method/00_input/SpeciesTree_rooted_node_labels.nwk")
OUT <- file.path(W, "analysis/iqtree_rerun/out"); dir.create(OUT, FALSE, TRUE)

ref <- read.tree(REF); ref$tip.label <- sub("[.]pep$", "", ref$tip.label); ref$node.label <- NULL
iq <- read.tree(a[1])
lab <- sub("[.]pep$", "", iq$tip.label)
if (length(a) >= 2) {
  m <- read.delim(a[2], header = FALSE, stringsAsFactors = FALSE)
  hit <- match(lab, m[[1]]); lab[!is.na(hit)] <- m[[2]][hit[!is.na(hit)]]
}
iq$tip.label <- lab
miss_iq <- setdiff(ref$tip.label, iq$tip.label); miss_ref <- setdiff(iq$tip.label, ref$tip.label)
if (length(miss_iq) || length(miss_ref)) {
  cat("Tip sets differ. In OrthoFinder tree only:", miss_iq, "\nIn IQ-TREE only:", miss_ref,
      "\nSupply a tip_map.tsv (iqtree_label<TAB>code) and rerun.\n"); quit(status = 2)
}
iq <- root(iq, "Lgi", resolve.root = TRUE)

# node support: IQ-TREE labels are "SH-aLRT/UFBoot" (with -alrt) or "UFBoot"
sup <- iq$node.label
spl <- function(tr) { p <- prop.part(tr); lapply(p, function(x) sort(attr(p, "labels")[x])) }
key <- function(s, all) { s2 <- sort(setdiff(all, s)); if (length(s) <= length(s2)) paste(s, collapse = ",") else paste(s2, collapse = ",") }
all <- sort(ref$tip.label)
k_ref <- vapply(spl(ref), key, "", all = all); k_iq <- vapply(spl(iq), key, "", all = all)
n <- Ntip(ref)
nontriv <- function(k) k[vapply(strsplit(k, ","), length, 1L) > 1]
kr <- unique(nontriv(k_ref)); ki <- unique(nontriv(k_iq))
tab <- data.frame(split = ki, in_orthofinder = ki %in% kr,
                  iqtree_support = if (is.null(sup)) NA else sup[match(ki, k_iq)], stringsAsFactors = FALSE)
write.table(tab, file.path(OUT, "split_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rf <- length(setdiff(kr, ki)) + length(setdiff(ki, kr))
sink(file.path(OUT, "topology_report.txt"))
cat("IQ-TREE file:", normalizePath(a[1]), "md5", tools::md5sum(a[1]), "\n")
cat("OrthoFinder file:", REF, "md5", tools::md5sum(REF), "\n")
cat("tips:", n, " rooted on Lgi\n")
cat("non-trivial splits: OrthoFinder", length(kr), " IQ-TREE", length(ki), "\n")
cat("Robinson-Foulds distance (unrooted, symmetric difference):", rf, "\n")
cat("identical topology:", rf == 0, "\n")
cat("\nIQ-TREE splits absent from OrthoFinder tree:\n"); print(tab[!tab$in_orthofinder, ], row.names = FALSE)
cat("\nOrthoFinder splits absent from IQ-TREE tree:\n"); print(setdiff(kr, ki))
cat("\nIQ-TREE node support (all splits):\n"); print(tab, row.names = FALSE)
sink()
cat(readLines(file.path(OUT, "topology_report.txt"))[1:6], sep = "\n")

# tree for the PGLS rerun: codes as tip labels, IQ-TREE branch lengths, no node labels
iq0 <- iq; iq0$node.label <- NULL
write.tree(iq0, file.path(OUT, "iqtree_species_tree_codes.nwk"))
