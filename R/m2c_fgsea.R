# M2c preranked GSEA per cell type (contract: audit/contracts/M2c_programs.md; decisions G1-G4).
#
# Usage: Rscript R/m2c_fgsea.R <ranks.csv> <gmt> [<gmt> ...] <out.csv> <seed>
# In:  ranks.csv  cell_type, gene, stat (one row per tested gene; stat = Wald statistic log2FC / lfcSE).
#      gmt        one or more GMT files (name<TAB>description<TAB>genes...); collection = file name without .gmt.
#                 Set genes are matched to ranked genes by exact symbol only (fgsea does no fuzzy matching).
# Out: out.csv    cell_type, collection, pathway, size, ES, NES, pval, padj, log2err, leadingEdge (";"-joined).
#                 Only sets with 15-500 genes among the cell type's ranked genes are tested.
# Routine, per cell type (alphabetical):
#   1) stats = named vector of stat, in the input order (Python sorts by stat desc, then gene, so ties are fixed)
#   2) set.seed(seed); fgseaMultilevel on all sets of all collections together
#   3) padj = BH over all sets of that cell type (both collections), recomputed explicitly
# Invalid input (NA/duplicate/non-finite stat, duplicate set name) stops with a non-zero exit before writing.

suppressPackageStartupMessages(library(fgsea))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: m2c_fgsea.R <ranks.csv> <gmt> [<gmt> ...] <out.csv> <seed>")
n <- length(args)
seed <- as.integer(args[n])
if (is.na(seed)) stop("seed must be an integer")
out_path <- args[n - 1]
gmt_paths <- args[2:(n - 2)]

ranks <- read.csv(args[1], check.names = FALSE, colClasses = c(cell_type = "character", gene = "character", stat = "numeric"))
if (!all(c("cell_type", "gene", "stat") %in% names(ranks))) stop("ranks.csv needs cell_type, gene, stat")
bad <- !is.finite(ranks$stat)
if (any(bad)) stop("non-finite stat: ", paste(ranks$cell_type[bad], ranks$gene[bad], collapse = "; "))
dup <- duplicated(ranks[, c("cell_type", "gene")])
if (any(dup)) stop("duplicated cell_type x gene: ", paste(ranks$cell_type[dup], ranks$gene[dup], collapse = "; "))

pathways <- list(); collection_of <- character()
for (path in gmt_paths) {
  sets <- gmtPathways(path)
  collection <- sub("\\.gmt$", "", basename(path))
  clash <- intersect(names(sets), names(pathways))
  if (length(clash)) stop("set names in more than one GMT: ", paste(clash, collapse = ", "))
  pathways <- c(pathways, sets)
  collection_of[names(sets)] <- collection
}

out <- list()
for (ct in sort(unique(ranks$cell_type))) {
  r <- ranks[ranks$cell_type == ct, ]
  stats <- setNames(r$stat, r$gene)
  set.seed(seed)
  res <- fgseaMultilevel(pathways = pathways, stats = stats, sampleSize = 101, minSize = 15, maxSize = 500,
                         eps = 0, scoreType = "std", nproc = 1, gseaParam = 1, BPPARAM = NULL,
                         nPermSimple = 1000, absEps = NULL)
  if (nrow(res) == 0) next
  if (any(is.na(res$pval))) stop("NA p-value in ", ct, ": ", paste(res$pathway[is.na(res$pval)], collapse = ", "))
  out[[ct]] <- data.frame(
    cell_type = ct, collection = unname(collection_of[res$pathway]), pathway = res$pathway, size = res$size,
    ES = res$ES, NES = res$NES, pval = res$pval,
    padj = p.adjust(res$pval, method = "BH"),  # BH within cell type, over all sets of all collections
    log2err = res$log2err,
    leadingEdge = vapply(res$leadingEdge, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE)
}
cols <- c("cell_type", "collection", "pathway", "size", "ES", "NES", "pval", "padj", "log2err", "leadingEdge")
result <- if (length(out)) do.call(rbind, out) else setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
write.csv(result, out_path, row.names = FALSE)
