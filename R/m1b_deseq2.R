# M1b pseudobulk DE (contract: audit/contracts/M1b_de.md; decisions O1-O6).
#
# Usage: Rscript R/m1b_deseq2.R <counts.csv> <samples.csv> <paired|unpaired> <results_out.csv> <fits_out.csv> <min_basemean>
# In:  counts.csv  genes x samples (first column = gene), integer counts, columns = sample ids.
#      samples.csv sample, cell_type, donor, condition, kept (other columns ignored).
#      Condition labels are read from the data: exactly 2, one of them "ctrl"/"control" (the reference).
# Out: results_out.csv  cell_type, gene, baseMean, log2FC, lfcSE, pvalue, qvalue, tested, cooks_flag, log2FC_shrunk
#        (every gene of every cell type; log2FC = treatment vs control, unshrunk MLE).
#        tested     = baseMean > min_basemean; qvalue only for tested genes with a p-value, else NA.
#        cooks_flag = DESeq2 replaced the gene's outlier count (mcols(dds)$replace, cell types with >= 7
#                     samples per condition) OR its p-value was set NA by results()' Cook's filter
#                     (p-value NA with cooksCutoff = TRUE but not with cooksCutoff = FALSE).
#      fits_out.csv     cell_type, fit_type_used ("parametric", or "local" when DESeq2 substituted it), n_samples, n_genes_tested.
# Routine, per cell type, kept samples only:
#   1) DESeqDataSet, design ~ donor + condition (paired) or ~ condition (unpaired); control = reference level
#   2) estimateSizeFactors(type = "ratio"); DESeq(Wald, parametric, betaPrior = FALSE, minReplicatesForReplace = 7)
#   3) results(cooksCutoff = TRUE (default cutoff), independentFiltering = FALSE: padj is not used and not written)
#   4) log2FC_shrunk = lfcShrink(type = "normal") (O5; V9b only)
#   5) qvalue(lambda = 0.05..0.95 by 0.05, smoother, smooth.df = 3) on tested genes' p-values (O3, per cell type)
# Invalid input, a failed fit, or an unexplained NA stops with a non-zero exit status before anything is written.

suppressPackageStartupMessages({
  library(DESeq2)
  library(qvalue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("usage: m1b_deseq2.R <counts.csv> <samples.csv> <paired|unpaired> <results_out.csv> <fits_out.csv> <min_basemean>")
design_name <- args[3]
if (!design_name %in% c("paired", "unpaired")) stop("design must be 'paired' or 'unpaired'")
min_basemean <- as.numeric(args[6])
if (!is.finite(min_basemean)) stop("min_basemean must be a finite number")

counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
samples <- read.csv(args[2], check.names = FALSE, colClasses = "character")

needed <- c("sample", "cell_type", "donor", "condition", "kept")
if (!all(needed %in% names(samples))) stop("samples.csv missing columns: ", paste(setdiff(needed, names(samples)), collapse = ", "))
if (anyDuplicated(samples$sample)) stop("duplicated sample ids in samples.csv")
if (anyDuplicated(rownames(counts))) stop("duplicated gene names in counts.csv")
if (!all(samples$kept %in% c("True", "False"))) stop("kept must be True/False")
samples <- samples[samples$kept == "True", ]
if (!all(samples$sample %in% colnames(counts))) stop("kept sample ids missing from counts.csv")
if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)))
  stop("counts must be finite non-negative integers")
storage.mode(counts) <- "integer"

labels <- sort(unique(samples$condition))
control <- intersect(labels, c("ctrl", "control"))
if (length(labels) != 2 || length(control) != 1)
  stop("need exactly one control label (ctrl/control) and one treatment; found ", paste(labels, collapse = ", "))
treatment <- setdiff(labels, control)
rownames(samples) <- samples$sample

FALLBACK <- "local regression fit was automatically substituted"
results_list <- list(); fits <- list()
for (ct in sort(unique(samples$cell_type))) {
  meta <- samples[samples$cell_type == ct, ]
  pairs <- table(factor(meta$donor), factor(meta$condition, levels = c(control, treatment)))
  if (any(pairs != 1))
    stop(ct, ": pairing violated: every kept donor needs exactly one ", control, " and one ", treatment, " sample")
  meta$donor <- factor(meta$donor)
  meta$condition <- relevel(factor(meta$condition), ref = control)
  design <- if (design_name == "paired") ~ donor + condition else ~ condition
  dds <- DESeqDataSetFromMatrix(counts[, meta$sample, drop = FALSE], colData = meta, design = design)
  dds <- estimateSizeFactors(dds, type = "ratio")
  fit_used <- "parametric"
  # Record DESeq2's automatic parametric -> local fallback; the message is still printed (not muffled)
  dds <- withCallingHandlers(
    DESeq(dds, test = "Wald", fitType = "parametric", betaPrior = FALSE, minReplicatesForReplace = 7,
          parallel = FALSE, quiet = FALSE),
    message = function(m) if (grepl(FALLBACK, conditionMessage(m), fixed = TRUE)) fit_used <<- "local")
  contrast <- c("condition", treatment, control)
  res <- results(dds, contrast = contrast, cooksCutoff = TRUE, independentFiltering = FALSE)
  res_nocooks <- results(dds, contrast = contrast, cooksCutoff = FALSE, independentFiltering = FALSE)
  all_zero <- mcols(dds)$allZero
  if (any(is.na(res$log2FoldChange) & !all_zero)) stop(ct, ": log2FC is NA for a gene with non-zero counts")
  if (any(is.na(res_nocooks$pvalue) & !all_zero)) stop(ct, ": p-value NA for a non-zero gene without Cook's filtering (unexplained)")
  cooks_na <- is.na(res$pvalue) & !all_zero
  replaced <- if (is.null(mcols(dds)$replace)) rep(FALSE, nrow(dds)) else mcols(dds)$replace %in% TRUE
  shrunk <- lfcShrink(dds, contrast = contrast, type = "normal", quiet = TRUE)
  stopifnot(identical(rownames(shrunk), rownames(res)))

  tested <- res$baseMean > min_basemean
  q <- rep(NA_real_, nrow(res))
  use <- tested & !is.na(res$pvalue)
  q[use] <- qvalue(res$pvalue[use], lambda = seq(0.05, 0.95, 0.05), pi0.method = "smoother", smooth.df = 3)$qvalues

  results_list[[ct]] <- data.frame(cell_type = ct, gene = rownames(res), baseMean = res$baseMean,
                                   log2FC = res$log2FoldChange, lfcSE = res$lfcSE, pvalue = res$pvalue, qvalue = q,
                                   tested = tested, cooks_flag = replaced | cooks_na,
                                   log2FC_shrunk = shrunk$log2FoldChange, row.names = NULL)
  fits[[ct]] <- data.frame(cell_type = ct, fit_type_used = fit_used, n_samples = ncol(dds), n_genes_tested = sum(tested))
}
fits <- do.call(rbind, fits)
print(fits, row.names = FALSE)
write.csv(do.call(rbind, results_list), args[4], row.names = FALSE)
write.csv(fits, args[5], row.names = FALSE)
