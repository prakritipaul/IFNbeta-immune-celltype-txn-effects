# M1a variance partitioning (contract: audit/contracts/M1a_pseudobulk.md; decisions A2-A3).
#
# Usage: Rscript R/m1a_varpart.R <counts.csv> <samples.csv> <out.csv> <min_basemean>
# In:  counts.csv  genes x samples (first column = gene), integer counts, columns = sample ids.
#      samples.csv sample, cell_type, donor, condition, n_cells, total_counts, kept.
# Out: out.csv  cell_type, gene, baseMean, donor, condition, log_n_cells, log_total_counts, Residuals
#      (variance fractions; each row sums to 1).
#      <out>_vst.csv   cell_type, gene, sample, vst (the transformed values the model was fit on; long format).
#      <out>_fits.csv  cell_type, fit_type_used (DESeq2 dispersion trend: "parametric", or "local" when
#                      DESeq2 substituted it automatically; decision: allowed, but logged).
# Routine, per cell type, kept samples only:
#   1) DESeq2 median-of-ratios size factors on all genes
#   2) keep genes with baseMean (mean normalised count) > min_basemean
#   3) varianceStabilizingTransformation(blind = TRUE, fitType = "parametric") on the kept genes (A3 revised)
#   4) variancePartition: ~ (1|donor) + (1|condition) + log(n_cells) + log(total_counts).
# Invalid input stops with a non-zero exit status before anything is written.

suppressPackageStartupMessages({
  library(DESeq2)
  library(variancePartition)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: m1a_varpart.R <counts.csv> <samples.csv> <out.csv> <min_basemean>")
min_basemean <- as.numeric(args[4])
if (!is.finite(min_basemean)) stop("min_basemean must be a finite number")
if (!grepl("\\.csv$", args[3])) stop("out path must end in .csv (the _vst and _fits files are named from it)")

counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
samples <- read.csv(args[2], check.names = FALSE, colClasses = "character")

needed <- c("sample", "cell_type", "donor", "condition", "n_cells", "total_counts", "kept")
if (!all(needed %in% names(samples))) stop("samples.csv missing columns: ", paste(setdiff(needed, names(samples)), collapse = ", "))
if (anyDuplicated(samples$sample)) stop("duplicated sample ids in samples.csv")
if (!setequal(samples$sample, colnames(counts))) stop("sample ids differ between counts.csv and samples.csv")
if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)))
  stop("counts must be finite non-negative integers")
# pandas writes True/False; anything else is an error, not a guess
if (!all(samples$kept %in% c("True", "False"))) stop("kept must be True/False")
samples$kept <- samples$kept == "True"
samples$n_cells <- as.numeric(samples$n_cells)
samples$total_counts <- as.numeric(samples$total_counts)
if (any(!is.finite(samples$n_cells) | samples$n_cells <= 0 | !is.finite(samples$total_counts) | samples$total_counts <= 0))
  stop("n_cells and total_counts must be positive numbers")

samples <- samples[samples$kept, ]
rownames(samples) <- samples$sample
samples$log_n_cells <- log(samples$n_cells)
samples$log_total_counts <- log(samples$total_counts)
storage.mode(counts) <- "integer"

form <- ~ (1 | donor) + (1 | condition) + log_n_cells + log_total_counts
results <- list(); vst_long <- list(); fits <- list()
FALLBACK <- "local regression fit was automatically substituted"
for (ct in sort(unique(samples$cell_type))) {
  meta <- samples[samples$cell_type == ct, ]
  meta$donor <- factor(meta$donor)
  meta$condition <- factor(meta$condition)
  # Paired design: exactly 2 condition labels (read from the data), one row of each per kept donor
  pairs <- table(meta$donor, meta$condition)
  if (ncol(pairs) != 2 || any(pairs != 1))
    stop(ct, ": pairing violated: every kept donor needs exactly one sample in each of the 2 conditions (",
         paste(colnames(pairs), collapse = ", "), ")")
  dds <- DESeqDataSetFromMatrix(counts[, meta$sample, drop = FALSE], colData = meta, design = ~ 1)
  dds <- estimateSizeFactors(dds, type = "ratio")
  base_mean <- rowMeans(counts(dds, normalized = TRUE))
  dds <- dds[base_mean > min_basemean, ]
  fit_used <- "parametric"
  # Record DESeq2's automatic parametric -> local fallback; the message is still printed (not muffled)
  vsd <- withCallingHandlers(
    varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric"),
    message = function(m) if (grepl(FALLBACK, conditionMessage(m), fixed = TRUE)) fit_used <<- "local")
  fits[[ct]] <- data.frame(cell_type = ct, fit_type_used = fit_used)
  v <- assay(vsd)
  vst_long[[ct]] <- data.frame(cell_type = ct, gene = rep(rownames(v), ncol(v)),
                               sample = rep(colnames(v), each = nrow(v)), vst = as.vector(v))
  vp <- fitExtractVarPartModel(v, form, meta, REML = FALSE, useWeights = TRUE,
                               showWarnings = TRUE, quiet = FALSE, BPPARAM = SerialParam())
  vp <- as.data.frame(vp)
  results[[ct]] <- data.frame(cell_type = ct, gene = rownames(vp), baseMean = base_mean[rownames(vp)],
                              vp[, c("donor", "condition", "log_n_cells", "log_total_counts", "Residuals")],
                              check.names = FALSE, row.names = NULL)
}
out <- do.call(rbind, results)
fractions <- out[, c("donor", "condition", "log_n_cells", "log_total_counts", "Residuals")]
if (any(!is.finite(as.matrix(fractions))) || any(abs(rowSums(fractions) - 1) > 1e-6))
  stop("variance fractions are not finite or do not sum to 1")
fits <- do.call(rbind, fits)
print(fits, row.names = FALSE)
write.csv(out, args[3], row.names = FALSE)
write.csv(do.call(rbind, vst_long), sub("\\.csv$", "_vst.csv", args[3]), row.names = FALSE)
write.csv(fits, sub("\\.csv$", "_fits.csv", args[3]), row.names = FALSE)
