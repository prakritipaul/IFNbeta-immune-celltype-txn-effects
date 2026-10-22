# M1a collinearity check: variancePartition canCorPairs between the variance-partitioning terms,
# per cell type, on kept samples (the same sample metadata as R/m1a_varpart.R; no expression data needed).
#
# Usage: Rscript R/m1a_cancor.R <samples.csv> <out.csv>
# Out: out.csv  cell_type, term1, term2, cancor (canonical correlation, 0-1; full matrix, long format).

suppressPackageStartupMessages(library(variancePartition))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: m1a_cancor.R <samples.csv> <out.csv>")
samples <- read.csv(args[1], check.names = FALSE, colClasses = "character")
needed <- c("sample", "cell_type", "donor", "condition", "n_cells", "total_counts", "kept")
if (!all(needed %in% names(samples))) stop("samples.csv missing columns: ", paste(setdiff(needed, names(samples)), collapse = ", "))
if (!all(samples$kept %in% c("True", "False"))) stop("kept must be True/False")
samples <- samples[samples$kept == "True", ]
samples$log_n_cells <- log(as.numeric(samples$n_cells))
samples$log_total_counts <- log(as.numeric(samples$total_counts))
if (any(!is.finite(samples$log_n_cells) | !is.finite(samples$log_total_counts))) stop("n_cells and total_counts must be positive numbers")

form <- ~ donor + condition + log_n_cells + log_total_counts
out <- list()
for (ct in sort(unique(samples$cell_type))) {
  meta <- samples[samples$cell_type == ct, ]
  meta$donor <- factor(meta$donor)
  meta$condition <- factor(meta$condition)
  # Paired design: exactly 2 condition labels (read from the data), one row of each per kept donor
  pairs <- table(meta$donor, meta$condition)
  if (ncol(pairs) != 2 || any(pairs != 1))
    stop(ct, ": pairing violated: every kept donor needs exactly one sample in each of the 2 conditions (",
         paste(colnames(pairs), collapse = ", "), ")")
  cc <- canCorPairs(form, meta, showWarnings = TRUE)
  out[[ct]] <- data.frame(cell_type = ct, term1 = rep(rownames(cc), ncol(cc)),
                          term2 = rep(colnames(cc), each = nrow(cc)), cancor = as.vector(cc))
}
out <- do.call(rbind, out)
write.csv(out, args[2], row.names = FALSE)
