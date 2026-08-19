#!/usr/bin/env Rscript
# Correlation of ENCODE bulk vs scPBMC pseudobulk over 5kb bins, as R corrplots.
# Produces BOTH: Spearman (rank-based, on raw signal) and Pearson (on log1p signal,
# the standard for genomic coverage). Reads deeptools multiBigwigSummary --outRawCounts.
suppressPackageStartupMessages(library(corrplot))

out  <- "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/encode_correlation"
tabf <- file.path(out, "counts_5kb.tab")

# --- read the deeptools raw-counts table (header: #'chr' 'start' 'end' 'label'...) ---
hdr    <- readLines(tabf, n = 1)
labels <- strsplit(gsub("#|'", "", hdr), "\t")[[1]]
d      <- read.table(tabf, header = FALSE, skip = 1, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
colnames(d) <- labels

mat <- as.matrix(d[, -(1:3)])                 # drop chr/start/end
mode(mat) <- "numeric"
keep <- rowSums(is.na(mat)) == 0 & rowSums(mat) > 0
mat  <- mat[keep, , drop = FALSE]
cat("bins used:", nrow(mat), "| tracks:", ncol(mat), "\n")

pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(200)
draw_corr <- function(cormat, file_base, title) {
  d <- function() corrplot(cormat, method = "color", order = "hclust", col = pal,
                           addCoef.col = "grey25", number.cex = 0.7,
                           tl.col = "black", tl.srt = 45, tl.cex = 0.8,
                           mar = c(0, 0, 1, 0), title = title)
  png(paste0(file_base, ".png"), width = 2000, height = 2000, res = 220); d(); dev.off()
  pdf(paste0(file_base, ".pdf"), width = 9, height = 9);                  d(); dev.off()
}

# --- Spearman (rank-based, on raw signal) ---
cs <- cor(mat, method = "spearman")
write.csv(round(cs, 4), file.path(out, "correlation_matrix_spearman.csv"))
ttl_s <- "ENCODE bulk vs scPBMC pseudobulk (Spearman, 5kb bins)"
draw_corr(cs, file.path(out, "corrplot"),          ttl_s)   # keep corrplot.png (Snakefile output)
draw_corr(cs, file.path(out, "corrplot_spearman"), ttl_s)

# --- Pearson (on log1p signal -- standard for genomic coverage) ---
cp <- cor(log1p(mat), method = "pearson")
write.csv(round(cp, 4), file.path(out, "correlation_matrix_pearson.csv"))
draw_corr(cp, file.path(out, "corrplot_pearson"),
          "ENCODE bulk vs scPBMC pseudobulk (Pearson on log1p, 5kb bins)")

cat("Wrote Spearman + Pearson matrices and corrplots to:", out, "\n")
