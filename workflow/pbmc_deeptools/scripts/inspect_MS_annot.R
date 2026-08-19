#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat); library(Signac)})
rds <- "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/first_eval/objects/merged_MS.rds"
x <- readRDS(rds)
ac <- x[["H3K27ac"]]; me <- x[["H3K27me3"]]

cat("===== H3K27ac meta.data columns =====\n"); print(colnames(ac@meta.data))
cat("\n===== H3K27me3 meta.data columns =====\n"); print(colnames(me@meta.data))

# print a table for every H3K27ac column that looks like an annotation
cand <- grep("predict|human|MS|CZI|annot|ident|celltype|cell.type|label",
             colnames(ac@meta.data), value = TRUE, ignore.case = TRUE)
cat("\n===== candidate annotation cols in H3K27ac =====\n"); print(cand)
for (cc in cand) {
  v <- ac@meta.data[[cc]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(as.character(v))
    if (length(u) <= 40) {
      cat("\n-- table(", cc, ") --\n"); print(table(v, useNA = "ifany"))
    } else cat("\n--", cc, ": ", length(u), "unique (numeric-like/too many)\n")
  }
}

# cell-name structure + matched-cell overlap (strip the mark token)
cat("\n===== cell names =====\n")
cat("ac[1:2]:", paste(head(colnames(ac), 2), collapse = " | "), "\n")
cat("me[1:2]:", paste(head(colnames(me), 2), collapse = " | "), "\n")
mkkey <- function(cn) gsub("H3K27ac_|H3K27me3_", "", cn)
ka <- mkkey(colnames(ac)); km <- mkkey(colnames(me))
cat("ac key[1:2]:", paste(head(ka, 2), collapse = " | "), "\n")
cat("me key[1:2]:", paste(head(km, 2), collapse = " | "), "\n")
cat("ac cells:", ncol(ac), " me cells:", ncol(me), "\n")
cat("dup keys  ac:", sum(duplicated(ka)), " me:", sum(duplicated(km)), "\n")
cat("me cells matched to an ac cell:", sum(km %in% ka), "of", length(km), "\n")
cat("ac cells matched to an me cell:", sum(ka %in% km), "of", length(ka), "\n")
