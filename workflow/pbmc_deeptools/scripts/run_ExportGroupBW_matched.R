#!/usr/bin/env Rscript
# Per-cell-type bigwigs for a multimodal MS object, label shared across marks:
#   1. group the REFERENCE mark (default H3K27ac) by an annotation column
#      (default predicted.id_Human_MS) and export one RC bigwig per group.
#   2. TRANSFER that label to the matched cells of the QUERY mark (default
#      H3K27me3) via a modality-agnostic cell key (sample.patient.16bp cell
#      barcode -- shared across modalities, differs only in the mark token),
#      keep only matched cells, and export by the SAME label.
# Loads the (large) rds once. Sequential export (no plan(multisession); that
# broke export.bw + OOM'd on the big group). Labels sanitized for filenames.
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(GenomicRanges)
  library(rtracklayer); library(pbapply)
  library(future); library(future.apply)   # ExportGroupBW calls nbrOfWorkers()
})
source("/home/mattia/Multi_nanoCTRNA/workflow/pbmc_deeptools/scripts/ExportGroupBW.R")

a <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) { i <- which(a == paste0("--", k)); if (length(i)) a[i + 1] else d }
rds      <- getarg("rds")
refmark  <- getarg("refmark", "H3K27ac")
qmark    <- getarg("qmark", "H3K27me3")
group    <- getarg("group", "predicted.id_Human_MS")
outref   <- getarg("outref")
outq     <- getarg("outq")
tileSize <- as.integer(getarg("tilesize", "100"))
cutoff   <- as.integer(getarg("cutoff", "4"))
minCells <- as.integer(getarg("mincells", "5"))
sanitize <- toupper(getarg("sanitize", "TRUE")) %in% c("TRUE", "T")
stopifnot(!is.null(rds), !is.null(outref), !is.null(outq))

san <- function(v) { v <- as.character(v); if (sanitize) { v <- gsub("[^A-Za-z0-9]+", "_", v); v <- gsub("^_|_$", "", v) }; v }
mkkey <- function(cn) gsub(paste0(refmark, "_|", qmark, "_"), "", cn)  # strip mark token

cat("Loading", rds, "\n")
x <- readRDS(rds)
ref <- x[[refmark]]; qy <- x[[qmark]]
stopifnot(group %in% colnames(ref@meta.data))
options(future.globals.maxSize = 8 * 1024^3)

# ---- 1. reference mark: group by annotation ----
ref@meta.data[[group]] <- san(ref@meta.data[[group]])
DefaultAssay(ref) <- "bins"
cat("\n== REF", refmark, "groups ==\n"); print(table(ref@meta.data[[group]], useNA = "ifany"))
dir.create(outref, recursive = TRUE, showWarnings = FALSE)
ExportGroupBW(ref, assay = "bins", group.by = group, normMethod = "RC",
              tileSize = tileSize, minCells = minCells, cutoff = cutoff,
              outdir = outref, verbose = TRUE)

# ---- 2. query mark: transfer label to matched cells, export ----
ka  <- mkkey(colnames(ref)); km <- mkkey(colnames(qy))
map <- setNames(ref@meta.data[[group]], ka)   # ka has no dups (checked)
lab <- unname(map[km])                         # NA where no matching ref cell
cat("\nmatched", qmark, "cells:", sum(!is.na(lab)), "of", ncol(qy), "\n")
qy@meta.data[[group]] <- lab
qy <- subset(qy, cells = colnames(qy)[!is.na(lab)])
qy@meta.data[[group]] <- as.character(qy@meta.data[[group]])
DefaultAssay(qy) <- "bins"
cat("\n== QUERY", qmark, "groups (matched cells) ==\n"); print(table(qy@meta.data[[group]], useNA = "ifany"))
dir.create(outq, recursive = TRUE, showWarnings = FALSE)
ExportGroupBW(qy, assay = "bins", group.by = group, normMethod = "RC",
              tileSize = tileSize, minCells = minCells, cutoff = cutoff,
              outdir = outq, verbose = TRUE)

cat("\nDONE\n")
cat("ref bigwigs:\n"); print(list.files(outref, pattern = "\\.bw$"))
cat("query bigwigs:\n"); print(list.files(outq, pattern = "\\.bw$"))
