#!/usr/bin/env Rscript
# Export one bigwig per group (e.g. transferred cell type) from a merged
# multiNanoCT Seurat object, using Signac PR#1383 ExportGroupBW (RC-normalized,
# fragment-level, minCells/cutoff filtered). Fragments already live in the
# ChromatinAssay of the merged object (one Fragment obj per experiment), so no
# fragment rebuild is needed -- SplitFragments iterates them.
#
# Usage:
#   Rscript run_ExportGroupBW.R \
#     --rds   /.../merged_MS.rds \
#     --mark  H3K27ac \            # element of the rds list
#     --assay bins \               # ChromatinAssay name
#     --group predicted.id_CZI_Human \  # falls back to predicted.id if absent
#     --outdir /.../MS_H3K27ac_bigwig_rc \
#     --workers 4 --tilesize 100 --cutoff 4 --mincells 5
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(GenomicRanges)
  library(rtracklayer); library(future); library(future.apply); library(pbapply)
})

# --- tiny arg parser (--key value) ---
a <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, default = NULL) {
  i <- which(a == paste0("--", k)); if (length(i)) a[i + 1] else default
}
rds      <- getarg("rds")
mark     <- getarg("mark")
assay    <- getarg("assay", "bins")
group    <- getarg("group", "predicted.id_CZI_Human")
outdir   <- getarg("outdir")
# NB: default sequential. plan(multisession) here broke export.bw (per-worker
# cwd/tempdir issues + the huge MOL group OOM'ing a concurrent worker); run
# groups one at a time for reliability.
workers  <- as.integer(getarg("workers", "1"))
tileSize <- as.integer(getarg("tilesize", "100"))
cutoff   <- as.integer(getarg("cutoff", "4"))
minCells <- as.integer(getarg("mincells", "5"))
stopifnot(!is.null(rds), !is.null(mark), !is.null(outdir))

source("/home/mattia/Multi_nanoCTRNA/workflow/pbmc_deeptools/scripts/ExportGroupBW.R")

cat("Loading", rds, "-> [[", mark, "]]  (large, be patient)\n")
obj <- readRDS(rds)[[mark]]
cat("cells:", ncol(obj), "| assays:", paste(Assays(obj), collapse = ", "), "\n")

# resolve group column (prefer requested; else generic predicted.id)
if (!group %in% colnames(obj@meta.data)) {
  if ("predicted.id" %in% colnames(obj@meta.data)) {
    cat(sprintf("!! '%s' not found; falling back to 'predicted.id'\n", group))
    group <- "predicted.id"
  } else stop("No usable group column ('", group, "' / 'predicted.id') in meta.data")
}
obj@meta.data[[group]] <- as.character(obj@meta.data[[group]])
cat("Grouping by:", group, "\n")
print(table(obj@meta.data[[group]], useNA = "ifany"))

DefaultAssay(obj) <- assay
if (workers > 1) plan(multisession, workers = workers)
options(future.globals.maxSize = 8 * 1024^3)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Exporting bigwigs to", outdir, "\n")
bw <- ExportGroupBW(
  object     = obj,
  assay      = assay,
  group.by   = group,
  normMethod = "RC",
  tileSize   = tileSize,
  minCells   = minCells,
  cutoff     = cutoff,
  outdir     = outdir,
  verbose    = TRUE
)
cat("Done. Wrote", length(unlist(bw)), "bigwig(s):\n")
print(list.files(outdir, pattern = "\\.bw$", full.names = TRUE))
