# functions_firstEval.R
# ---------------------------------------------------------------------------
# Pure function library for a streamlined, sample/modality-AGNOSTIC, multi-sample
# single-cell nano-CUT&Tag "first evaluation" analysis (Signac/Seurat, bins
# feature space). Experiments are discovered from the analysis folder (or, as a
# fallback, a config yaml); the epigenomic path (bins) and an OPTIONAL RNA path
# (cellranger GEX) are kept separate so the same report runs with or without RNA
# and with whatever samples/modalities are present. A separate Rmd sources this
# file and calls the functions below; names / args / return types are the contract.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(yaml)
  library(ggVennDiagram)
  library(ComplexUpset)
})

# ---------------------------------------------------------------------------
# Internal helper: derive the base mark from a modality name.
# First tries to recognise a known assay/mark token anywhere in the name -- so
# descriptive names that also encode a patient/replicate (e.g. "H3K27ac_MS381",
# "H3K27ac-1") still collapse to the mark ("H3K27ac"). Falls back to stripping a
# trailing "-<n>"/"_<n>" replicate suffix for anything unrecognised.
# ---------------------------------------------------------------------------
.KNOWN_MARKS <- c("H3K27ac", "H3K27me3", "H3K4me3", "H3K4me1",
                  "H3K9me3", "H3K36me3", "H3K27me1", "ATAC", "RNA")
.base_mark <- function(modality) {
  hit <- .KNOWN_MARKS[vapply(.KNOWN_MARKS,
                             function(k) grepl(k, modality, fixed = TRUE),
                             logical(1))]
  if (length(hit) > 0) return(hit[1])
  sub("[-_][0-9]+$", "", modality)
}

#' Build the experiment table from a config yaml.
#'
#' Iterates every sample and every antibody/modality barcode declared under
#' `config$samples[[sample]]$barcodes` and returns one row per experiment.
#'
#' @param config_path path to the config yaml.
#' @return data.frame with columns sample, modality, barcode, mark, dir, exp_id.
read_experiment_table <- function(config_path) {
  config <- yaml::read_yaml(config_path)
  rows <- list()
  for (sample in names(config$samples)) {
    barcodes <- config$samples[[sample]]$barcodes
    if (is.null(barcodes)) next
    for (modality in names(barcodes)) {
      barcode <- as.character(barcodes[[modality]])
      rows[[length(rows) + 1]] <- data.frame(
        sample   = sample,
        modality = modality,
        barcode  = barcode,
        mark     = .base_mark(modality),
        dir      = paste0(modality, "_", barcode),
        exp_id   = paste0(sample, ".", modality),
        stringsAsFactors = FALSE
      )
    }
  }
  exp_df <- do.call(rbind, rows)
  rownames(exp_df) <- NULL
  exp_df
}

#' Discover experiments by scanning the analysis working directory.
#'
#' Sample/modality-agnostic: walks `<wd>/<sample>/<modality>_<barcode>/` and
#' returns every experiment folder it finds, regardless of which config produced
#' it. A subdirectory is treated as an experiment when its name parses as
#' `<modality>_<barcode>` (barcode = the final underscore-delimited token, a run
#' of A/C/G/T/N) AND it contains recognisable pipeline output. The modality is
#' everything before that final barcode token, so replicate names like
#' `H3K27ac-1` and RNA (`RNA_AAAAGGGG`) are handled. Base mark is derived the
#' same way as `read_experiment_table`.
#'
#' @param wd analysis root directory.
#' @param require what a folder must contain to count as an experiment:
#'   "any" (metadata OR fragments OR a cellranger dir; default), "metadata"
#'   (cell_picking/metadata.csv), or "fragments" (a fragments*.tsv.gz).
#' @return data.frame with columns sample, modality, barcode, mark, dir, exp_id
#'   (0 rows if nothing found).
discover_experiments <- function(wd, require = c("any", "metadata", "fragments")) {
  require <- match.arg(require)
  empty <- data.frame(sample = character(), modality = character(),
                      barcode = character(), mark = character(),
                      dir = character(), exp_id = character(),
                      stringsAsFactors = FALSE)
  if (!dir.exists(wd)) return(empty)

  samples <- list.dirs(wd, recursive = FALSE, full.names = FALSE)
  rows <- list()
  for (sample in samples) {
    subdirs <- list.dirs(file.path(wd, sample), recursive = FALSE, full.names = FALSE)
    for (d in subdirs) {
      barcode  <- sub("^.*_", "", d)         # final underscore token
      modality <- sub("_[^_]+$", "", d)      # everything before it
      if (identical(modality, d)) next       # no underscore -> not an experiment
      if (!grepl("^[ACGTN]{5,}$", barcode)) next  # final token is not a barcode

      base   <- file.path(wd, sample, d)
      has_meta <- file.exists(file.path(base, "cell_picking", "metadata.csv"))
      has_frag <- length(Sys.glob(file.path(base, "cellranger", "outs",
                                            "fragments*.tsv.gz"))) > 0
      has_cr   <- dir.exists(file.path(base, "cellranger"))
      ok <- switch(require,
                   any       = has_meta || has_frag || has_cr,
                   metadata  = has_meta,
                   fragments = has_frag)
      if (!ok) next

      rows[[length(rows) + 1]] <- data.frame(
        sample = sample, modality = modality, barcode = barcode,
        mark = .base_mark(modality), dir = d,
        exp_id = paste0(sample, ".", modality), stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) return(empty)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Build the experiment table from either the filesystem or a config yaml.
#'
#' @param source "discover" (scan `wd`, default) or "config" (read `config_path`).
#' @return experiment data.frame (see discover_experiments / read_experiment_table).
build_experiment_table <- function(wd = NULL, config_path = NULL,
                                    source = c("discover", "config")) {
  source <- match.arg(source)
  if (source == "config") {
    if (is.null(config_path)) stop("build_experiment_table: config_path required for source='config'")
    return(read_experiment_table(config_path))
  }
  if (is.null(wd)) stop("build_experiment_table: wd required for source='discover'")
  discover_experiments(wd)
}

#' Resolve the fragments file path for one experiment.
#'
#' Returns the first existing of `fragments_raw.tsv.gz` then `fragments.tsv.gz`
#' under `<wd>/<sample>/<dir>/cellranger/outs`. If neither exists, returns
#' NA_character_ (does not error).
#'
#' @return character path or NA_character_.
resolve_fragments <- function(wd, sample, dir) {
  candidates <- file.path(wd, sample, dir, "cellranger/outs",
                          c("fragments_raw.tsv.gz", "fragments.tsv.gz"))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

#' Load per-experiment metadata into a named list of data.frames.
#'
#' For each experiment in `exp_df`, reads
#' `<wd>/<sample>/<dir>/cell_picking/metadata.csv`, drops the leading unnamed
#' index column, sets rownames to `barcode`, and annotates sample / mark /
#' experiment / modality (= mark). Experiments whose metadata.csv is missing are
#' skipped with a warning. If `filter_passed` is TRUE, keeps only rows where
#' `passedMB` is TRUE.
#'
#' @return named list of data.frames keyed by exp_id (only those that loaded).
load_metadata_list <- function(exp_df, wd, filter_passed = TRUE) {
  meta.ls <- list()
  for (i in seq_len(nrow(exp_df))) {
    sample   <- exp_df$sample[i]
    dir      <- exp_df$dir[i]
    mark     <- exp_df$mark[i]
    exp_id   <- exp_df$exp_id[i]

    path <- file.path(wd, sample, dir, "cell_picking", "metadata.csv")
    if (!file.exists(path)) {
      warning(sprintf("metadata.csv not found for %s (%s) - skipping", exp_id, path))
      next
    }

    df <- read.csv(path, check.names = FALSE)
    # drop the leading unnamed row-index column
    df <- df[, -1, drop = FALSE]
    # set rownames to barcode
    rownames(df) <- df$barcode

    # optional filtering to passed cells
    if (filter_passed && "passedMB" %in% colnames(df)) {
      df <- df[which(df$passedMB %in% TRUE), , drop = FALSE]
    }

    # annotate
    df$sample     <- sample
    df$mark       <- mark
    df$experiment <- exp_id
    df$modality   <- mark

    meta.ls[[exp_id]] <- df
  }
  meta.ls
}

#' Build per-experiment bins Seurat objects.
#'
#' For each experiment present in `meta.ls`, builds a genome-bin count matrix
#' from the fragments file over the metadata barcodes and wraps it in a Seurat
#' object with a "bins" ChromatinAssay. Each experiment is wrapped in tryCatch;
#' failures warn() and are skipped.
#'
#' @param bs_genome a BSgenome object (used for seqlengths).
#' @return named list of Seurat objects keyed by exp_id.
# Internal: lapply that runs in parallel via future.apply when available (using
# whatever future::plan() the caller set), else a plain serial lapply. Used to
# fan the independent per-experiment / per-group work across cores.
.par_lapply <- function(X, FUN) {
  if (requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

build_bins_objects <- function(exp_df, meta.ls, wd, bs_genome,
                               binsize = 5000, min_cell = 1, min_feat = 1,
                               counts_cache = NULL) {
  ids <- names(meta.ls)
  build_one <- function(exp_id) {
    idx    <- match(exp_id, exp_df$exp_id)
    sample   <- exp_df$sample[idx]
    dir      <- exp_df$dir[idx]
    mark     <- exp_df$mark[idx]
    modality <- exp_df$modality[idx]
    # patient = the modality part after the mark (e.g. H3K27ac_MS381 -> MS381);
    # falls back to the sample when the modality carries no patient suffix.
    patient  <- sub(paste0("^", mark, "[-_]?"), "", modality)
    if (!nzchar(patient) || identical(patient, modality)) patient <- sample
    md       <- meta.ls[[exp_id]]

    tryCatch({
      frag_path <- resolve_fragments(wd, sample, dir)
      if (is.na(frag_path)) stop(sprintf("no fragments file found for %s", exp_id))
      cells <- rownames(md)
      fragments <- CreateFragmentObject(path = frag_path, cells = cells)
      # Reuse the cached 5kb bin matrix (the expensive GenomeBinMatrix step) when
      # it exists AND covers exactly the current cells; otherwise (re)compute it.
      counts <- NULL
      if (!is.null(counts_cache) && !is.null(counts_cache[[exp_id]])) {
        cc <- counts_cache[[exp_id]]
        if (all(cells %in% colnames(cc))) counts <- cc[, cells, drop = FALSE]
      }
      if (is.null(counts)) {
        counts <- GenomeBinMatrix(
          fragments = fragments, cells = cells, binsize = binsize,
          genome = seqlengths(bs_genome), process_n = 2e6)
      }
      chrom_assay <- CreateChromatinAssay(
        counts = counts, fragments = fragments, genome = "hg38",
        min.cells = min_cell, min.features = min_feat)
      seurat_obj <- CreateSeuratObject(
        counts = chrom_assay, assay = "bins", meta.data = md, project = sample)
      seurat_obj$sample     <- sample
      seurat_obj$mark       <- mark
      seurat_obj$experiment <- exp_id
      seurat_obj$modality   <- mark
      seurat_obj$patient    <- patient
      seurat_obj
    }, error = function(e) {
      warning(sprintf("build_bins_objects failed for %s: %s", exp_id, conditionMessage(e)))
      NULL
    })
  }

  # experiments are independent -> fan out across cores
  objs <- .par_lapply(ids, build_one)
  names(objs) <- ids
  objs[!vapply(objs, is.null, logical(1))]
}

#' Extract the raw bin count matrices (sparse) from a list of bins objects.
#'
#' Returns a named list of dgCMatrix (features x cells) keyed by exp_id, so the
#' expensive GenomeBinMatrix step can be reused without rebuilding. Robust to the
#' Seurat v5 layer / legacy slot API.
#'
#' @return named list of sparse count matrices (NULL entries dropped).
extract_counts_list <- function(obj.ls, assay = "bins") {
  cl <- lapply(obj.ls, function(o) {
    tryCatch(
      SeuratObject::GetAssayData(o, assay = assay, layer = "counts"),
      error = function(e) tryCatch(
        SeuratObject::GetAssayData(o, assay = assay, slot = "counts"),
        error = function(e2) NULL))
  })
  cl[!vapply(cl, is.null, logical(1))]
}

#' Merge per-experiment objects into one object per base mark.
#'
#' Groups the objects in `obj.ls` by their base mark (via `exp_df`) and
#' merges each group with `add.cell.ids` set to the contributing exp_ids. Marks
#' with a single object are returned unchanged.
#'
#' @return named list of Seurat objects keyed by unique mark.
merge_by_mark <- function(obj.ls, exp_df) {
  merged.ls <- list()
  # only consider exp_ids that actually produced an object
  present <- exp_df[exp_df$exp_id %in% names(obj.ls), , drop = FALSE]

  for (mark in unique(present$mark)) {
    ids <- present$exp_id[present$mark == mark]
    objs <- obj.ls[ids]

    if (length(objs) == 1) {
      merged.ls[[mark]] <- objs[[1]]
    } else {
      merged.ls[[mark]] <- merge(
        x            = objs[[1]],
        y            = objs[-1],
        add.cell.ids = ids
      )
    }
  }
  merged.ls
}

#' Assign a group label to each sample (for joint-within-type UMAPs).
#'
#' Applies `pattern` (a regex) to each sample name; the matched substring is the
#' group (e.g. "MS"), so replicates/conditions of the same type share a
#' group while different types stay separate. Samples that do not match become
#' their own group (i.e. no forced joining). The default groups by a leading
#' "MS" token.
#'
#' @return a character vector of group labels, parallel to `samples`.
sample_group <- function(samples, pattern = "^(MS)") {
  vapply(samples, function(s) {
    m <- regmatches(s, regexpr(pattern, s, perl = TRUE))
    if (length(m) == 0 || !nzchar(m[1])) s else m[1]
  }, character(1), USE.NAMES = FALSE)
}

#' Merge per-experiment objects by (group x base mark).
#'
#' Like `merge_by_mark`, but only pools objects that share BOTH a sample group
#' (via `sample_group`) AND a base mark -- so a joint UMAP is built per type
#' (e.g. MS H3K27ac) and never mixes unrelated sample types.
#' Keys are "<group>|<mark>".
#'
#' @return named list of Seurat objects keyed by "<group>|<mark>".
merge_by_group_mark <- function(obj.ls, exp_df, group_pattern = "^(MS)") {
  present <- exp_df[exp_df$exp_id %in% names(obj.ls), , drop = FALSE]
  present$group <- sample_group(present$sample, group_pattern)
  present$gm    <- paste(present$group, present$mark, sep = "|")

  merged.ls <- list()
  for (key in unique(present$gm)) {
    ids  <- present$exp_id[present$gm == key]
    objs <- obj.ls[ids]
    merged.ls[[key]] <- if (length(objs) == 1) objs[[1]]
                        else merge(objs[[1]], y = objs[-1], add.cell.ids = ids)
  }
  merged.ls
}

#' Normalize and reduce dimensionality of a bins object.
#'
#' Runs TF-IDF, top-feature selection, and SVD (LSI).
#'
#' @return the input object with TFIDF/top-features/SVD applied.
normalize_dr <- function(obj, min_cutoff = 5) {
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = min_cutoff)
  obj <- RunSVD(obj)
  obj
}

#' UMAP + graph-based clustering on the LSI reduction.
#'
#' @return the input object with UMAP, neighbours, and clusters computed.
run_umap_cluster <- function(obj, dims = 2:30, res = 0.6, algorithm = 3) {
  obj <- RunUMAP(obj, reduction = "lsi", dims = dims)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, resolution = res, algorithm = algorithm, verbose = FALSE)
  obj
}

# ---------------------------------------------------------------------------
# Internal helper: bind a list of metadata data.frames into a single tidy
# frame keeping only the columns needed for plotting plus the requested
# feature. Robust to differing columns across experiments.
# ---------------------------------------------------------------------------
.bind_meta <- function(meta.ls, extra_cols = character(0)) {
  keep <- unique(c("sample", "mark", "experiment", extra_cols))
  parts <- list()
  for (exp_id in names(meta.ls)) {
    df <- meta.ls[[exp_id]]
    # guarantee annotation columns
    if (!"experiment" %in% colnames(df)) df$experiment <- exp_id
    have <- intersect(keep, colnames(df))
    sub <- df[, have, drop = FALSE]
    # add any missing requested columns as NA so rbind works
    for (col in setdiff(keep, have)) sub[[col]] <- NA
    parts[[exp_id]] <- sub[, keep, drop = FALSE]
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Violin (+ inner boxplot) of a per-cell QC feature across experiments.
#'
#' x = experiment, fill = sample, faceted by mark. Optionally log-transforms y.
#'
#' @return a ggplot object.
qc_violin <- function(meta.ls, feature, ylab = feature, log_y = FALSE) {
  df <- .bind_meta(meta.ls, extra_cols = feature)
  df$experiment <- factor(df$experiment)
  df$sample     <- factor(df$sample)

  p <- ggplot(df, aes(x = experiment, y = .data[[feature]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    facet_wrap(~mark, scales = "free_x") +
    ylab(ylab) +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (log_y) {
    p <- p + scale_y_log10()
  }
  p
}

#' Barplot of the number of cells per experiment.
#'
#' x = experiment, fill = mark, faceted by sample.
#'
#' @return a ggplot object.
ncells_barplot <- function(meta.ls) {
  rows <- list()
  for (exp_id in names(meta.ls)) {
    df <- meta.ls[[exp_id]]
    rows[[exp_id]] <- data.frame(
      experiment = exp_id,
      sample     = if ("sample" %in% colnames(df)) df$sample[1] else NA_character_,
      mark       = if ("mark" %in% colnames(df)) df$mark[1] else NA_character_,
      ncells     = nrow(df),
      stringsAsFactors = FALSE
    )
  }
  counts <- do.call(rbind, rows)
  counts$experiment <- factor(counts$experiment)

  ggplot(counts, aes(x = experiment, y = ncells, fill = mark)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    facet_wrap(~sample, scales = "free_x") +
    ylab("Number of cells") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Overlap of shared cell barcodes across marks within one sample.
#'
#' Generic replacement for `commonCellHistonMarks`. Takes the barcode set of
#' every experiment belonging to `sample_name` (optionally restricted to
#' `marks`), and renders a Venn diagram (2-4 sets, via ggVennDiagram) or an
#' UpSet plot (>4 sets, via ComplexUpset) of the barcode overlaps. Set names are
#' the experiment modality names.
#'
#' @return a ggplot object.
overlap_plot <- function(meta.ls, sample_name, marks = NULL) {
  # experiments belonging to this sample
  sets <- list()
  for (exp_id in names(meta.ls)) {
    df <- meta.ls[[exp_id]]
    smp <- if ("sample" %in% colnames(df)) df$sample[1] else NA_character_
    if (is.na(smp) || smp != sample_name) next
    mk <- if ("mark" %in% colnames(df)) df$mark[1] else NA_character_
    if (!is.null(marks) && !(mk %in% marks)) next
    # set name = experiment modality name
    modality <- if ("modality" %in% colnames(df)) df$modality[1] else mk
    set_name <- sub(paste0("^", sample_name, "\\."), "", exp_id)
    sets[[set_name]] <- rownames(df)
  }

  if (length(sets) == 0) {
    stop(sprintf("overlap_plot: no experiments found for sample '%s'", sample_name))
  }

  if (length(sets) <= 4) {
    p <- ggVennDiagram::ggVennDiagram(sets) +
      scale_fill_gradient(low = "white", high = "white") +
      theme(legend.position = "none") +
      ggtitle(sample_name) +
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  } else {
    # build membership data.frame for ComplexUpset
    all_bc <- unique(unlist(sets, use.names = FALSE))
    membership <- data.frame(barcode = all_bc, stringsAsFactors = FALSE)
    for (set_name in names(sets)) {
      membership[[set_name]] <- membership$barcode %in% sets[[set_name]]
    }
    p <- ComplexUpset::upset(
      membership,
      intersect   = names(sets),
      min_size    = 1,
      width_ratio = 0.2
    ) + ggtitle(sample_name)
  }
  p
}

#' Save a ggplot to PNG, creating the target directory if needed.
#'
#' @return (invisibly) the path written.
save_png <- function(plot, name, dir, width = 8, height = 6, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".png"))
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

# ---------------------------------------------------------------------------
# OPTIONAL RNA MODALITY (only used when the report is run with include_rna=TRUE)
# RNA experiments are those whose base mark is "RNA". They are NOT part of the
# fragments/bins path -- they are loaded from the cellranger GEX
# `filtered_feature_bc_matrix.h5` and analysed with the standard
# log-normalize -> PCA -> UMAP/cluster workflow. This keeps the epigenomic and
# RNA paths cleanly separate so the same report runs with or without RNA.
# ---------------------------------------------------------------------------

#' Split an experiment table into epigenomic and RNA subsets.
#'
#' @return list(epi = <exp_df rows with mark != "RNA">, rna = <rows mark == "RNA">).
split_experiments <- function(exp_df) {
  is_rna <- exp_df$mark == "RNA"
  list(epi = exp_df[!is_rna, , drop = FALSE],
       rna = exp_df[ is_rna, , drop = FALSE])
}

#' Load per-experiment RNA Seurat objects from cellranger filtered h5.
#'
#' Only processes rows of `exp_df` whose mark == "RNA". Reads
#' `<wd>/<sample>/<dir>/cellranger/outs/filtered_feature_bc_matrix.h5`, creates
#' an RNA assay, computes percent.MT, and annotates sample/mark/experiment/
#' modality. Missing or failed experiments warn() and are skipped.
#'
#' @return named list of Seurat objects keyed by exp_id.
load_rna_objects <- function(exp_df, wd, min_cells = 0, min_features = 0,
                             mt_pattern = "^MT-") {
  rna_df <- exp_df[exp_df$mark == "RNA", , drop = FALSE]
  obj.ls <- list()
  for (i in seq_len(nrow(rna_df))) {
    sample <- rna_df$sample[i]; dir <- rna_df$dir[i]; exp_id <- rna_df$exp_id[i]
    path <- file.path(wd, sample, dir, "cellranger/outs",
                      "filtered_feature_bc_matrix.h5")
    obj <- tryCatch({
      if (!file.exists(path)) stop(sprintf("h5 not found: %s", path))
      counts <- Read10X_h5(path)
      o <- CreateSeuratObject(counts = counts, assay = "RNA", project = sample,
                              min.cells = min_cells, min.features = min_features)
      o[["percent.MT"]] <- PercentageFeatureSet(o, pattern = mt_pattern)
      o$sample <- sample; o$mark <- "RNA"; o$experiment <- exp_id; o$modality <- "RNA"
      o
    }, error = function(e) {
      warning(sprintf("load_rna_objects failed for %s: %s", exp_id, conditionMessage(e)))
      NULL
    })
    if (!is.null(obj)) obj.ls[[exp_id]] <- obj
  }
  obj.ls
}

#' Filter RNA cells by count / mitochondrial thresholds.
#'
#' @return the subset object, or the input unchanged (with a warning) if no
#'   cells pass.
qc_filter_rna <- function(obj, max_count = 25000, min_count = 200, max_mt = 5) {
  keep <- obj$nCount_RNA > min_count & obj$nCount_RNA < max_count &
          obj$percent.MT < max_mt
  if (sum(keep) == 0) {
    warning("qc_filter_rna: no cells pass thresholds; returning input unchanged")
    return(obj)
  }
  subset(obj, cells = colnames(obj)[keep])
}

#' Standard RNA normalization + PCA (LogNormalize workflow).
#'
#' @return the input object with normalization, variable features, scaling, PCA.
normalize_rna <- function(obj, nfeatures = 2000) {
  # merging in Seurat v5 leaves split count layers; join them before the
  # single-matrix LogNormalize -> PCA workflow (we merge, not integrate).
  obj <- tryCatch(SeuratObject::JoinLayers(obj), error = function(e) obj)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj
}

#' RNA UMAP + graph-based clustering on the PCA reduction.
#'
#' @return the input object with neighbours, clusters and UMAP computed.
run_umap_cluster_rna <- function(obj, dims = 1:15, res = 0.8) {
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims)
  obj
}

# ---------------------------------------------------------------------------
# CROSS-DATASET QC BENCHMARK (reproduces the per-cell comparison panels of
# qc/QC.R: fragments/cell, fragment/feature, FRiP, reads/cell, duplication
# rate -- for the current run and, optionally, external reference datasets).
# ---------------------------------------------------------------------------

# internal: first colname of `md` matching any regex in `patterns`, else NA.
.first_col <- function(md, patterns) {
  cn <- colnames(md)
  for (p in patterns) {
    hit <- grep(p, cn, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

#' Extract a standardized per-cell QC table from a Seurat object / metadata df.
#'
#' Auto-detects assay-specific column names (`nCount_*bins` / `nFeature_*bins`)
#' and falls back to metadata.csv fields, so it works across datasets that name
#' columns differently (nanoCT, 10x multiome, Droplet Paired-Tag, scCUT&Tag-pro).
#'
#' @param x a Seurat object or a per-cell metadata data.frame.
#' @param dataset dataset label.
#' @param mark base mark label; if NULL, taken from a `mark` column if present.
#' @return data.frame(dataset, mark, fragments, features, frip, reads, duplicate,
#'   dup_rate, frag_per_feat) with one row per cell, or NULL if empty.
collect_qc_metrics <- function(x, dataset, mark = NULL, sample = NA_character_) {
  md <- if (inherits(x, "Seurat")) x@meta.data else as.data.frame(x)
  if (nrow(md) == 0) return(NULL)
  if (is.null(mark))
    mark <- if ("mark" %in% colnames(md)) as.character(md$mark) else NA_character_

  frag_col <- .first_col(md, c("^nCount_.*bins$", "^nCount_bins$", "^all_unique_MB$"))
  feat_col <- .first_col(md, c("^nFeature_.*bins$", "^nFeature_bins$"))
  frip_col <- .first_col(md, c("^peak_ratio_MB$", "^peak_region_ratio$", "^FRiP$"))
  read_col <- .first_col(md, c("^total$", "^passed_filters$"))
  pf_col   <- .first_col(md, c("^passed_filters$"))
  dup_col  <- .first_col(md, c("^duplicate$"))

  num <- function(col) if (!is.na(col) && col %in% colnames(md)) as.numeric(md[[col]]) else NA_real_
  fragments <- num(frag_col); features <- num(feat_col); frip <- num(frip_col)
  reads <- num(read_col); duplicate <- num(dup_col)

  data.frame(
    dataset   = dataset,
    sample    = sample,
    mark      = mark,
    fragments = fragments,
    features  = features,
    frip      = frip,
    reads     = reads,
    passed_filters = num(pf_col),      # cross-dataset-comparable unique fragments
    duplicate = duplicate,
    dup_rate  = duplicate / reads,
    frag_per_feat = fragments / features,
    stringsAsFactors = FALSE
  )
}

#' Load reference QC metrics from a YAML describing external `.rds` datasets.
#'
#' YAML format:
#'   references:
#'     - label: nanoCT
#'       path: /path/combined.obj.ls_all_mod_nanoCT.rds
#'       marks: [ATAC_TATAGCCT, H3K27ac_ATAGAGGC]   # optional subset of list names
#' Each `.rds` must be a named list of Seurat objects (names like
#' `<mark>_<barcode>`); the base mark is taken from each element's first token.
#'
#' @return a long data.frame, or NULL if the yaml is missing/empty/unreadable.
load_reference_qc <- function(reference_yaml) {
  if (is.null(reference_yaml) || !nzchar(reference_yaml) ||
      !file.exists(reference_yaml)) return(NULL)
  spec <- yaml::read_yaml(reference_yaml)
  refs <- spec$references
  if (is.null(refs) || length(refs) == 0) return(NULL)

  parts <- list()
  for (r in refs) {
    label <- r$label; path <- r$path
    obj_list <- tryCatch(readRDS(path), error = function(e) {
      warning(sprintf("load_reference_qc: cannot read %s: %s", path, conditionMessage(e)))
      NULL })
    if (is.null(obj_list)) next
    elem_names <- names(obj_list)
    if (!is.null(r$marks)) elem_names <- intersect(elem_names, r$marks)
    for (nm in elem_names) {
      mk <- .base_mark(nm)   # token-match the mark anywhere in the element name
                             # (handles "H3K27ac_CCTATCCT" AND "MS1_H3K27ac")
      m  <- tryCatch(collect_qc_metrics(obj_list[[nm]], dataset = label, mark = mk),
                     error = function(e) { warning(conditionMessage(e)); NULL })
      if (!is.null(m)) parts[[paste(label, nm, sep = ".")]] <- m
    }
  }
  if (length(parts) == 0) return(NULL)
  do.call(rbind, parts)
}

#' Combined QC table for the current run (+ optional reference datasets).
#'
#' Current-run metrics come from `obj.ls` (built bins objects carry both the
#' metadata.csv fields and nCount_bins/nFeature_bins). References are appended
#' from `reference_yaml` for cross-technology comparison like qc/QC.R.
#'
#' @return a long data.frame (rows = cells across all datasets/marks).
qc_benchmark_table <- function(obj.ls, exp_df, dataset_label = "this_run",
                               reference_yaml = NULL, group_pattern = "^(MS)") {
  parts <- list()
  for (id in names(obj.ls)) {
    idx <- match(id, exp_df$exp_id)
    parts[[id]] <- collect_qc_metrics(obj.ls[[id]], dataset = dataset_label,
                                      mark = exp_df$mark[idx], sample = exp_df$sample[idx])
  }
  cur <- do.call(rbind, parts)
  if (!is.null(cur)) cur$group <- sample_group(cur$sample, group_pattern)
  ref <- load_reference_qc(reference_yaml)
  if (!is.null(ref)) ref$group <- ref$dataset   # each reference dataset is its own group
  out <- rbind(cur, ref)
  rownames(out) <- NULL
  out
}

#' QC.R-style comparison violin for one metric, across datasets, faceted by mark.
#'
#' @param df output of qc_benchmark_table.
#' @param metric one of fragments/features/frip/reads/dup_rate/frag_per_feat.
#' @return a ggplot object.
qc_benchmark_violin <- function(df, metric, ylab = metric, log_y = FALSE, group = NULL) {
  d <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (!is.null(group)) d <- d[d$group %in% group, , drop = FALSE]
  if (nrow(d) == 0) stop(sprintf("qc_benchmark_violin: no finite values for '%s'", metric))
  # per-group view compares that group's samples; combined view shows each
  # current-run sample individually AND each reference dataset as its own violin
  # (so it never collapses the whole run into a single box).
  if (!is.null(group)) {
    xcol <- "sample"
  } else {
    d$series <- ifelse(!is.na(d$sample) & nzchar(as.character(d$sample)),
                       as.character(d$sample), as.character(d$dataset))
    xcol <- "series"
  }
  d[[xcol]] <- factor(d[[xcol]])
  p <- ggplot(d, aes(x = .data[[xcol]], y = .data[[metric]], fill = .data[[xcol]])) +
    geom_violin(scale = "width", trim = TRUE, color = "black", alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    stat_summary(fun = median, geom = "text",
                 aes(label = round(after_stat(y), 2)), vjust = -0.6, size = 3) +
    facet_wrap(~mark, scales = "free_x") +
    ylab(ylab) + xlab("") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(group)) p <- p + ggtitle(paste(group, collapse = "/"))
  if (log_y) p <- p + scale_y_log10(labels = scales::comma)
  p
}

# ---------------------------------------------------------------------------
# METADATA-BASED per-cell QC (no Seurat objects needed) -- the "beautiful
# violins" (logUMI, unique fragments, reads/cell, duplication rate, FRiP) read
# straight from metadata.csv, plus a per-experiment summary and an Excel export.
# Every plot can be drawn combined (all samples) or restricted to one group.
# ---------------------------------------------------------------------------

#' Per-cell QC table straight from the loaded metadata (no objects required).
#'
#' @return data.frame: experiment, sample, mark, group + per-cell logUMI,
#'   uniqFrag (all_unique_MB), reads (total), duplicate, dup_rate, FRiP
#'   (peak_ratio_MB). NULL if nothing loaded.
qc_cell_table <- function(meta.ls, group_pattern = "^(MS)") {
  col <- function(md, nm) if (nm %in% names(md)) as.numeric(md[[nm]]) else NA_real_
  # fragments-per-cell counted in the bins assay (QC.R's nCount_bins); present
  # only when the table is built from Seurat objects' @meta.data, NA otherwise.
  bincol <- function(md) {
    hit <- grep("^nCount_.*bins$|^nCount_bins$", names(md), value = TRUE)
    if (length(hit)) as.numeric(md[[hit[1]]]) else NA_real_
  }
  parts <- list()
  for (id in names(meta.ls)) {
    md <- meta.ls[[id]]
    if (nrow(md) == 0) next
    reads <- col(md, "total"); dup <- col(md, "duplicate")
    parts[[id]] <- data.frame(
      experiment  = id,
      sample      = if ("sample" %in% names(md)) md$sample else NA_character_,
      mark        = if ("mark" %in% names(md)) md$mark else NA_character_,
      logUMI         = col(md, "logUMI"),
      nCount_bins    = bincol(md),                     # unique fragments in 5kb bins (from objects)
      passed_filters = col(md, "passed_filters"),      # cellranger unique usable fragments
      reads          = reads,
      duplicate      = dup,
      dup_rate       = dup / reads,
      FRiP           = col(md, "peak_ratio_MB"),
      TSS_enrichment = col(md, "TSS.enrichment"),      # present only when compute_tss = TRUE
      passedMB       = if ("passedMB" %in% names(md)) as.logical(md$passedMB) else NA,
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, parts)
  if (is.null(df)) return(df)
  df$group <- sample_group(df$sample, group_pattern)
  rownames(df) <- NULL
  df
}

#' Per-experiment QC summary (means + medians of every metric) for export.
#'
#' @return data.frame, one row per experiment.
qc_summary_table <- function(cell_df) {
  if (is.null(cell_df) || nrow(cell_df) == 0) return(NULL)
  metrics <- intersect(c("logUMI", "nCount_bins", "passed_filters", "reads",
                         "duplicate", "dup_rate", "FRiP", "TSS_enrichment"),
                       names(cell_df))
  agg <- lapply(split(cell_df, cell_df$experiment), function(d) {
    base <- data.frame(experiment = d$experiment[1], sample = d$sample[1],
                       mark = d$mark[1], group = d$group[1], cells = nrow(d),
                       stringsAsFactors = FALSE)
    for (m in metrics) {
      base[[paste0(m, "_mean")]]   <- round(mean(d[[m]], na.rm = TRUE), 4)
      base[[paste0(m, "_median")]] <- round(median(d[[m]], na.rm = TRUE), 4)
    }
    base
  })
  out <- do.call(rbind, agg); rownames(out) <- NULL
  out[order(out$group, out$mark, out$sample), ]
}

#' Violin (+ box, median label) of a metadata QC metric; optionally one group.
#'
#' @return a ggplot object.
qc_metric_violin <- function(cell_df, metric, ylab = metric, log_y = FALSE, group = NULL) {
  d <- cell_df[is.finite(cell_df[[metric]]), , drop = FALSE]
  if (!is.null(group)) d <- d[d$group %in% group, , drop = FALSE]
  if (nrow(d) == 0) stop(sprintf("qc_metric_violin: no data for '%s'", metric))
  d$experiment <- factor(d$experiment); d$sample <- factor(d$sample)
  p <- ggplot(d, aes(x = experiment, y = .data[[metric]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", alpha = 0.85) +
    stat_summary(fun = median, geom = "text",
                 aes(label = round(after_stat(y), 2)), vjust = -0.6, size = 3) +
    facet_wrap(~mark, scales = "free_x") +
    ylab(ylab) + xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(group)) p <- p + ggtitle(paste(group, collapse = "/"))
  if (log_y) p <- p + scale_y_log10(labels = scales::comma)
  p
}

#' Cell-selection scatter (the pick_cells plot): log10 unique reads (logUMI, the
#' metric the GMM cell caller uses) vs FRiP, each barcode coloured by whether it
#' passed (passedMB).
#'
#' Needs an UNFILTERED per-cell table (both passed and dropped cells), i.e. built
#' from metadata loaded with filter_passed = FALSE. One panel per experiment.
#'
#' @param cell_all_df per-cell table with logUMI, FRiP, passedMB, experiment.
#' @param group optional group label to restrict to.
#' @return a ggplot object.
plot_cell_selection <- function(cell_all_df, group = NULL) {
  d <- cell_all_df
  if (!is.null(group)) d <- d[d$group %in% group, , drop = FALSE]
  d <- d[is.finite(d$logUMI) & is.finite(d$FRiP) & !is.na(d$passedMB), , drop = FALSE]
  if (nrow(d) == 0) stop("plot_cell_selection: no cells with logUMI/FRiP/passedMB")
  d$passedMB <- factor(ifelse(d$passedMB, "passedMB TRUE", "passedMB FALSE"),
                       levels = c("passedMB TRUE", "passedMB FALSE"))
  p <- ggplot(d, aes(x = logUMI, y = FRiP, color = passedMB)) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_color_manual(values = c("passedMB TRUE" = "#1F78B4",
                                  "passedMB FALSE" = "grey70")) +
    facet_wrap(~experiment, scales = "free") +
    xlab("log10 unique reads per cell (logUMI)") + ylab("FRiP (peak_ratio_MB)") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_classic() +
    theme(legend.position = "bottom")
  if (!is.null(group)) p <- p + ggtitle(paste(group, collapse = "/"))
  p
}

#' Write the QC summary to an .xlsx: an "All" sheet plus one sheet per group.
#'
#' @return (invisibly) the path written.
write_qc_excel <- function(summary_df, path) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("write_qc_excel needs the openxlsx package")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "All")
  openxlsx::writeData(wb, "All", summary_df)
  for (g in unique(summary_df$group)) {
    sh <- substr(gsub("[^A-Za-z0-9]+", "_", g), 1, 28)
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, summary_df[summary_df$group == g, , drop = FALSE])
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

#' Compute TSS enrichment on a bins object using an EnsDb annotation.
#'
#' @return the object with TSS.enrichment metadata, or unchanged on error.
add_tss_enrichment <- function(obj, ensdb, genome_name = "hg38") {
  tryCatch({
    annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb)
    GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"
    GenomeInfoDb::genome(annotations) <- genome_name
    Signac::Annotation(obj) <- annotations
    Signac::TSSEnrichment(obj, fast = FALSE)
  }, error = function(e) {
    warning("add_tss_enrichment: ", conditionMessage(e)); obj
  })
}
