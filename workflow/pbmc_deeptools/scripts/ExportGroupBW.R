#!/usr/bin/env Rscript
# ExportGroupBW + CreateBWGroup, sourced from Signac PR #1383 (not in Signac 1.17.1).
# Exports one bigwig per group.by level by splitting fragments per group and
# tiling coverage. Requires: Signac (SplitFragments, DefaultAssay, Idents),
# rtracklayer, GenomicRanges, IRanges, S4Vectors, BiocGenerics, Matrix,
# future/future.apply, pbapply. SetIfNull is Signac-internal, defined here.
SetIfNull <- function(x, y) if (is.null(x)) y else x

ExportGroupBW <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = NULL,
    verbose = TRUE
) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer."); return(NULL)
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  DefaultAssay(object = object) <- assay
  group.by <- SetIfNull(x = group.by, y = "ident")
  Idents(object = object) <- group.by
  idents <- SetIfNull(x = idents, y = levels(x = object))
  GroupsNames <- names(x = table(object[[group.by]])[table(object[[group.by]]) > minCells])
  GroupsNames <- GroupsNames[GroupsNames %in% idents]
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present and will be overwritten!", x))
      file.remove(fn)
    }
  })
  Signac::SplitFragments(
    object = object, assay = assay, group.by = group.by, idents = idents,
    outdir = outdir, file.suffix = "", append = TRUE,
    buffer_length = 256L, verbose = verbose
  )
  if (!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c("rc", "ncells", "none")) {
      normBy <- normMethod
    } else {
      normBy <- object[[normMethod, drop = FALSE]]
    }
  }
  if (!is.null(x = chromosome)) seqlevels(object) <- chromosome
  availableChr <- names(x = seqlengths(object))
  chromLengths <- seqlengths(object)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(start = rep(1, length(x = availableChr)),
                     end = as.numeric(x = chromLengths))
  )
  if (verbose) message("Creating tiles")
  tiles <- unlist(x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize))
  if (verbose) message("Creating bigwig files at ", outdir)
  if (future::nbrOfWorkers() > 1) {
    mylapply <- future.apply::future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
  }
  covFiles <- mylapply(
    GroupsNames, FUN = CreateBWGroup, availableChr, chromLengths, tiles,
    normBy, tileSize, normMethod, cutoff, outdir
  )
  return(covFiles)
}

CreateBWGroup <- function(
    groupNamei, availableChr, chromLengths, tiles, normBy,
    tileSize, normMethod, cutoff, outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer."); return(NULL)
  }
  normMethod <- tolower(x = normMethod)
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"), format = "bed")
  cellGroupi <- unique(x = fragi$name)
  covFile <- file.path(
    outdir, paste0(groupNamei, "-TileSize-", tileSize, "-normMethod-", normMethod, ".bw"))
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k], ]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
    } else {
      nTiles <- chromLengths[availableChr[k]] / tileSize
      if (nTiles %% 1 != 0) nTiles <- trunc(x = nTiles) + 1
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      mat <- Matrix::sparseMatrix(
        i = c(trunc(x = start(x = fragik) / tileSize),
              trunc(x = end(x = fragik) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2 * length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      if (!is.null(x = cutoff)) mat@x[mat@x > cutoff] <- cutoff
      mat <- Matrix::rowSums(x = mat)
      tilesk$reads <- mat
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else if (!is.null(x = normBy)) {
          tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}
