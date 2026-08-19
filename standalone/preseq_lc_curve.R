#!/usr/bin/env Rscript
# =============================================================================
# preseq_lc_curve.R
# -----------------------------------------------------------------------------
# Library-complexity (preseq `lc_extrap`) curves, folder-agnostic. Discovers
# every `.../cellranger/outs/yield.txt` under an analysis working directory,
# parses sample + modality from the path, and plots the LC curves with ONE LINE
# PER MODALITY, faceted per sample (so each sample panel shows its 2-3 marks).
#
# Complements the hand-tuned standalone/lc_curve.R (which is manual). This one
# needs no editing -- point it at a run folder and it follows whatever is there.
#
# Usage:
#   conda run -n nanoctarna-analysis Rscript standalone/preseq_lc_curve.R \
#       wd=/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS
#   ... preseq_lc_curve.R wd=/path/run fig_dir=/path/figs max_reads_M=300 out=preseq_lc.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --- parse key=value CLI args over defaults ----------------------------------
defaults <- list(
  wd          = "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS",
  fig_dir     = "",          # default: <wd>/first_eval/figures
  out         = "preseq_lc_curves.png",
  max_reads_M = "500"        # drop extrapolated points beyond this (millions)
)
args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) >= 2) defaults[[kv[1]]] <- paste(kv[-1], collapse = "=")
}
wd          <- defaults$wd
fig_dir     <- if (nzchar(defaults$fig_dir)) defaults$fig_dir else file.path(wd, "first_eval", "figures")
max_reads_M <- suppressWarnings(as.numeric(defaults$max_reads_M))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

base_mark <- function(m) sub("[-_][0-9]+$", "", m)   # H3K27ac-1 -> H3K27ac

# --- discover yield.txt under wd ---------------------------------------------
files <- list.files(wd, pattern = "^yield\\.txt$", recursive = TRUE, full.names = TRUE)
files <- files[grepl("/cellranger/outs/yield\\.txt$", files)]
if (length(files) == 0)
  stop("No cellranger/outs/yield.txt found under ", wd,
       " -- run the preseq rule first.")
cat("Found", length(files), "yield.txt file(s).\n")

wd_norm <- normalizePath(wd, mustWork = FALSE)

dat <- do.call(rbind, lapply(files, function(f) {
  df <- tryCatch(read.delim(f, header = TRUE, sep = "\t", check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  rel   <- sub(paste0("^", wd_norm, "/?"), "", normalizePath(f, mustWork = FALSE))
  parts <- strsplit(rel, "/")[[1]]
  # accept ONLY the exact layout <sample>/<modality_barcode>/cellranger/outs/yield.txt
  # so stray yield.txt nested elsewhere in the tree can't leak in
  if (length(parts) != 5 || parts[3] != "cellranger" || parts[4] != "outs") return(NULL)
  sample   <- parts[1]
  expdir   <- parts[2]
  modality <- sub("_[^_]+$", "", expdir)     # strip trailing _<barcode>
  df$sample     <- sample
  df$modality   <- modality
  df$mark       <- base_mark(modality)
  df$experiment <- paste0(sample, ".", modality)
  df$expdir     <- expdir            # <modality>_<barcode> folder, for locating FastQC
  df
}))
if (is.null(dat)) stop("Could not parse any yield.txt")

# preseq lc_extrap columns: TOTAL_READS, EXPECTED_DISTINCT (+ CI columns)
names(dat)[names(dat) == "TOTAL_READS"]       <- "total_reads"
names(dat)[names(dat) == "EXPECTED_DISTINCT"] <- "expected_distinct"
dat <- dat %>%
  mutate(total_M = total_reads / 1e6, distinct_M = expected_distinct / 1e6) %>%
  filter(total_M <= max_reads_M) %>%
  arrange(experiment, total_M)

# --- actual reads sequenced per experiment (from FastQC of debarcoded fastqs) -
# Each barcode folder has a fastqc/ dir with one *_fastqc.zip per lane per read.
# Sum "Total Sequences" over lanes for a SINGLE read type (R1/R2/R3 have equal
# counts, so picking one avoids 3x over-counting) -> read pairs sequenced.
fastqc_total_reads <- function(exp_base) {
  fqdir <- file.path(exp_base, "fastqc")
  if (!dir.exists(fqdir)) return(NA_real_)
  zips <- list.files(fqdir, pattern = "_fastqc\\.zip$", full.names = TRUE)
  if (length(zips) == 0) return(NA_real_)
  for (rd in c("_R1_", "_R3_", "_R2_")) {            # prefer one read type
    sel <- zips[grepl(rd, basename(zips))]
    if (length(sel) > 0) { zips <- sel; break }
  }
  total <- 0; got <- FALSE
  for (z in zips) {
    inner <- paste0(sub("\\.zip$", "", basename(z)), "/fastqc_data.txt")
    lines <- tryCatch(readLines(unz(z, inner), warn = FALSE), error = function(e) NULL)
    if (is.null(lines)) next
    ts <- lines[grepl("^Total Sequences", lines)]
    if (length(ts) == 0) next
    n <- suppressWarnings(as.numeric(sub(".*\\t", "", ts[1])))
    if (!is.na(n)) { total <- total + n; got <- TRUE }
  }
  if (got) total else NA_real_
}

exp_meta <- unique(dat[, c("experiment", "sample", "modality", "expdir")])
exp_meta$actual_reads <- vapply(seq_len(nrow(exp_meta)), function(i)
  fastqc_total_reads(file.path(wd, exp_meta$sample[i], exp_meta$expdir[i])), numeric(1))
exp_meta$actual_M <- exp_meta$actual_reads / 1e6
# expected-distinct at the sequenced depth, interpolated from each LC curve
exp_meta$distinct_M <- vapply(seq_len(nrow(exp_meta)), function(i) {
  if (is.na(exp_meta$actual_M[i])) return(NA_real_)
  d <- dat[dat$experiment == exp_meta$experiment[i], ]
  approx(d$total_M, d$distinct_M, xout = exp_meta$actual_M[i], rule = 2)$y
}, numeric(1))
actual_pts <- exp_meta[is.finite(exp_meta$actual_M) & is.finite(exp_meta$distinct_M), , drop = FALSE]
if (nrow(actual_pts) > 0) {
  actual_pts$label <- paste0(round(actual_pts$actual_M, 1), "M")
  cat("Sequenced depth marked for", nrow(actual_pts), "experiment(s) from FastQC.\n")
} else {
  cat("No FastQC read counts found -- drawing curves without sequenced-depth markers.\n")
}

# --- plot 1: faceted per sample, one line per modality -----------------------
p_facet <- ggplot(dat, aes(total_M, distinct_M, color = modality, group = experiment)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ sample, scales = "free") +
  labs(title = "Library complexity (preseq lc_extrap)",
       subtitle = if (nrow(actual_pts) > 0) "● = reads actually sequenced (FastQC)" else NULL,
       x = "Total reads (millions)", y = "Expected distinct reads (millions)",
       color = "Modality") +
  theme_classic() +
  theme(legend.position = "bottom")

if (nrow(actual_pts) > 0) {
  p_facet <- p_facet +
    geom_point(data = actual_pts,
               aes(x = actual_M, y = distinct_M, color = modality),
               shape = 21, fill = "white", size = 3, stroke = 1.2, inherit.aes = FALSE) +
    geom_text(data = actual_pts,
              aes(x = actual_M, y = distinct_M, label = label, color = modality),
              vjust = -1, size = 3, show.legend = FALSE, inherit.aes = FALSE)
}

out_facet <- file.path(fig_dir, defaults$out)
ggsave(out_facet, p_facet, width = 12, height = 8, dpi = 300, bg = "white")
cat("Wrote", out_facet, "(", length(unique(dat$experiment)), "curves,",
    length(unique(dat$sample)), "samples )\n")

# --- plot 2: single overlay of all experiments -------------------------------
p_all <- ggplot(dat, aes(total_M, distinct_M, color = experiment)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Library complexity - all experiments",
       subtitle = if (nrow(actual_pts) > 0) "● = reads actually sequenced (FastQC)" else NULL,
       x = "Total reads (millions)", y = "Expected distinct reads (millions)",
       color = "Experiment") +
  theme_classic()

if (nrow(actual_pts) > 0) {
  p_all <- p_all +
    geom_point(data = actual_pts,
               aes(x = actual_M, y = distinct_M, color = experiment),
               shape = 21, fill = "white", size = 2.5, stroke = 1.1, inherit.aes = FALSE)
}

out_all <- file.path(fig_dir, sub("\\.png$", "_overlay.png", defaults$out))
ggsave(out_all, p_all, width = 12, height = 8, dpi = 300, bg = "white")
cat("Wrote", out_all, "\n")
