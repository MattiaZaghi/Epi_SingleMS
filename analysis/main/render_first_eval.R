#!/usr/bin/env Rscript
# =============================================================================
# render_first_eval.R
# -----------------------------------------------------------------------------
# Small runner that renders the parametrized R Markdown report
#   analysis/main/first_eval_MGI.Rmd
# to a SINGLE self-contained HTML file, passing params and making sure the
# high-res PNG output directory exists first.
#
# Usage:
#   Rscript analysis/main/render_first_eval.R
#   Rscript analysis/main/render_first_eval.R resolution=0.4 filter_passed=FALSE
#
# Any default param below can be overridden on the command line as key=value
# pairs (numerics/logicals are coerced automatically).
# =============================================================================

# --- Path to the Rmd we render (kept next to this runner) --------------------
rmd <- "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/analysis/main/first_eval_MGI.Rmd"

# --- Default parameters passed to the Rmd via params= ------------------------
# These mirror the `params:` block expected inside first_eval_MGI.Rmd.
default_params <- list(
  wd                = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/",
  experiment_source = "discover",  # "discover" (scan wd) or "config"
  config_path       = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/config/config_human_MGI.yaml",
  genome            = "hg38",
  binsize           = 5000,       # bin feature-space size (bp)
  filter_passed     = TRUE,       # keep only passedMB==TRUE cells
  qc_low            = 0.01,       # lower logUMI quantile for optional QC trim
  qc_high           = 0.99,       # upper logUMI quantile for optional QC trim
  dims              = 30,         # LSI dims (2:dims used downstream)
  resolution        = 0.6,        # clustering resolution
  workers           = 4,          # parallel cores for independent per-experiment work
  group_pattern     = "^(MS)",    # joint UMAP only within samples matching the same group
  dataset_label     = "this_run", # label for the current dataset in QC-benchmark plots
  reference_yaml    = "",         # optional external-datasets YAML (qc/QC.R-style benchmark)
  compute_tss       = TRUE,      # compute TSS enrichment per experiment (heavy) + add to QC
  reuse_counts      = TRUE,       # reuse cached 5kb count matrices instead of rebuilding
  fig_scale         = 1.25,       # scale factor applied to every saved PNG
  include_rna       = FALSE,      # analyse RNA experiments (mark == "RNA") if present
  rna_dims          = 15,         # RNA PCA dims for UMAP/clustering
  rna_res           = 0.8,        # RNA clustering resolution
  rna_max_count     = 25000,      # RNA QC: max nCount_RNA
  rna_min_count     = 200,        # RNA QC: min nCount_RNA
  rna_max_mt        = 5,          # RNA QC: max percent.MT
  fig_dir           = "/date/gcb/gcb_wq/nanoCTAR_pipeline_barcode_correction/first_eval_figures"
)

# --- Coerce a raw string value into a sensible R type ------------------------
# Order matters: logicals first, then numerics, else leave as a string.
coerce_value <- function(x) {
  # case-insensitive booleans -- also handles Python's "True"/"False" that the
  # Snakefile passes for config booleans (str(True) == "True").
  if (toupper(x) %in% c("TRUE", "FALSE", "T", "F")) {
    return(toupper(x) %in% c("TRUE", "T"))
  }
  # Treat as numeric only if it fully parses as a finite number.
  num <- suppressWarnings(as.numeric(x))
  if (!is.na(num) && grepl("^\\s*-?[0-9.eE+-]+\\s*$", x)) {
    return(num)
  }
  x
}

# --- Parse `key=value` command-line overrides --------------------------------
parse_overrides <- function(args) {
  overrides <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      warning("Ignoring argument without '=': ", arg)
      next
    }
    kv  <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- kv[1]
    # Re-join in case the value itself contains '=' (e.g. a path or query).
    val <- paste(kv[-1], collapse = "=")
    overrides[[key]] <- coerce_value(val)
  }
  overrides
}

# --- Build the final param list (overrides merged over defaults) -------------
cli_args   <- commandArgs(trailingOnly = TRUE)
overrides  <- parse_overrides(cli_args)

# `out_html` is a RUNNER-level override (a fixed output path), not an Rmd param.
# Snakemake/Condor pass it so the output filename is deterministic; interactive
# runs omit it and get a timestamped name instead.
out_override <- overrides[["out_html"]]
overrides[["out_html"]] <- NULL

unknown <- setdiff(names(overrides), names(default_params))
if (length(unknown) > 0) {
  warning("Unknown param(s) not in defaults, passing through anyway: ",
          paste(unknown, collapse = ", "))
}

params <- modifyList(default_params, overrides)

# --- Ensure the high-res PNG output directory exists -------------------------
dir.create(params$fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Self-contained HTML output path -----------------------------------------
# A fixed path when `out_html=` is supplied (Snakemake/Condor); otherwise a
# timestamped name next to the Rmd (interactive runs). Timestamp lives in the
# runner, never the Rmd, so the report stays reproducible.
if (!is.null(out_override)) {
  out_html <- out_override
} else {
  stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_html <- file.path(dirname(rmd), paste0("first_eval_MGI_", stamp, ".html"))
}
dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)

# --- Render ------------------------------------------------------------------
rmarkdown::render(
  input           = rmd,
  output_file     = out_html,
  output_options  = list(self_contained = TRUE),
  params          = params,
  envir           = new.env(parent = globalenv())
)

# --- Report resolved params and final HTML path ------------------------------
cat("\n=== Resolved params ===\n")
for (nm in names(params)) {
  cat(sprintf("  %-13s = %s\n", nm, format(params[[nm]])))
}
cat("\n=== Output HTML ===\n")
cat("  ", out_html, "\n", sep = "")
