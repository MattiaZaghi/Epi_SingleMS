# =============================================================================
# WHAT YOU NEED TO LOAD TO USE THIS
# =============================================================================
#
# 0) PACKAGES
#      install.packages(c("shiny", "miniUI", "plotly"))   # required
#      install.packages("gatepoints")                     # only for gateCellsBase()
#
# 1) TRY IT WITH NO DATA AT ALL  <-- start here
#      source("scripts/gateCells.R")
#      gated <- gateCellsDemo()          # synthetic FRiP vs UMI, opens the gadget
#      sapply(gated, function(d) sum(d$passedManual))
#
#    gateCellsDemo() builds a fake `metadata.ls` with the same column names and
#    the same two-cloud shape as the real data (a high-UMI/high-FRiP cell cloud
#    and a low-UMI/low-FRiP debris cloud), so you can check the interaction
#    works before pointing it at anything real.
#
# 2) THEN WITH YOUR REAL DATA
#    The only object you need is `metadata.ls` exactly as built in section 1.2
#    of scmultiNanoCT_bins.Rmd:
#
#      metadata.ls <- list()
#      for (smpl in samples) {
#        for (mod in modalities_epi) {
#          metadata.ls[[paste0(mod, "_", smpl)]] <-
#            read.csv(paste0(smpl, "/", mod, "/cell_picking/metadata.csv"),
#                     stringsAsFactors = FALSE)
#          rownames(metadata.ls[[paste0(mod, "_", smpl)]]) <-
#            metadata.ls[[paste0(mod, "_", smpl)]]$barcode
#        }
#      }
#
#      checkGateInput(metadata.ls)               # verify before launching
#      metadata.ls <- gateCells(metadata.ls, ref_col = "passedMB")
#
#    REQUIRED columns in each data.frame (these come straight out of
#    pick_cells.R / metadata.csv, so nothing extra to compute):
#      all_unique_MB   numeric  unique reads per barcode      -> x axis (log10)
#      peak_ratio_MB   numeric  fraction of reads in peaks    -> y axis (FRiP)
#      rownames        unique cell id (the vignette sets these to `barcode`)
#
#    OPTIONAL:
#      passedMB        logical  the Mclust call, drawn as point colour so you
#                               can see the automatic gate while drawing yours
#      passedMZ        logical  the manual logUMI/FRiP cutoff, same purpose
#
#    A Seurat object works too -- pass its metadata:
#      md <- gateCells(obj@meta.data, ref_col = "passedMB")
#      obj$passedManual <- md$passedManual
#
# =============================================================================


#' Check that metadata is ready for gateCells()
#'
#' Reports, per panel, whether the required columns exist, whether ids are
#' unique, and how many cells would be dropped by the log10 x axis.
#'
#' @return Invisibly, TRUE if every panel is usable.
checkGateInput <- function(mdata,
                           x_col  = "all_unique_MB",
                           y_col  = "peak_ratio_MB",
                           id_col = NULL,
                           log_x  = TRUE) {
  panels <- if (is.data.frame(mdata)) list(data = mdata) else mdata
  ok <- TRUE
  for (nm in names(panels)) {
    d <- panels[[nm]]
    msgs <- character(0)
    for (cl in c(x_col, y_col)) {
      if (!cl %in% colnames(d)) {
        msgs <- c(msgs, paste0("MISSING column '", cl, "'")); ok <- FALSE
      } else if (!is.numeric(d[[cl]])) {
        msgs <- c(msgs, paste0("'", cl, "' is not numeric")); ok <- FALSE
      }
    }
    ids <- if (is.null(id_col)) rownames(d) else as.character(d[[id_col]])
    if (is.null(ids) || anyDuplicated(ids)) {
      msgs <- c(msgs, "cell ids missing or not unique"); ok <- FALSE
    }
    if (x_col %in% colnames(d) && y_col %in% colnames(d)) {
      bad <- !is.finite(d[[x_col]]) | !is.finite(d[[y_col]]) |
        (log_x & d[[x_col]] <= 0)
      if (any(bad)) msgs <- c(msgs, sprintf("%d cell(s) not plottable (NA/<=0 on log x)",
                                            sum(bad)))
    }
    ref <- intersect(c("passedMB", "passedMZ"), colnames(d))
    cat(sprintf("%-26s %6d cells | ref cols: %-18s | %s\n",
                nm, nrow(d),
                if (length(ref)) paste(ref, collapse = ",") else "none",
                if (length(msgs)) paste(msgs, collapse = "; ") else "OK"))
  }
  invisible(ok)
}


#' Build a synthetic metadata.ls for testing gateCells()
#'
#' Same column names and same two-cloud geometry as the real cell_picking
#' metadata, so the gadget can be exercised without the HPC data mounted.
#' `passedMB` is a deliberately imperfect automatic call, so there is something
#' for a manual gate to disagree with.
makeDemoMetadata <- function(n_cells    = 3000,
                             samples    = c("MS3BL1_CT", "MS3BL2_CT"),
                             modalities = c("ATAC", "H3K27ac", "H3K27me3"),
                             seed       = 1) {
  set.seed(seed)
  out <- list()
  for (smpl in samples) {
    for (mod in modalities) {
      n_real <- round(n_cells * 0.45)
      n_deb  <- n_cells - n_real
      
      # real cells: high unique reads, high FRiP
      umi_r  <- 10^stats::rnorm(n_real, 3.6, 0.25)
      frip_r <- pmin(pmax(stats::rnorm(n_real, 0.45, 0.10), 0.02), 0.95)
      # ambient / debris: low unique reads, low FRiP
      umi_d  <- 10^stats::rnorm(n_deb, 2.3, 0.35)
      frip_d <- pmin(pmax(stats::rnorm(n_deb, 0.12, 0.06), 0.001), 0.90)
      
      umi  <- c(umi_r, umi_d)
      frip <- c(frip_r, frip_d)
      truth <- rep(c(TRUE, FALSE), c(n_real, n_deb))
      
      # imperfect automatic call: a hard threshold, which mislabels the
      # cells sitting between the two clouds
      passedMB <- umi > 10^3.1 & frip > 0.25
      
      bc <- sprintf("%s-1", vapply(seq_len(n_cells), function(i)
        paste0(sample(c("A", "C", "G", "T"), 16, replace = TRUE), collapse = ""),
        character(1)))
      bc <- make.unique(bc, sep = "_")
      
      d <- data.frame(
        barcode       = bc,
        all_unique_MB = umi,
        peak_MB       = round(umi * frip),
        peak_ratio_MB = frip,
        logUMI        = log10(umi),
        passedMB      = passedMB,
        passedMZ      = log10(umi) > 0 & frip > 0.1,
        true_cell     = truth,              # ground truth, demo only
        sample        = smpl,
        modality      = mod,
        stringsAsFactors = FALSE
      )
      rownames(d) <- d$barcode
      out[[paste0(mod, "_", smpl)]] <- d[sample.int(n_cells), , drop = FALSE]
    }
  }
  out
}


#' Launch gateCells() on synthetic data -- the fastest way to test it
#'
#' @examples
#' \dontrun{
#' source("scripts/gateCells.R")
#' gated <- gateCellsDemo()
#' # how well did your hand-drawn gate recover the ground truth?
#' with(gated[[1]], table(true_cell, passedManual))
#' }
gateCellsDemo <- function(...) {
  md <- makeDemoMetadata(...)
  cat("Synthetic metadata.ls:\n")
  checkGateInput(md)
  cat("\nDraw a lasso around the upper-right cloud (high UMI, high FRiP).\n",
      "Then check your gate against the ground truth, e.g. if you assigned the\n",
      "result to `gated`:  with(gated[[1]], table(true_cell, passedManual))\n\n")
  gateCells(md, ref_col = "passedMB")
}


#' Interactive manual gating of cells on a FRiP vs UMI scatter
#'
#' Reproduces the axes of `plotPassedCells()` (x = all_unique_MB on a log10
#' scale, y = peak_ratio_MB) as an interactive plotly scatter inside an RStudio
#' gadget, and lets you draw a lasso / box around the cells you want to keep.
#' Returns the same object with a logical `passedManual` column added, so it
#' drops straight into the place of `passedMB` / `passedMZ`.
#'
#' Works on a single metadata data.frame or on a named list of them (e.g. the
#' `metadata.ls` built in scmultiNanoCT_bins.Rmd). With a list you get a
#' drop-down to switch between sample/modality panels; each panel keeps its own
#' selection, so you can gate them one at a time and collect everything in one
#' pass.
#'
#' Requires: shiny, miniUI, plotly.
#'
#' @param mdata      data.frame, or named list of data.frames.
#' @param x_col,y_col Columns to plot. Defaults match plotPassedCells().
#' @param id_col     Column holding the cell id. NULL (default) uses rownames,
#'                   which is what the vignette sets to `barcode`.
#' @param log_x      Log10 x axis (default TRUE, as in plotPassedCells).
#' @param ref_col    Existing logical call to show as point colour for
#'                   reference, e.g. "passedMB" (the Mclust call) or "passedMZ".
#'                   Set NULL for a single colour.
#' @param out_col    Name of the logical column written back (default
#'                   "passedManual").
#' @param point_size Marker size.
#' @param max_points Safety valve for very large panels. If a panel has more
#'                   cells than this it is randomly downsampled *and only the
#'                   displayed cells can be gated* -- the rest are returned
#'                   FALSE. Default Inf (no downsampling): the WebGL renderer
#'                   handles the ~35k cells per panel in this project without
#'                   trouble. Lower it only if the gadget feels sluggish, and
#'                   be aware of what it costs you.
#'
#' @return Same structure as `mdata`, with `out_col` added.
#'
#' @examples
#' \dontrun{
#' source("scripts/gateCells.R")
#'
#' # single panel
#' md <- gateCells(metadata.ls[["ATAC_TATAGCCT"]], ref_col = "passedMB")
#' table(md$passedMB, md$passedManual)
#'
#' # all panels in one session
#' metadata.ls <- gateCells(metadata.ls, ref_col = "passedMB")
#' sapply(metadata.ls, function(d) sum(d$passedManual))
#' }
gateCells <- function(mdata,
                      x_col      = "all_unique_MB",
                      y_col      = "peak_ratio_MB",
                      id_col     = NULL,
                      log_x      = TRUE,
                      ref_col    = "passedMB",
                      out_col    = "passedManual",
                      point_size = 4,
                      max_points = Inf) {
  
  for (p in c("shiny", "miniUI", "plotly")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required: install.packages('", p, "')", call. = FALSE)
    }
  }
  
  was_df <- is.data.frame(mdata)
  panels <- if (was_df) list(data = mdata) else mdata
  if (is.null(names(panels)) || any(names(panels) == "")) {
    stop("`mdata` must be a data.frame or a *named* list of data.frames.", call. = FALSE)
  }
  
  # ---- prepare a plotting frame per panel -----------------------------------
  prep <- lapply(names(panels), function(nm) {
    d <- panels[[nm]]
    for (cl in c(x_col, y_col)) {
      if (!cl %in% colnames(d)) {
        stop("Panel '", nm, "' has no column '", cl, "'.", call. = FALSE)
      }
    }
    ids <- if (is.null(id_col)) rownames(d) else as.character(d[[id_col]])
    if (anyDuplicated(ids)) stop("Non-unique cell ids in panel '", nm, "'.", call. = FALSE)
    
    ref <- if (!is.null(ref_col) && ref_col %in% colnames(d)) {
      factor(ifelse(as.logical(d[[ref_col]]), paste0(ref_col, ": TRUE"),
                    paste0(ref_col, ": FALSE")))
    } else {
      factor(rep("cells", nrow(d)))
    }
    
    out <- data.frame(.id = ids,
                      x   = as.numeric(d[[x_col]]),
                      y   = as.numeric(d[[y_col]]),
                      ref = ref,
                      stringsAsFactors = FALSE)
    # log axis cannot show non-positive x
    out <- out[is.finite(out$x) & is.finite(out$y), , drop = FALSE]
    if (log_x) out <- out[out$x > 0, , drop = FALSE]
    out
  })
  names(prep) <- names(panels)
  
  # ---- UI --------------------------------------------------------------------
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Gate cells — lasso or box select"),
    miniUI::miniContentPanel(
      shiny::fillCol(
        flex = c(NA, 1, NA),
        shiny::fillRow(
          height = "60px",
          shiny::selectInput("panel", "Panel", choices = names(prep), width = "260px"),
          shiny::radioButtons("mode", "Selection", inline = TRUE,
                              choices = c("Add" = "add", "Remove" = "remove",
                                          "Replace" = "replace"),
                              selected = "add"),
          shiny::actionButton("clear",  "Clear panel"),
          shiny::actionButton("invert", "Invert panel")
        ),
        plotly::plotlyOutput("plot", height = "100%"),
        shiny::verbatimTextOutput("info")
      )
    )
  )
  
  # ---- server ----------------------------------------------------------------
  server <- function(input, output, session) {
    
    selected <- shiny::reactiveValues()
    for (nm in names(prep)) selected[[nm]] <- character(0)
    
    cur <- shiny::reactive({
      d <- prep[[input$panel]]
      # NB: cells dropped here cannot be gated -- they come back FALSE.
      if (is.finite(max_points) && nrow(d) > max_points) {
        set.seed(1)
        d <- d[sort(sample.int(nrow(d), max_points)), , drop = FALSE]
        shiny::showNotification(
          sprintf("Panel downsampled to %d of %d cells; only these can be gated.",
                  max_points, nrow(prep[[input$panel]])),
          type = "warning", duration = 8)
      }
      d
    })
    
    output$plot <- plotly::renderPlotly({
      d <- cur()
      sel <- selected[[input$panel]]
      d$state <- ifelse(d$.id %in% sel, "selected", as.character(d$ref))
      
      pal <- c("selected" = "#111111")
      lv  <- setdiff(unique(d$state), "selected")
      base_cols <- c("#F8766D", "#00BA38", "#56B4E9", "#E69F00")
      pal[lv] <- base_cols[seq_along(lv)]
      
      plotly::plot_ly(
        d, x = ~x, y = ~y,
        customdata = ~.id,
        color = ~state, colors = pal,
        type = "scattergl", mode = "markers",
        marker = list(size = point_size, opacity = 0.55,
                      line = list(width = 0)),
        hoverinfo = "text",
        text = ~paste0(.id, "<br>", x_col, ": ", signif(x, 4),
                       "<br>", y_col, ": ", signif(y, 3)),
        source = "gate"
      ) |>
        plotly::layout(
          dragmode = "lasso",
          xaxis = list(title = x_col, type = if (log_x) "log" else "linear"),
          yaxis = list(title = paste0(y_col, "  (FRiP)")),
          legend = list(orientation = "h", y = -0.15)
        ) |>
        plotly::config(displaylogo = FALSE,
                       modeBarButtonsToAdd = c("lasso2d", "select2d")) |>
        plotly::event_register("plotly_selected")
    })
    
    shiny::observeEvent(plotly::event_data("plotly_selected", source = "gate"), {
      ev <- plotly::event_data("plotly_selected", source = "gate")
      if (is.null(ev) || !nrow(ev)) return()
      hit <- as.character(ev$customdata)
      old <- selected[[input$panel]]
      selected[[input$panel]] <- switch(
        input$mode,
        add     = union(old, hit),
        remove  = setdiff(old, hit),
        replace = hit
      )
    })
    
    shiny::observeEvent(input$clear,  { selected[[input$panel]] <- character(0) })
    shiny::observeEvent(input$invert, {
      all_ids <- prep[[input$panel]]$.id
      selected[[input$panel]] <- setdiff(all_ids, selected[[input$panel]])
    })
    
    output$info <- shiny::renderText({
      n_tot <- nrow(prep[[input$panel]])
      n_sel <- length(selected[[input$panel]])
      tot <- sum(vapply(names(prep), function(n) length(selected[[n]]), integer(1)))
      sprintf("%s: %d / %d selected (%.1f%%)   |   all panels: %d selected",
              input$panel, n_sel, n_tot, 100 * n_sel / max(n_tot, 1), tot)
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(shiny::isolate(
        stats::setNames(lapply(names(prep), function(n) selected[[n]]), names(prep))
      ))
    })
    shiny::observeEvent(input$cancel, shiny::stopApp(NULL))
  }
  
  res <- shiny::runGadget(ui, server,
                          viewer = shiny::dialogViewer("Gate cells", width = 1000, height = 720))
  if (is.null(res)) {
    message("Gating cancelled -- returning input unchanged.")
    return(mdata)
  }
  
  # ---- write the gate back ----------------------------------------------------
  out <- lapply(names(panels), function(nm) {
    d <- panels[[nm]]
    ids <- if (is.null(id_col)) rownames(d) else as.character(d[[id_col]])
    d[[out_col]] <- ids %in% res[[nm]]
    d
  })
  names(out) <- names(panels)
  
  for (nm in names(out)) {
    message(sprintf("%-24s %6d / %6d cells gated", nm,
                    sum(out[[nm]][[out_col]]), nrow(out[[nm]])))
  }
  
  if (was_df) out[[1]] else out
}


#' Minimal base-graphics alternative (no Shiny)
#'
#' Uses gatepoints::fhs(): click round the cells you want, then right-click /
#' Esc to close the polygon. Blocks the console until you finish. Handy for a
#' single quick gate; `gateCells()` is better for many panels.
#'
#' @examples
#' \dontrun{
#' install.packages("gatepoints")
#' keep <- gateCellsBase(metadata.ls[["ATAC_TATAGCCT"]])
#' }
gateCellsBase <- function(mdata,
                          x_col = "all_unique_MB",
                          y_col = "peak_ratio_MB",
                          log_x = TRUE) {
  if (!requireNamespace("gatepoints", quietly = TRUE)) {
    stop("install.packages('gatepoints')", call. = FALSE)
  }
  xy <- data.frame(x = as.numeric(mdata[[x_col]]),
                   y = as.numeric(mdata[[y_col]]))
  rownames(xy) <- rownames(mdata)
  keep_finite <- is.finite(xy$x) & is.finite(xy$y) & (!log_x | xy$x > 0)
  xy <- xy[keep_finite, , drop = FALSE]
  if (log_x) xy$x <- log10(xy$x)
  
  plot(xy$x, xy$y, pch = 16, cex = 0.3, col = "#00000055",
       xlab = if (log_x) paste0("log10(", x_col, ")") else x_col,
       ylab = paste0(y_col, " (FRiP)"))
  message("Click to draw the gate; right-click (or Esc) to close it.")
  sel <- gatepoints::fhs(xy, mark = TRUE)
  rownames(mdata) %in% as.character(sel)
}
