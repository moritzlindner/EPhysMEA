#' PSTH heatmaps for EPhysContinuous objects and pre-binned matrices
#'
#' \code{ggPSTHHeatmap()} draws PSTH-like heatmaps from an
#' \code{\link[EPhysData:EPhysContinuous-class]{EPhysContinuous-class}} object
#' using the same auto-layout logic as \code{\link{ggEventRaster}}.
#'
#' Valid input configurations for \code{ggPSTHHeatmap()} are:
#' \enumerate{
#'   \item \code{length(\link[EPhysData:Channels]{Channels}(X)) == 1}: one channel
#'   across all \code{\link[EPhysData:Metadata]{Metadata}(X)} rows
#'   (heatmap rows correspond to metadata rows / runs).
#'   \item \code{nrow(\link[EPhysData:Metadata]{Metadata}(X)) == 1}: all
#'   \code{\link[EPhysData:Channels]{Channels}(X)} within a single metadata row
#'   (heatmap rows correspond to channels).
#' }
#'
#' \code{X@Data} is assumed to be arranged as \code{[time × run × channel]}, as
#' produced by \code{\link{Bin}} and optionally \code{\link{AverageTrials}}.
#'
#' \code{ggPSTHHeatmap.primitive()} exposes the same plotting backend for a
#' pre-binned numeric matrix plus matching metadata. In that interface, time-bin
#' duration and time-unit labeling are supplied explicitly, and no
#' \code{EPhysContinuous} object is required.
#'
#' @param X For \code{ggPSTHHeatmap()} only: an
#'   \code{\link[EPhysData:EPhysContinuous-class]{EPhysContinuous-class}} object
#'   in a valid configuration.
#' @param tlim For \code{ggPSTHHeatmap()} only: optional numeric length-2 vector
#'   \code{c(tmin, tmax)} specifying the time window, in the units returned by
#'   \code{\link[EPhysData:TimeUnits]{TimeUnits}(X)}.
#' @param value Character scalar controlling fill labeling. \code{"freq"} yields
#'   the legend title \dQuote{Rate (1/s)}; \code{"count"} yields \dQuote{Count}.
#'   In \code{ggPSTHHeatmap()}, \code{"freq"} additionally converts counts to
#'   rates by dividing by bin duration. In \code{ggPSTHHeatmap.primitive()}, the
#'   supplied matrix is plotted as-is, so callers should provide values already in
#'   the desired units.
#' @param normalize Logical. In \code{ggPSTHHeatmap()}, if \code{TRUE}, each row is
#'   min-max normalised over time to the interval \code{[0, 1]}. In
#'   \code{ggPSTHHeatmap.primitive()}, the supplied matrix is plotted as-is, so any
#'   desired normalisation should be performed before calling the function.
#' @param show_stimulus For \code{ggPSTHHeatmap()} only: logical; if \code{TRUE},
#'   stack a stimulus panel underneath using \code{\link{ggStimulusPlot}} when
#'   stimulus information is available.
#' @param stimulus_height For \code{ggPSTHHeatmap()} only: relative height of the
#'   stimulus panel in
#'   \code{\link[cowplot:plot_grid]{cowplot::plot_grid()}}.
#' @param row_labels,label_col,step_gap,step_order_decreasing For
#'   \code{ggPSTHHeatmap()} only: same meaning as in \code{\link{ggEventRaster}}.
#' @param NoRowLabels Logical; if \code{TRUE}, hide y-axis labels and ticks.
#' @param Legend Logical; if \code{FALSE}, hide the heatmap legend.
#' @param PlotTitle Character scalar or \code{NULL} used as plot title.
#' @param bin_mtx For \code{ggPSTHHeatmap.primitive()} only: numeric matrix whose
#'   rows correspond to observations / units and whose columns correspond to time
#'   bins.
#' @param meta For \code{ggPSTHHeatmap.primitive()} only: \code{data.frame} with
#'   \code{nrow(meta) == nrow(bin_mtx)}. Used for row ordering, row labeling and
#'   optional faceting.
#' @param time_unit For \code{ggPSTHHeatmap.primitive()} only: character scalar or
#'   \code{NULL} giving the x-axis time unit label.
#' @param bin_dur For \code{ggPSTHHeatmap.primitive()} only: numeric scalar giving
#'   the duration of one time bin.
#' @param xlim_use For \code{ggPSTHHeatmap.primitive()} only: numeric length-2
#'   vector or \code{NULL}, passed to
#'   \code{\link[ggplot2:coord_cartesian]{ggplot2::coord_cartesian()}}.
#' @param sort_by For \code{ggPSTHHeatmap.primitive()} only: character vector of
#'   column names in \code{meta}. Rows are ordered by these variables, within each
#'   facet when \code{facet_row} is set. Row labels are formed by concatenating
#'   these variables with \dQuote{ · }.
#' @param facet_row For \code{ggPSTHHeatmap.primitive()} only: \code{NULL} or
#'   character scalar naming a column in \code{meta}. If provided, the heatmap is
#'   faceted vertically by this variable.
#' @param ... Further arguments reserved for future use.#'
#'
#' @return
#' \code{ggPSTHHeatmap.primitive()} returns a
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} object.
#'
#' \code{ggPSTHHeatmap()} returns a
#' \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} object and, when
#' \code{show_stimulus = TRUE} and a stimulus panel can be created, may instead
#' return a
#' \code{\link[cowplot:plot_grid]{cowplot::plot_grid()}} composite.
#'
#' @seealso \code{\link{ggPSTHHeatmap.primitive}}, \code{\link{ggPSTHHeatmap.draw}},
#'   \code{\link{ggEventRaster}}, \code{\link{ggStimulusPlot}},
#'   \code{\link{Bin}}, \code{\link{AverageTrials}}
#'
#' @importClassesFrom EPhysData EPhysContinuous
#' @importFrom EPhysData Metadata Channels TimeTrace TimeUnits
#' @importFrom ggplot2 ggplot aes geom_raster scale_y_continuous scale_y_discrete
#'   scale_x_continuous scale_fill_viridis_c expansion coord_cartesian
#'   labs theme element_blank margin facet_grid waiver
#' @importFrom ggpubr theme_pubr
#' @importFrom cowplot plot_grid
#' @export
setGeneric("ggPSTHHeatmap", function(X,
                                     tlim = NULL,
                                     value = c("freq", "count"),
                                     normalize = FALSE,
                                     show_stimulus = TRUE,
                                     stimulus_height = 0.6,
                                     row_labels = NULL,
                                     label_col = NULL,
                                     step_gap = 0,
                                     step_order_decreasing = FALSE,
                                     NoRowLabels = FALSE,
                                     Legend = TRUE,
                                     PlotTitle = "PSTH heatmap") standardGeneric("ggPSTHHeatmap"))

#' @export
setMethod("ggPSTHHeatmap", signature(X = "EPhysContinuous"),
          function(X,
                   tlim = NULL,
                   value = c("freq", "count"),
                   normalize = FALSE,
                   show_stimulus = TRUE,
                   stimulus_height = 0.6,
                   row_labels = NULL,
                   label_col = NULL,
                   step_gap = 0,
                   step_order_decreasing = FALSE,
                   NoRowLabels = FALSE,
                   Legend = TRUE,
                   PlotTitle = "PSTH heatmap") {

            value <- match.arg(value)

            lay <- .ephys_autolayout(
              X,
              row_labels = row_labels,
              label_col = label_col,
              step_gap = step_gap,
              step_order_decreasing = step_order_decreasing
            )

            arr <- X@Data
            stopifnot(is.array(arr), length(dim(arr)) == 3L) # [time × run × channel]

            tt <- to_num(TimeTrace(X))
            if (!length(tt)) stop("TimeTrace(X) is empty.")
            # bin duration from bin starts
            bin_dur <- if (length(tt) >= 2L) stats::median(diff(tt), na.rm = TRUE) else 1
            if (!is.finite(bin_dur) || bin_dur <= 0) bin_dur <- 1
            time_centers <- tt + bin_dur / 2

            # x-limits
            xlim_use <- if (!is.null(tlim)) {
              stopifnot(is.numeric(tlim), length(tlim) == 2L, all(is.finite(tlim)))
              tlim
            } else {
              range(time_centers, finite = TRUE)
            }

            keep_t <- which(time_centers >= xlim_use[1] & time_centers <= xlim_use[2])
            if (!length(keep_t)) stop("No time bins fall within tlim.")

            # ---- extract plotting matrix + unit ids ----
            if (isTRUE(lay$caseA)) {
              # one channel across Metadata rows / runs  -> mat: [unit(run) × time]
              ch_idx <- 1L
              tmp <- arr[keep_t, , ch_idx, drop = TRUE]  # [time × run] or vector
              if (is.null(dim(tmp))) tmp <- matrix(tmp, nrow = length(keep_t), ncol = 1)
              tmp <- as.matrix(tmp)

              mat <- t(tmp)  # [run × time]

              # IMPORTANT: keys must match lay$y_map names (row indices as characters)
              unit_id <- as.character(seq_len(nrow(mat)))

            } else {
              # one Metadata row, many channels -> mat: [channel × time]
              run_idx <- 1L
              tmp <- arr[keep_t, run_idx, , drop = TRUE] # [time × channel] or vector
              if (is.null(dim(tmp))) tmp <- matrix(tmp, nrow = length(keep_t), ncol = 1)
              tmp <- as.matrix(tmp)

              mat <- t(tmp)  # [channel × time]

              # IMPORTANT: keys must match lay$y_map names (channel names)
              unit_id <- as.character(EPhysData::Channels(X))
              rownames(mat) <- unit_id
            }

            # ---- counts -> freq ----
            if (value == "freq") mat <- mat / bin_dur

            # ---- optional per-row normalization ----
            if (isTRUE(normalize)) {
              rmin <- apply(mat, 1, min, na.rm = TRUE)
              rmax <- apply(mat, 1, max, na.rm = TRUE)
              den  <- rmax - rmin
              den[!is.finite(den) | den <= 0] <- NA_real_
              mat  <- (mat - rmin) / den
              mat[is.na(mat)] <- 0
            }

            # ---- long df ----
            nU <- nrow(mat)
            nT <- ncol(mat)

            df <- data.frame(
              unit = rep(unit_id, each = nT),
              t    = rep(time_centers[keep_t], times = nU),
              v    = as.vector(t(mat)),
              stringsAsFactors = FALSE
            )

            df$y <- unname(lay$y_map[as.character(df$unit)])

            time_unit <- tryCatch(as.character(EPhysData::TimeUnits(X))[1], error = function(e) NULL)
            if (!length(time_unit) || is.na(time_unit) || !nzchar(time_unit)) {
              time_unit <- NULL
            }

            p <- ggPSTHHeatmap.draw(
              df          = df,
              lay         = lay,
              bin_dur     = bin_dur,
              xlim_use    = xlim_use,
              time_unit   = time_unit,
              PlotTitle   = PlotTitle,
              normalize   = normalize,
              value       = value,
              NoRowLabels = NoRowLabels,
              Legend      = Legend
            )

            if (isTRUE(show_stimulus)) {
              stim_plot <- tryCatch({
                sp <- ggStimulusPlot(X, NoRowLabels = isTRUE(NoRowLabels))
                if (is.null(sp)) return(NULL)
                sp + ggplot2::coord_cartesian(xlim = xlim_use, expand = FALSE)
              }, error = function(e) NULL)

              if (!is.null(stim_plot) && requireNamespace("cowplot", quietly = TRUE)) {
                p2 <- p + ggplot2::theme(
                  axis.text.x  = ggplot2::element_blank(),
                  axis.ticks.x = ggplot2::element_blank(),
                  axis.title.x = ggplot2::element_blank()
                )
                out <- cowplot::plot_grid(
                  p2, stim_plot, ncol = 1, align = "v", axis = "lr",
                  rel_heights = c(1, stimulus_height)
                )
                #attr(out, "heatmap_data") <- df
                return(out)
              }
            }

            p
          })


#' Draw helper for PSTH heatmaps
#'
#' Internal plotting helper used by \code{\link{ggPSTHHeatmap}} and
#' \code{\link{ggPSTHHeatmap.primitive}}.
#'
#' Expects a long-format data frame with columns \code{t} (time), \code{y}
#' (row position or row identifier) and \code{v} (value mapped to fill). If
#' \code{facet_row} is provided, the plot is faceted vertically with
#' \code{\link[ggplot2:facet_grid]{ggplot2::facet_grid()}} by that variable.
#'
#' @param df A \code{data.frame} with columns \code{t}, \code{y}, \code{v}, and
#'   optionally the column named in \code{facet_row}.
#' @param lay A layout list; must contain \code{y_axis_title}. If \code{df$y} is
#'   numeric, \code{lay$y_breaks} and \code{lay$y_labels} are used for the y-axis.
#' @param bin_dur Numeric scalar giving the time-bin duration.
#' @param xlim_use Numeric length-2 vector or \code{NULL}. Passed to
#'   \code{\link[ggplot2:coord_cartesian]{ggplot2::coord_cartesian()}}.
#' @param time_unit Character scalar or \code{NULL}. If non-empty, the x-axis
#'   label becomes \dQuote{Time (<unit>)}; otherwise it is simply \dQuote{Time}.
#' @param PlotTitle Character scalar or \code{NULL} used as plot title.
#' @param normalize Logical. Used to choose the fill legend title.
#' @param value Character scalar. Used to choose the fill legend title;
#'   \code{"freq"} yields \dQuote{Rate (1/s)}, otherwise \dQuote{Count}.
#' @param NoRowLabels Logical; if \code{TRUE}, hides y-axis tick labels and ticks.
#' @param Legend Logical; if \code{FALSE}, hides the legend.
#' @param facet_row \code{NULL} or character scalar naming a column in \code{df}
#'   used for vertical faceting.
#' @param y_label_fun \code{NULL} or function. Optional post-processing function
#'   applied to y-axis labels.
#'
#' @return A \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} object.
#'
#' @seealso \code{\link{ggPSTHHeatmap}}, \code{\link{ggPSTHHeatmap.primitive}}
#' @keywords internal
ggPSTHHeatmap.draw <- function(df,
                               lay,
                               bin_dur,
                               xlim_use,
                               time_unit = NULL,
                               PlotTitle,
                               normalize = FALSE,
                               value = "count",
                               NoRowLabels = FALSE,
                               Legend = TRUE,
                               facet_row = NULL,
                               y_label_fun = NULL) {

  if (!is.null(time_unit)) {
    stopifnot(is.character(time_unit), length(time_unit) == 1L)
  }

  x_lab <- if (!is.null(time_unit) && nzchar(time_unit)) {
    paste0("Time (", time_unit, ")")
  } else {
    "Time"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = y, fill = v)) +
    ggplot2::geom_raster(interpolate = FALSE) +#ggplot2::geom_tile(width = bin_dur, height = 1, colour = NA, linewidth = 0, na.rm = TRUE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = xlim_use, expand = FALSE) +
    ggplot2::labs(
      x = x_lab,
      y = lay$y_axis_title,
      title = PlotTitle
    ) +
    ggpubr::theme_pubr(base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      legend.position = "bottom"
    )

  if (is.factor(df$y) || is.character(df$y)) {
    p <- p + ggplot2::scale_y_discrete(
      labels = if (is.function(y_label_fun)) y_label_fun else ggplot2::waiver(),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    )
  } else {
    p <- p + ggplot2::scale_y_continuous(
      breaks = lay$y_breaks,
      labels = lay$y_labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    )
  }

  p <- p + ggplot2::scale_fill_viridis_c(
    name = if (isTRUE(normalize)) "Normalized"
    else if (identical(value, "freq")) "Rate (1/s)"
    else "Count"
  )

  if (!is.null(facet_row)) {
    if (!facet_row %in% names(df)) {
      stop("facet_row '", facet_row, "' not found in df.")
    }
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste0("`",facet_row, "` ~ .")),
      scales = "free_y",
      space  = "free_y"
    )
  }

  if (isTRUE(NoRowLabels)) {
    p <- p + ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  if (!isTRUE(Legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}


#' Primitive PSTH heatmap renderer
#'
#' \code{ggPSTHHeatmap.primitive()} is a low-level interface to the same plotting
#' backend used by \code{\link{ggPSTHHeatmap}}. It works directly on a pre-binned
#' matrix plus matching metadata and therefore does not require an
#' \code{\link[EPhysData:EPhysContinuous-class]{EPhysContinuous-class}} object.
#'
#' The matrix supplied in \code{bin_mtx} is plotted as provided. In particular,
#' \code{value = "freq"} and \code{normalize = TRUE} affect legend labeling but do
#' not transform the matrix internally; callers should apply any desired count-to-
#' rate conversion or normalisation before calling this function.
#'
#' @details
#' Row order is determined from \code{meta} via \code{sort_by}, optionally within
#' each \code{facet_row} level. Row labels are generated by concatenating the
#' \code{sort_by} columns with \dQuote{ · }. If \code{colnames(bin_mtx)} are
#' present, they must be numeric and are interpreted as time-bin centers;
#' otherwise time-bin centers are reconstructed from \code{bin_dur}.
#'
#' @seealso \code{\link{ggPSTHHeatmap}}, \code{\link{ggPSTHHeatmap.draw}}
#' @rdname ggPSTHHeatmap
#' @export
ggPSTHHeatmap.primitive <- function(bin_mtx,
                                    meta,
                                    time_unit = NULL,
                                    PlotTitle = NULL,
                                    bin_dur,
                                    xlim_use = NULL,
                                    normalize = FALSE,
                                    value = "count",
                                    NoRowLabels = FALSE,
                                    Legend = TRUE,
                                    sort_by = NULL,
                                    facet_row = NULL,
                                    ...) {

  stopifnot(is.matrix(bin_mtx))
  stopifnot(is.data.frame(meta))
  stopifnot(nrow(bin_mtx) == nrow(meta))

  if (!is.null(time_unit)) {
    stopifnot(is.character(time_unit), length(time_unit) == 1L)
  }

  if (!is.null(facet_row) && !facet_row %in% names(meta)) {
    stop("facet_row '", facet_row, "' not found in meta.")
  }

  if (!is.null(colnames(bin_mtx))) {
    time_centers <- as.numeric(colnames(bin_mtx))
    if (anyNA(time_centers)) {
      stop("colnames(bin_mtx) must be numeric if provided.")
    }
  } else {
    time_centers <- ((seq_len(ncol(bin_mtx)) - 1) * bin_dur) + bin_dur / 2
  }
  nU <- nrow(bin_mtx)
  nT <- ncol(bin_mtx)

  sort_spec <- .parse_sort_by_spec(
    sort_by = sort_by,
    meta_names = names(meta),
    exclude = facet_row
  )

  make_label <- function(m) {
    if (!length(sort_spec$cols)) return(as.character(seq_len(nrow(m))))
    apply(m[, sort_spec$cols, drop = FALSE], 1, function(z) paste(z, collapse = " · "))
  }

  ord <- .order_meta_rows(
    meta = meta,
    sort_spec = sort_spec,
    split_by = facet_row
  )

  meta_o <- meta[ord, , drop = FALSE]
  bin_o  <- bin_mtx[ord, , drop = FALSE]

  row_lab <- make_label(meta_o)

  row_id <- seq_len(nrow(meta_o))
  y_id <- paste0(row_id, "||", row_lab)
  y_label_fun <- function(x) sub("^[^|]+\\|\\|", "", x)

  y_fac <- factor(y_id, levels = unique(y_id))

  df <- data.frame(
    t = rep(time_centers, times = nU),
    y = rep(y_fac, each = nT),
    v = as.vector(t(bin_o)),
    stringsAsFactors = FALSE
  )

  if (!is.null(facet_row)) {
    df[[facet_row]] <- rep(meta_o[[facet_row]], each = nT)
  }

  lay <- list(y_axis_title = "Units")

  ggPSTHHeatmap.draw(
    df          = df,
    lay         = lay,
    bin_dur     = bin_dur,
    xlim_use    = xlim_use,
    time_unit   = time_unit,
    PlotTitle   = PlotTitle,
    normalize   = normalize,
    value       = value,
    NoRowLabels = NoRowLabels,
    Legend      = Legend,
    facet_row   = facet_row,
    y_label_fun = y_label_fun
  )
}
