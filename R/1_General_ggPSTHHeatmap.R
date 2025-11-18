

# ---- Main plotting method (by channel; single recording) -----------------
#' ggPSTHHeatmap for EPhysContinuous with Stimulus Overlay (+ optional right-side dot strip)
#'
#' Plots a heatmap of the trial-averaged PSTH from an \code{EPhysContinuous}
#' object, and—if a \code{StimulusTrace} is present—a stimulus trace aligned in time.
#' Optionally adds a right-side dot strip (one dot per channel) showing a
#' per-channel position value.
#'
#' @inheritParams AddStimulus
#' @param normalize  Logical; if \code{TRUE}, each channel’s PSTH is min–max normalized to [0,1].
#' @param sort_by    Optional ordering for channels. One of:
#'   \itemize{
#'     \item \code{NULL} (no reordering; default),
#'     \item a character vector of channel names (unspecified are appended in original order),
#'     \item an integer vector of channel indices (unspecified are appended),
#'     \item the string \code{"variance"} (order by descending PSTH variance).
#'   }
#' @param NoRowLabels Logical; hide y-axis labels/ticks if \code{TRUE}.
#' @param Legend     Logical; hide heatmap legend if \code{FALSE}.
#' @param PlotTitle  Title above the heatmap.
#' @param trim       Numeric in [0, 0.5]; fraction of bins trimmed from each end before averaging.
#' @param sample_channels FALSE (default) or numerif if only to plot \code{sample_channels} random samples.
#' @param dot_position Either:
#'   \itemize{
#'     \item \code{NULL}: no dotbar,
#'     \item a single numeric: same object-position for all dots,
#'     \item a numeric vector of length \code{length(Channels(object))}: per-channel object-positions,
#'           assumed in the same order as \code{Channels(object)}; re-ordered if \code{sort_by} changes the order.
#'   }
#' @param dot_size   Numeric (scalar or vector length \code{length(Channels(object))}). If scalar,
#'                   a fixed point size is used. If a vector, sizes are mapped with a compressed range
#'                   so the maximum appears ~20\% smaller than ggplot2's default max.
#' @param dot_color  Character (scalar or vector length \code{length(Channels(object))}).
#'                   Colors are taken as-is (no legend).
#' @param dotbar_width Numeric; right-side dotbar width weight relative to 100.
#'                     The top/bottom rows use \code{rel_widths = c(100 - dotbar_width, dotbar_width)}.
#'
#' @return A \pkg{cowplot} object.
#' @import ggplot2
#' @importFrom cowplot plot_grid theme_nothing
#' @importClassesFrom EPhysData EPhysContinuous
#' @importFrom EPhysData Metadata Channels ChannelMetadata TimeTrace StimulusTrace Subset as.data.frame HasStimulus
#' @name ggPSTHHeatmap
#' @export
setGeneric("ggPSTHHeatmap",
           function(X,
                    normalize    = FALSE,
                    sort_by      = NULL,
                    NoRowLabels  = FALSE,
                    Legend       = TRUE,
                    PlotTitle    = "Average PSTH Heatmap",
                    trim         = 0,
                    sample_channels = FALSE,
                    dot_position = 1,
                    dot_size     = 1,
                    dot_color    = "grey50",
                    dotbar_width = 2)
           standardGeneric("ggPSTHHeatmap"))

#' @describeIn ggPSTHHeatmap Method for EPhysContinuous
#' @export
setMethod("ggPSTHHeatmap",
          signature(X = "EPhysContinuous"),
          function(X,
                   normalize    = FALSE,
                   sort_by      = NULL,
                   NoRowLabels  = FALSE,
                   Legend       = TRUE,
                   PlotTitle    = "Average PSTH Heatmap",
                   trim         = 0,
                   sample_channels = FALSE,
                   dot_position = 1,
                   dot_size     = 1,
                   dot_color    = "grey50",
                   dotbar_width = 2) {

            if (is.numeric(sample_channels)) {
              keep <- sample(Channels(X), size = sample_channels)
              X <- Subset(X, Channels = keep)
            } else {
              if (!is.logical(sample_channels)){
                stop("sample_channels must be false or numeric")
              } else {
                if(isTRUE(sample_channels)){
                  stop("sample_channels must be false or numeric")
                }
              }
            }

            rel_data_height<-rel_data_height(length(Channels(X)))

            ## ---- Prepare data (average trials) ----
            if(length(unique(Metadata(X)$RecordingID))!=1){
              stop("Data seem to contain Trials from multiple recordings")
            }
            avg_obj <- AverageTrials(X, trim = trim)
            donorm<-as.logical(normalize)
            if (isTRUE(donorm)) {
              avg_obj <- Normalize(avg_obj)
            }
            df <- as.data.frame(avg_obj, ReturnAs = "freq")
            keep_cols <- c("Channel", "Time", "Value")
            if ("Stimulus" %in% names(df)) keep_cols <- c(keep_cols, "Stimulus")
            df <- df[, keep_cols, drop = FALSE]

            ## Base channel order = X order
            ch_orig <- as.character(Channels(X))
            df$Channel <- factor(as.character(df$Channel), levels = ch_orig)

            ## ---- Sorting ----
            if (!is.null(sort_by)) {
              ord <- SortChannels(avg_obj, sort_by = sort_by )
              df$Channel <- factor(as.character(df$Channel), levels = ord)
            }
            ch_final <- levels(df$Channel)  # final channel order

            time_u<-deparse_unit(df$Time)
            val_u<-deparse_unit(df$Value)
            df$Value<-drop_units(df$Value)
            df$Time<-drop_units(df$Time)

            ## ---- Heatmap ----


            heat <- ggplot(df, aes(x = Time, y = Channel, fill = Value)) +
              geom_tile() +
              labs(x = paste0("Time (", time_u, ")"),
                   y = "Channel",
                   title = PlotTitle) +
              theme_minimal() +
              theme(
                panel.grid  = element_blank(),
                plot.margin = margin(0,0,0,0)
              ) +
              scale_y_discrete(expand = c(0,0))+
              scale_x_continuous(expand = c(0,0))
            if(isTRUE(normalize)){
              heat <- heat +
                scale_fill_viridis_c(name = "Normalized Rate")
            } else {
              heat <- heat +
              scale_fill_viridis_c(name = "Rate (Hz)")
            }
            if (isTRUE(NoRowLabels)) {
              heat <- heat + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
            }
            if (!isTRUE(Legend)) {
              heat <- heat + theme(legend.position = "none")
            }

            ## ---- Stimulus (optional) ----
            if (HasStimulus(X)) {
              stim_plot <- ggStimulusPlot(X, NoRowLabels = NoRowLabels)
              heat <- heat + theme(axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank())
            }

            ## ---- Right-side dot strip (optional) ----
            add_dotbar <-
              any(length(dot_position) > 1, length(dot_size) > 1, length(dot_color) > 1)
            if (add_dotbar) {
              dot_info <- coerce_dot_inputs_channels(X, ch_final, dot_position, dot_size, dot_color)
              dot_plot <- build_dotbar_plot(ch_final, dot_info, y_label = "Channel")

              top_row <- cowplot::plot_grid(
                heat, dot_plot,
                ncol = 2, align = "h",
                rel_widths = c(100 - dotbar_width, dotbar_width)
              )

              if (HasStimulus(X)) {
                spacer <- ggplot2::ggplot() + cowplot::theme_nothing()
                bottom_row <- cowplot::plot_grid(
                  stim_plot, spacer,
                  ncol = 2, align = "h",
                  rel_widths = c(100 - dotbar_width, dotbar_width)
                )
                return(cowplot::plot_grid(top_row, bottom_row, ncol = 1, align = "v",
                                          rel_heights = c(rel_data_height(length(ch_final)), 1)))
              } else {
                return(top_row)
              }
            }

            ## ---- No dotbar: just heat (and optional stim below) ----
            if (HasStimulus(X)) {
              upper <- heat + theme(axis.title.x = element_blank(),
                                    axis.text.x  = element_blank(),
                                    axis.ticks.x = element_blank())
              return(plot_grid(upper, stim_plot, ncol = 1, align = "v",
                               rel_heights = c(rel_data_height(length(ch_final)), 1)))
            } else {
              return(heat)
            }
          })




# ---- Internal: per-recording PSTH builder --------------------------------
#' Build per-recording PSTH (assumes a single channel)
#'
#' Splits by RecordingID, averages trials (via \code{AverageTrials()}),
#' converts to data.frame (ReturnAs = "freq") and returns a stacked data.frame
#' with \strong{units}-typed \code{Time}, \code{Value}, and (if present) \code{Stimulus}.
#' No averaging across channels is performed; instead this method \emph{requires}
#' that \code{X} contains exactly one channel. The caller must subset beforehand.
#'
#' @param X         An \code{EPhysContinuous} object (single channel).
#' @param trim      Numeric in [0, 0.5]; trial-end trimming passed to \code{AverageTrials()}.
#' @param normalize Logical; if \code{TRUE}, per-recording min–max normalization on
#'                  \code{Value}. The result keeps \code{units} class with unit \code{1}.
#' @return data.frame with columns:
#'   \item{RecordingID}{Factor in the final order set by the caller.}
#'   \item{Time}{\code{units} vector.}
#'   \item{Value}{\code{units} vector (Hz unless \code{normalize=TRUE} → unitless 1).}
#'   \item{Stimulus}{\code{units} vector, if present in \code{as.data.frame()}.}
#' @name build_rec_psth_df
#' @keywords internal
setGeneric("build_rec_psth_df", function(X, trim = 0.1, normalize = FALSE)
  standardGeneric("build_rec_psth_df"))

#' @describeIn build_rec_psth_df Method for EPhysContinuous
setMethod("build_rec_psth_df",
          signature(X = "EPhysContinuous"),
          function(X, trim = 0.1, normalize = FALSE) {

            md <- Metadata(X)

            # enforce single-channel invariant
            ch <- Channels(X)
            if (length(ch) != 1L)
              stop("`build_rec_psth_df` expects exactly ONE channel in `X`. ",
                   "Subset at the caller (e.g., Subset(X, Channels = ...)).")

            rec_all <- md$RecordingID
            rec_ids <- unique(rec_all)  # original first-appearance order

            out <- vector("list", length(rec_ids))

            for (i in seq_along(rec_ids)) {
              rid <- rec_ids[i]
              xi  <- Subset(X, RecordingID = rid)           # same single channel

              avg <- AverageTrials(xi, trim = trim)    # average trials only
              dfi <- as.data.frame(avg, ReturnAs = "freq")

              # Keep units; do NOT drop here.
              keep <- c("Time", "Value")
              if ("Stimulus" %in% names(dfi)) keep <- c(keep, "Stimulus")
              dfi <- dfi[, keep, drop = FALSE]

              # Optional normalization (min–max) on units-valued column:
              if (isTRUE(normalize)) {
                vnum <- drop_units(dfi$Value)
                rng  <- range(vnum, finite = TRUE)
                if (diff(rng) > 0) vnum <- (vnum - rng[1]) / diff(rng) else vnum <- vnum * 0
                dfi$Value <- set_units(vnum, 1)  # unitless but still class 'units'
              }

              dfi$RecordingID <- rid
              out[[i]] <- dfi
            }

            df <- do.call(rbind, out)
            # RecordingID will be factored (and ordered) by the caller.
            df
          })

# ---- Main plotting method (by recording; single channel) -----------------

setGeneric("ggPSTHHeatmapByRecording",
           function(X,
                    normalize    = FALSE,
                    sort_by      = NULL,
                    NoRowLabels  = FALSE,
                    Legend       = TRUE,
                    PlotTitle    = "Average PSTH per Recording",
                    trim         = 0.1,
                    dot_position = NULL,
                    dot_size     = 1,
                    dot_color    = "grey50",
                    dotbar_width = 2)
             standardGeneric("ggPSTHHeatmapByRecording"))

#' ggPSTHHeatmapByRecording for EPhysContinuous (+ optional stimulus & right dot strip)
#'
#' Heatmap with one row per \code{RecordingID}, using the single channel present
#' in \code{X}. Trials are averaged via \code{AverageTrials()}, then plotted.
#' Optionally normalizes per recording, sorts rows via \code{SortRun()}, and
#' adds a right-side dot strip keyed per recording.
#'
#' @inheritParams as.data.frame.EPhysContinuous-method
#' @param normalize  Logical; if \code{TRUE}, per-recording min–max normalization (unitless \code{1}).
#' @param sort_by    Passed to \code{SortRun(X, sort_by)} to derive row order.
#' @param NoRowLabels,Legend,PlotTitle,trim See channel-based variant.
#' @param dot_position,dot_size,dot_color,dotbar_width See channel-based variant; all interpreted per recording.
#' @return A \pkg{cowplot} object.
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_minimal theme element_blank
#' @importFrom ggplot2 margin scale_y_discrete scale_x_continuous scale_fill_viridis_c
#' @importFrom cowplot plot_grid theme_nothing
#' @importFrom units deparse_unit drop_units
#' @importFrom cowplot plot_grid theme_nothing
#' @describeIn ggPSTHHeatmap ggPSTHHeatmapByRecording
#' @export
setMethod("ggPSTHHeatmapByRecording",
          signature(X = "EPhysContinuous"),
          function(X,
                   normalize    = FALSE,
                   sort_by      = NULL,
                   NoRowLabels  = FALSE,
                   Legend       = TRUE,
                   PlotTitle    = "Average PSTH per Recording",
                   trim         = 0,
                   dot_position = NULL,
                   dot_size     = 1,
                   dot_color    = "grey50",
                   dotbar_width = 2) {

            # Build per-recording df with UNITS columns (single channel invariant)
            df <- build_rec_psth_df(X, trim = trim, normalize = normalize)

            # Decide final RecordingID order using SortRun()
            md <- Metadata(X)
            run_order_idx <- SortRun(as(X, "EPhysContainer"), sort_by = sort_by)
            ord_rec <- unique(md$RecordingID[run_order_idx])
            if (!length(ord_rec)) {
              # fallback: original first-appearance order in Metadata
              ord_rec <- unique(md$RecordingID)
            }
            # Keep only those present (after sampling later)
            # (We set levels after optional sampling.)

            # Finalize factor levels
            present <- unique(df$RecordingID)
            final_levels <- intersect(ord_rec, present)
            if (!length(final_levels)) final_levels <- present
            df$RecordingID <- factor(df$RecordingID, levels = final_levels)

            # Capture units for axis/legend labels before dropping for ggplot
            time_u <- deparse_unit(df$Time)
            val_u  <- deparse_unit(df$Value)

            # ggplot does not handle units -> drop only for plotting aesthetics
            df$Time  <- drop_units(df$Time)
            df$Value <- drop_units(df$Value)
            if ("Stimulus" %in% names(df)) {
              df$Stimulus <- drop_units(df$Stimulus)
            }

            # Heatmap
            heat <- ggplot(df, aes(x = Time, y = RecordingID, fill = Value)) +
              geom_tile() +
              labs(x = paste0("Time (", time_u, ")"),
                            y = "Recording",
                            title = PlotTitle) +
              theme_minimal() +
              theme(
                panel.grid  = element_blank(),
                plot.margin = margin(0,0,0,0)
              ) +
              scale_y_discrete(expand = c(0,0))+
              scale_x_continuous(expand = c(0,0))

            if (isTRUE(normalize)) {
              heat <- heat + scale_fill_viridis_c(name = "Normalized Rate")
            } else {
              lab <- if (is.null(val_u) || val_u == "1") "Rate" else paste0("Rate (", val_u, ")")
              heat <- heat + scale_fill_viridis_c(name = lab)
            }

            if (isTRUE(NoRowLabels)) {
              heat <- heat + theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank())
            }
            if (!isTRUE(Legend)) {
              heat <- heat + theme(legend.position = "none")
            }

            # Stimulus (from full object; time-aligned)
            if (HasStimulus(X)) {
              stim_plot <- ggStimulusPlot(X, NoRowLabels = NoRowLabels)
              heat <- heat + theme(axis.text.x = element_blank(),
                                            axis.ticks.x = element_blank())
            }

            # Right-side dot strip?
            add_dotbar <- !is.null(dot_position) &&
              (
                length(dot_size) > 1L ||
                  length(dot_color) > 1L ||
                  (is.character(dot_size) && length(dot_size) == 1L) ||
                  (is.character(dot_color) && length(dot_color) == 1L)
              )

            if (add_dotbar) {
              rec_final <- levels(df$RecordingID)
              dot_info  <- .coerce_dot_inputs(X, rec_final, dot_position, dot_size, dot_color)
              dot_plot  <- build_dotbar_plot(rec_final, dot_info, y_label = "Recording")

              top_row <- plot_grid(
                heat, dot_plot,
                ncol = 2, align = "h",
                rel_widths = c(100 - dotbar_width, dotbar_width)
              )

              if (HasStimulus(X)) {
                spacer <- ggplot() + theme_nothing()
                bottom_row <- plot_grid(
                  stim_plot, spacer,
                  ncol = 2, align = "h",
                  rel_widths = c(100 - dotbar_width, dotbar_width)
                )
                return(plot_grid(top_row, bottom_row, ncol = 1, align = "v",
                                          rel_heights = c(rel_data_height(length(rec_final)), 1)))
              } else {
                return(top_row)
              }
            }

            # No dotbar
            if (HasStimulus(X)) {
              upper <- heat + theme(axis.title.x = element_blank(),
                                             axis.text.x  = element_blank(),
                                             axis.ticks.x = element_blank())
              return(plot_grid(upper, stim_plot, ncol = 1, align = "v",
                                        rel_heights = c(rel_data_height(length(levels(df$RecordingID))), 1)))
            } else {
              return(heat)
            }
          })



# ---- Helpers for dot-strip inputs (unchanged logic; works per RecordingID) ----

meta_first_per_recording <- function(X, col) {
  md <- Metadata(X)
  rids <- md$RecordingID
  split_vals <- split(md[[col]], rids)
  sapply(split_vals, function(v) {
    w <- which(!is.na(v))
    if (length(w)) v[w[1]] else NA
  })
}

.coerce_dot_inputs <- function(X, rec_ids, dot_position, dot_size, dot_color) {
  # Position
  if (is.null(dot_position)) {
    pos <- NULL
  } else if (is.character(dot_position) && length(dot_position) == 1L) {
    pos <- meta_first_per_recording(X, dot_position)
    pos <- as.numeric(pos[rec_ids])
  } else if (length(dot_position) == 1L && is.numeric(dot_position)) {
    pos <- rep(dot_position, length(rec_ids))
  } else {
    if (!is.numeric(dot_position))
      stop("dot_position must be numeric, NULL, or a Metadata column name.")
    pos <- if (!is.null(names(dot_position))) as.numeric(dot_position[rec_ids]) else {
      if (length(dot_position) != length(rec_ids))
        stop("dot_position must be length 1 or length(n_recordings).")
      as.numeric(dot_position)
    }
  }

  # Size
  if (is.character(dot_size) && length(dot_size) == 1L) {
    siz <- meta_first_per_recording(X, dot_size)
    siz <- as.numeric(siz[rec_ids])
  } else if (length(dot_size) == 1L && is.numeric(dot_size)) {
    siz <- rep(as.numeric(dot_size), length(rec_ids))
  } else {
    if (!is.numeric(dot_size))
      stop("dot_size must be numeric or a Metadata column name.")
    siz <- if (!is.null(names(dot_size))) as.numeric(dot_size[rec_ids]) else {
      if (length(dot_size) != length(rec_ids))
        stop("dot_size must be length 1 or length(n_recordings).")
      as.numeric(dot_size)
    }
  }

  # Color
  if (is.character(dot_color) && length(dot_color) == 1L && dot_color %in% names(Metadata(X))) {
    colv <- meta_first_per_recording(X, dot_color)
    colv <- as.character(colv[rec_ids])
  } else if (length(dot_color) == 1L) {
    colv <- rep(as.character(dot_color), length(rec_ids))
  } else {
    if (length(dot_color) != length(rec_ids))
      stop("dot_color must be length 1 or length(n_recordings).")
    colv <- as.character(dot_color)
  }

  list(pos = pos, size = siz, color = colv)
}

coerce_dot_inputs_channels <- function(X, ch_final, dot_position, dot_size, dot_color) {
  ch_orig <- as.character(Channels(X))

  # Convenience for pulling from ChannelMetadata by name
  pull_chmeta <- function(x) ChannelMetadata(X, x)

  # POSITION
  if (is.null(dot_position)) {
    pos <- NULL
  } else if (is.character(dot_position) && length(dot_position) == 1L &&
             dot_position %in% colnames(ChannelMetadata(X))) {
    pos <- as.numeric(pull_chmeta(dot_position))
  } else if (length(dot_position) == 1L && is.numeric(dot_position)) {
    pos <- rep(as.numeric(dot_position), length(ch_orig))
  } else {
    if (!is.numeric(dot_position))
      stop("dot_position must be numeric or a ChannelMetadata column name.")
    if (length(dot_position) != length(ch_orig))
      stop("dot_position must be length 1 or length(Channels(X)).")
    pos <- as.numeric(dot_position)
  }

  # SIZE
  if (is.character(dot_size) && length(dot_size) == 1L &&
      dot_size %in% colnames(ChannelMetadata(X))) {
    siz <- as.numeric(pull_chmeta(dot_size))
  } else if (length(dot_size) == 1L && is.numeric(dot_size)) {
    siz <- rep(as.numeric(dot_size), length(ch_orig))
  } else {
    if (!is.numeric(dot_size))
      stop("dot_size must be numeric or a ChannelMetadata column name.")
    if (length(dot_size) != length(ch_orig))
      stop("dot_size must be length 1 or length(Channels(X)).")
    siz <- as.numeric(dot_size)
  }

  # COLOR
  if (is.character(dot_color) && length(dot_color) == 1L &&
      dot_color %in% colnames(ChannelMetadata(X))) {
    colv <- as.character(pull_chmeta(dot_color))
  } else if (length(dot_color) == 1L) {
    colv <- rep(as.character(dot_color), length(ch_orig))
  } else {
    if (length(dot_color) != length(ch_orig))
      stop("dot_color must be length 1 or length(Channels(X)).")
    colv <- as.character(dot_color)
  }

  # Reorder to FINAL channel order
  idx <- match(ch_final, ch_orig)
  list(
    pos   = if (is.null(pos)) NULL else pos[idx],
    size  = siz[idx],
    color = colv[idx]
  )
}


#' Build right-side dotbar plot (generic helper)
#'
#' @param y_values Character vector of y-axis values (levels shown top-to-bottom).
#' @param dot_info List with numeric vectors \code{pos}, \code{size} and a
#'   character vector \code{color}; all length \code{length(y_values)}.
#' @param y_label  Character label for the y-axis (e.g., "Channel" or "Recording").
#' @return A ggplot object.
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous labs theme margin
#' @importFrom cowplot theme_nothing
build_dotbar_plot <- function(y_values, dot_info, y_label = "Recording") {
  side_df <- data.frame(
    Y       = factor(y_values, levels = y_values),
    DotX    = as.numeric(dot_info$pos),
    DotSize = as.numeric(dot_info$size),
    DotColor= dot_info$color,
    stringsAsFactors = FALSE
  )
  ggplot(side_df, aes(x = DotX, y = Y, size = DotSize, color = DotColor)) +
    geom_point(na.rm = TRUE) +
    scale_size_continuous(range = c(0, 3)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_nothing() +
    labs(y = y_label) +
    theme(plot.margin = margin(0, 0, 0, 0))
}
rel_data_height <- function(n_rows) sqrt(sqrt(n_rows)) * 4
