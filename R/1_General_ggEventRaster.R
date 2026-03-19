#' Spike Raster for EPhysEvents (auto-layout)
#'
#' Expects \code{X} to be in one of two valid formats:
#' \enumerate{
#'   \item \code{length(\link[EPhysData]{Channels}(X)) == 1}:
#'         plots one channel across all \link[EPhysData]{Metadata} rows (rows = traces).
#'   \item \code{nrow(\link[EPhysData]{Metadata}(X)) == 1}:
#'         plots all \link[EPhysData]{Channels} within the single metadata row (rows = channels).
#' }
#'
#' Spikes are drawn as vertical ticks via \code{\link[ggplot2]{geom_linerange}}.
#' If a stimulus is present and \code{show_stimulus = TRUE}, the output of
#' \code{\link{ggStimulusPlot}} is stacked underneath via
#' \code{\link[cowplot]{plot_grid}} (\pkg{cowplot}).
#'
#' In Case A (one channel across metadata rows), rows can be labeled/grouped by a
#' metadata column via \code{label_col}. If \code{step_gap > 0}, the raster is
#' grouped by the chosen row labels, vertical gaps are inserted between groups,
#' and the y-axis uses group-midpoint ticks (block labels).
#'
#' @param X An \code{\link[EPhysData]{EPhysEvents-class}} object in a valid configuration (see above).
#' @param tlim Optional numeric length-2 \code{c(tmin, tmax)} time window (s) passed to
#'   \code{\link[ggplot2]{coord_cartesian}}.
#' @param tick_height Numeric height of each spike tick in y-units (row spacing). Default 0.8.
#' @param line_size Line width for \code{\link[ggplot2]{geom_linerange}}. Default 0.3.
#' @param alpha_raster Numeric in \eqn{[0,1]}. Opacity of spike ticks in
#'   \code{\link[ggplot2]{geom_linerange}}. Default 1.
#' @param show_stimulus Logical; if TRUE, compose a stimulus panel underneath using
#'   \code{\link{ggStimulusPlot}} and \code{\link[cowplot]{plot_grid}}. Default TRUE.
#' @param stimulus_height Relative height of the stimulus panel in
#'   \code{\link[cowplot]{plot_grid}}. Default 0.6.
#' @param row_labels Optional character vector of row labels for Case A; if provided,
#'   this overrides \code{label_col}.
#' @param label_col Optional character scalar naming a \code{\link[EPhysData]{Metadata}} column to use as
#'   row labels in Case A when \code{row_labels} is \code{NULL}. Ignored in Case B.
#' @param step_gap Non-negative numeric. If \code{> 0} (Case A only), inserts a
#'   vertical gap of this size (in y-units) between successive label groups and
#'   switches the y-axis to group-midpoint ticks (block labels) via
#'   \code{\link[ggplot2]{scale_y_continuous}}.
#' @param step_order_decreasing Logical. If TRUE (Case A only), reverses the order
#'   of label groups. Ordering is numeric if labels parse as numeric, otherwise
#'   lexicographic.
#' @param background NULL (default) or the name of a \code{\link[EPhysData]{Metadata}} column.
#'   If non-NULL, a background tile is drawn for each raster row with
#'   \code{\link[ggplot2]{geom_tile}}, and \code{alpha} mapped to that column.
#' @param background_range Optional numeric length-2 range used to draw the background tile
#'   along x (typically the stimulus-encoding rectangle window).
#' @param background_fill The fill colour for the background tiles (used in \code{\link[ggplot2]{geom_tile}}).
#'
#' @return A \code{ggplot} object (or a \code{cowplot} composite). The tidy spike data
#'   are attached as \code{attr(., "raster_data")}.
#'
#' @examples
#' \dontrun{
#' ggEventRaster(ephys, tlim = c(-0.5, 1), show_stimulus = FALSE)
#' ggEventRaster(ephys, tlim = c(0, 2), label_col = "Intensity", step_gap = 1,
#'              step_order_decreasing = TRUE, show_stimulus = FALSE)
#' }
#'
#' @seealso
#' \code{\link[EPhysData]{Metadata}}, \code{\link[EPhysData]{Channels}},
#' \code{\link{ggStimulusPlot}},
#' \code{\link[ggplot2]{geom_linerange}}, \code{\link[ggplot2]{geom_tile}},
#' \code{\link[cowplot]{plot_grid}}
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @importFrom ggpubr theme_pubr
#' @importFrom ggplot2 ggplot aes geom_blank geom_tile geom_linerange guide_axis
#'   scale_y_continuous expansion coord_cartesian labs theme element_blank
#' @name ggEventRaster
setGeneric("ggEventRaster", function(X, ...)
  standardGeneric("ggEventRaster"))

#' @export
setMethod("ggEventRaster", signature(X = "EPhysEvents"),
          function(X,
                   tlim                = NULL,
                   tick_height         = 0.8,
                   line_size           = 0.3,
                   show_stimulus       = TRUE,
                   stimulus_height     = 0.6,
                   row_labels          = NULL,
                   label_col           = NULL,
                   step_gap            = 0,
                   step_order_decreasing = FALSE,
                   background          = NULL,
                   background_range    = NULL,
                   background_fill     = "black",
                   alpha_raster        = 1) {

            stopifnot(is.numeric(step_gap), length(step_gap) == 1L, is.finite(step_gap), step_gap >= 0)

            md  <- Metadata(X)
            chs <- Channels(X)
            n_rows <- if (!is.null(md)) nrow(md) else length(X@Data)

            # optional background column check
            if (!is.null(background)) {
              if (is.null(md)) stop("`background` was specified, but Metadata(X) is NULL.")
              if (!background %in% names(md)) stop("background column '", background, "' not found in Metadata(X).")
              if (is.null(background_range) || length(background_range) != 2L || !all(is.finite(background_range))) {
                stop("`background_range` must be a finite length-2 numeric vector when `background` is specified.")
              }
            }

            # Determine layout case:
            caseA <- !is.null(chs) && length(chs) == 1L          # one channel across rows
            caseB <- !is.null(md)  && n_rows == 1L               # one row with many channels

            if (!(caseA || caseB)) {
              stop(
                "Object must satisfy one of:\n",
                "  • length(Channels(X)) == 1, or\n",
                "  • nrow(Metadata(X)) == 1."
              )
            }

            # Helpers
            channels_in_row <- function(x, i) {
              nm <- names(x@Data[[i]])
              if (is.null(nm)) seq_along(x@Data[[i]]) else nm
            }
            extract_spike_times <- function(cell) {
              if (is.null(cell)) return(numeric(0))
              if (is.numeric(cell)) return(as.numeric(cell))
              if (is.data.frame(cell)) {
                if ("Time" %in% names(cell)) return(as.numeric(cell$Time))
                if ("t"    %in% names(cell)) return(as.numeric(cell$t))
                return(as.numeric(cell[[1]]))
              }
              if (is.list(cell)) return(as.numeric(unlist(cell, use.names = FALSE)))
              stop("Unsupported spike container type: ", class(cell)[1])
            }

            df <- NULL
            y_tick_breaks <- NULL

            # ---- Case A: one channel across all Metadata rows
            if (caseA) {
              ch1 <- chs[[1]]
              missing_in <- integer(0)
              for (i in seq_len(n_rows)) {
                row_ch <- channels_in_row(X, i)
                if (!(ch1 %in% row_ch)) { missing_in <- c(missing_in, i); next }
                spikes <- extract_spike_times(X@Data[[i]][[ch1]])
                if (!is.null(tlim)) spikes <- spikes[spikes >= tlim[1] & spikes <= tlim[2]]

                df <- rbind(df, if (length(spikes)) {
                  data.frame(runuid = md$RunUID[i], step = md$Step[i], unit = i, label = NA_character_,
                             channel = ch1, t = spikes, stringsAsFactors = FALSE)
                } else {
                  data.frame(runuid = md$RunUID[i], step = md$Step[i], unit = i, label = NA_character_,
                             channel = ch1, t = NA_real_, stringsAsFactors = FALSE)
                })
              }
              if (length(missing_in)) {
                stop("Channel '", ch1, "' not present in rows: ", paste(missing_in, collapse = ", "), ".")
              }

            } else {
              # ---- Case B: one Metadata row with multiple channels
              i <- 1L
              ch_use <- channels_in_row(X, i)
              for (ch in ch_use) {
                spikes <- extract_spike_times(X@Data[[i]][[ch]])
                if (!is.null(tlim)) spikes <- spikes[spikes >= tlim[1] & spikes <= tlim[2]]

                df <- rbind(df, if (length(spikes)) {
                  data.frame(runuid = md$RunUID[i], step = md$Step[i], unit = ch, label = as.character(ch),
                             channel = ch, t = spikes, stringsAsFactors = FALSE)
                } else {
                  data.frame(runuid = md$RunUID[i], step = md$Step[i], unit = ch, label = as.character(ch),
                             channel = ch, t = NA_real_, stringsAsFactors = FALSE)
                })
              }
            }

            # ---- Build row labels (Case A) / channel labels (Case B)
            # ---- Shared auto-layout (labels + y mapping), replaces lines 178–270 ----
            lay <- .ephys_autolayout(
              X,
              row_labels = row_labels,
              label_col = label_col,
              step_gap = step_gap,
              step_order_decreasing = step_order_decreasing
            )

            df$label <- if (isTRUE(lay$caseA)) unname(lay$rowlabs[as.character(df$unit)]) else as.character(df$unit)
            df$label[is.na(df$label) | !nzchar(df$label)] <- "NA"

            df$y <- unname(lay$y_map[as.character(df$unit)])

            y_axis_title <- lay$y_axis_title
            y_breaks     <- lay$y_breaks
            y_labels     <- lay$y_labels
            y_tick_breaks <- lay$y_tick_breaks

            # ---- x-limits
            if (is.null(tlim)) {
              xlim_use <- if (nrow(df) && any(is.finite(df$t))) range(df$t, na.rm = TRUE) else c(0, 1)
              if (!all(is.finite(xlim_use))) xlim_use <- c(0, 1)
            } else {
              xlim_use <- tlim
            }

            # ---- Optional background tiles
            bg_df <- NULL
            if (!is.null(background)) {
              bg_meta <- md[, c("RunUID", background), drop = FALSE]
              names(bg_meta)[2] <- "bg_value"

              df <- merge(df, bg_meta,
                          by.x  = "runuid",
                          by.y  = "RunUID",
                          all.x = TRUE,
                          sort  = FALSE)

              bg_df <- unique(df[, c("y", "bg_value")])
              bg_df$x      <- mean(background_range)
              bg_df$width  <- diff(background_range)
              bg_df$height <- 1
            }

            # Blank scaffold so rows without spikes still reserve space
            blank_df <- data.frame(t = xlim_use[1], y = sort(unique(df$y)))

            # ---- Plot
            p_raster <- ggplot(df, aes(x = t)) +
              geom_blank(data = blank_df, aes(x = t, y = y))

            if (!is.null(bg_df)) {
              p_raster <- p_raster +
                geom_tile(
                  data = bg_df,
                  aes(x = x, y = y, alpha = bg_value),
                  fill        = background_fill,
                  width       = bg_df$width[1],
                  height      = bg_df$height[1],
                  inherit.aes = FALSE
                )
            }

            p_raster <- p_raster +
              geom_linerange(
                aes(ymin = y - tick_height/2, ymax = y + tick_height/2),
                alpha=alpha_raster,
                linewidth = line_size, na.rm = TRUE
              ) +
              scale_y_continuous(
                breaks = y_breaks,
                labels = y_labels,
                expand = expansion(mult = c(0.02, 0.02)),
                minor_breaks = y_tick_breaks
              ) +
              coord_cartesian(xlim = xlim_use) +
              labs(
                x = "Time (s)",
                y = y_axis_title,
                title = if (caseA)
                  paste0("Raster: Channel ", unique(df$channel[is.finite(df$t)])[1])
                else
                  "Raster: Single trace, per channel"
              ) +
              theme_pubr(base_size = 8) +
              guides(y = guide_axis(minor.ticks = TRUE)) +
              theme(
                # panel.grid.major.y = element_blank(),
                # panel.grid.minor   = element_blank(),
                # remove gridlines at minor breaks if you don’t want them
                panel.grid.minor.y = element_line(),
                axis.ticks.y.left = element_blank(),
                axis.minor.ticks.y.left = element_line(),
                axis.minor.ticks.length = ggplot2::rel(2)
              )

            p_out <- p_raster

            # ---- Optional stimulus panel via cowplot
            if (isTRUE(show_stimulus)) {
              stim_plot <- tryCatch({
                ggStimulusPlot(X, tlim = xlim_use)
              }, error = function(e) NULL)

              if (!is.null(stim_plot)) {
                if (requireNamespace("cowplot", quietly = TRUE)) {
                  p_out <- cowplot::plot_grid(
                    p_raster, stim_plot, ncol = 1,
                    rel_heights = c(1, stimulus_height), align = "v", axis = "lr"
                  )
                } else {
                  p_out <- p_raster
                }
              }
            }

            p_out
          })
