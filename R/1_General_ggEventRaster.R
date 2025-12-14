#' Spike Raster for EPhysEvents (auto-layout)
#'
#' Expects \code{object} to be in one of two valid formats:
#' \enumerate{
#'   \item \code{length(Channels(object)) == 1}: plots one channel across all Metadata rows (rows = traces).
#'   \item \code{nrow(Metadata(object)) == 1}: plots all channels within the single Metadata row (rows = channels).
#' }
#' Spikes are drawn as vertical ticks via \code{geom_linerange}. If a stimulus is
#' present and \code{show_stimulus = TRUE}, the output of \code{ggStimulusPlot()}
#' is stacked underneath via \pkg{cowplot}.
#'
#' @param object An \code{EPhysEvents} object in a valid configuration (see above).
#' @param tlim Optional numeric length-2 \code{c(tmin, tmax)} time window (s).
#' @param tick_height Numeric height of each spike tick in y-units (row spacing). Default 0.8.
#' @param line_size Line width for \code{geom_linerange}. Default 0.3.
#' @param show_stimulus Logical; if TRUE, try to compose a stimulus panel underneath. Default TRUE.
#' @param stimulus_height Relative height of the stimulus panel in \code{cowplot::plot_grid}. Default 0.6.
#' @param row_labels Optional character vector of row labels; otherwise derived from Metadata
#' @param background NULL (default) or the name of a Metadata column.
#'   If non-NULL, a background tile is drawn for each raster row, with
#'   `alpha` mapped to that column.
#' @param background_fill
#'   (RecordingID/rownames) or channel names, depending on layout.
#' @return A \code{ggplot} object (or a \code{cowplot} composite). The tidy spike data
#'   is attached as \code{attr(., "raster_data")}.
#'
#' @examples
#' \dontrun{
#' # Case 1: object has exactly one channel overall
#' ggEventRaster(ephys, tlim = c(-0.5, 1))
#'
#' # Case 2: object has exactly one metadata row (one trace)
#' ggEventRaster(ephys_one_row, tlim = c(0, 2), show_stimulus = FALSE)
#' }
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @importFrom ggpubr theme_pubr
#' @importFrom ggplot2 ggplot aes geom_blank geom_linerange scale_y_continuous expansion coord_cartesian labs theme element_blank
#' @name ggEventRaster
setGeneric("ggEventRaster", function(X, ...)
  standardGeneric("ggEventRaster"))

#' @export
setMethod("ggEventRaster", signature(X = "EPhysEvents"),
          function(X,
                   tlim            = NULL,
                   tick_height     = 0.8,
                   line_size       = 0.3,
                   show_stimulus   = TRUE,
                   stimulus_height = 0.6,
                   row_labels      = NULL,
                   background      = NULL,
                   background_fill = "black") {

            md  <- Metadata(X)
            chs <- Channels(X)
            n_rows <- if (!is.null(md)) nrow(md) else length(X@Data)

            # optional background column check
            if (!is.null(background)) {
              if (is.null(md)) {
                stop("`background` was specified, but Metadata(X) is NULL.")
              }
              if (!background %in% names(md)) {
                stop("background column '", background, "' not found in Metadata(X).")
              }
            }


            # Determine layout case:
            caseA <- !is.null(chs) && length(chs) == 1L                # one channel across rows
            caseB <- !is.null(md)  && n_rows == 1L                     # one row with many channels

            if (!(caseA || caseB)) {
              stop(
                "Object must satisfy  one of:\n",
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
                  data.frame(runuid=md$RunUID[i], step=md$Step[i], unit = i, label = i, channel = ch1, t = spikes, stringsAsFactors = FALSE)
                } else {
                  data.frame(runuid=md$RunUID[i], step=md$Step[i], unit = i, label = i, channel = ch1, t = NA_real_, stringsAsFactors = FALSE)
                })
              }
              if (length(missing_in)) {
                stop("Channel '", ch1, "' not present in rows: ",
                     paste(missing_in, collapse = ", "), ".")
              }

            } else {
              # ---- Case B: one Metadata row with multiple channels
              i <- 1L
              ch_use <- channels_in_row(X, i)
              for (ch in ch_use) {
                spikes <- extract_spike_times(X@Data[[i]][[ch]])
                if (!is.null(tlim)) spikes <- spikes[spikes >= tlim[1] & spikes <= tlim[2]]
                df <- rbind(df, if (length(spikes)) {
                  data.frame(runuid=md$RunUID[i], step=md$Step[i], unit = ch, label = ch, channel = ch, t = spikes, stringsAsFactors = FALSE)
                } else {
                  data.frame(runuid=md$RunUID[i], step=md$Step[i], unit = ch, label = ch, channel = ch, t = NA_real_, stringsAsFactors = FALSE)
                })
              }
            }

            # Build labels
            if (caseA) {
              rows_use <- seq_len(n_rows)
              if (is.null(row_labels)) {
                if (!is.null(md)) {
                  if (!is.null(rownames(md)) && length(rownames(md))) {
                    rowlabs <- setNames(rownames(md)[rows_use], rows_use)
                  } else if ("RecordingID" %in% names(md)) {
                    rowlabs <- setNames(as.character(md$RecordingID[rows_use]), rows_use)
                  } else {
                    rowlabs <- setNames(paste0("row_", rows_use), rows_use)
                  }
                } else {
                  rowlabs <- setNames(paste0("row_", rows_use), rows_use)
                }
              } else {
                stopifnot(length(row_labels) == length(rows_use))
                rowlabs <- setNames(row_labels, rows_use)
              }
            } else {
              rowlabs <- setNames(unique(df$unit), unique(df$unit))  # channel names as labels
            }

            # Map rows to integer y positions
            ## --- Map rows to integer y positions, with Step-order in Case A ----
            if (caseA) {
              rows_use <- unique(df$unit)
              step_by_row <- md$Step[rows_use]
              order_keys  <- rows_use[order(step_by_row, decreasing = T)]
            } else {
              order_keys <- unique(df$unit)
            }

            y_map <- setNames(seq_along(order_keys), order_keys)
            df$y  <- unname(y_map[as.character(df$unit)])


            # x-limits
            if (is.null(tlim)) {
              xlim_use <- if (nrow(df) && any(is.finite(df$t))) range(df$t, na.rm = TRUE) else c(0, 1)
              if (!all(is.finite(xlim_use))) xlim_use <- c(0, 1)
            } else {
              xlim_use <- tlim
            }

            # ---- Optional background tiles ---------------------------------------------
            bg_df <- NULL

            if (!is.null(background)) {
              # join the requested metadata column onto df by runuid
              bg_meta <- md[, c("RunUID", background), drop = FALSE]
              names(bg_meta)[2] <- "bg_value"

              df <- merge(df, bg_meta,
                          by.x  = "runuid",
                          by.y  = "RunUID",
                          all.x = TRUE,
                          sort  = FALSE)

              # one row per y for the tile layer
              bg_df <- unique(df[, c("y", "bg_value")])
              bg_df$x      <- mean(xlim_use)
              bg_df$width  <- diff(xlim_use)
              bg_df$height <- 1
            }

            # Blank scaffold so rows without spikes still reserve space
            blank_df <- data.frame(t = xlim_use[1], y = seq_along(order_keys))

            # ---- Plot
            p_raster <- ggplot(df, aes(x = t)) +
              geom_blank(data = blank_df, aes(x = t, y = y))

            # optional background tiles
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

            # spike ticks as before
            p_raster <- p_raster +
              geom_linerange(
                aes(ymin = y - tick_height/2, ymax = y + tick_height/2),
                linewidth = line_size, na.rm = TRUE
              ) +
              scale_y_continuous(
                breaks = seq_along(order_keys),
                labels = unname(if (caseA) rowlabs[as.character(order_keys)] else rowlabs),
                expand = expansion(mult = c(0.02, 0.02))
              ) +
              coord_cartesian(xlim = xlim_use) +
              labs(
                x = "Time (s)",
                y = if (caseA) "Trace (Metadata row)" else "Channel",
                title = if (caseA)
                  paste0("Raster: Channel ", unique(df$channel[is.finite(df$t)])[1])
                else
                  "Raster: Single trace, per channel"
              ) +
              theme_pubr(base_size = 8) +
              theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor   = element_blank()
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
                  # cowplot not installed: fall back to just the raster (no stop)
                  p_out <- p_raster
                }
              }
            }
            p_out
          })
