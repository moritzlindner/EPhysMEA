#' Spike Raster for EPhysEvents (auto-layout)
#'
#' Expects \code{X} to be in one of two valid formats:
#' \enumerate{
#'   \item \code{length(Channels(X)) == 1}: plots one channel across all Metadata rows (rows = traces).
#'   \item \code{nrow(Metadata(X)) == 1}: plots all channels within the single Metadata row (rows = channels).
#' }
#' Spikes are drawn as vertical ticks via \code{geom_linerange}. If a stimulus is
#' present and \code{show_stimulus = TRUE}, the output of \code{ggStimulusPlot()}
#' is stacked underneath via \pkg{cowplot}.
#'
#' In Case A (one channel across Metadata rows), rows can be labeled/grouped by a
#' Metadata column via \code{label_col}. If \code{step_gap > 0}, the raster is
#' grouped by the chosen row labels, vertical gaps are inserted between groups,
#' and the y-axis uses group-midpoint ticks (block labels).
#'
#' @param X An \code{EPhysEvents} object in a valid configuration (see above).
#' @param tlim Optional numeric length-2 \code{c(tmin, tmax)} time window (s).
#' @param tick_height Numeric height of each spike tick in y-units (row spacing). Default 0.8.
#' @param line_size Line width for \code{geom_linerange}. Default 0.3.
#' @param show_stimulus Logical; if TRUE, try to compose a stimulus panel underneath. Default TRUE.
#' @param stimulus_height Relative height of the stimulus panel in \code{cowplot::plot_grid}. Default 0.6.
#' @param row_labels Optional character vector of row labels for Case A; if provided,
#'   this overrides \code{label_col}.
#' @param label_col Optional character scalar naming a Metadata column to use as
#'   row labels in Case A when \code{row_labels} is \code{NULL}. Ignored in Case B.
#' @param step_gap Non-negative numeric. If \code{> 0} (Case A only), inserts a
#'   vertical gap of this size (in y-units) between successive label groups and
#'   switches the y-axis to group-midpoint ticks (block labels).
#' @param step_order_decreasing Logical. If TRUE (Case A only), reverses the order
#'   of label groups. Ordering is numeric if labels parse as numeric, otherwise
#'   lexicographic.
#' @param background NULL (default) or the name of a Metadata column. If non-NULL,
#'   a background tile is drawn for each raster row, with \code{alpha} mapped to that column.
#' @param background_range Optional numeric length-2 range used to draw the background tile
#'   along x (typically the stimulus-encoding rectangle window).
#' @param background_fill The fill colour for the background tiles.
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
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @importFrom ggpubr theme_pubr
#' @importFrom ggplot2 ggplot aes geom_blank geom_tile geom_linerange
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
                   background_fill     = "black") {

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
            y_axis_title <- if (caseA) "Trace (Metadata row)" else "Channel"

            if (caseA) {
              rows_use <- seq_len(n_rows)

              if (!is.null(row_labels)) {
                stopifnot(length(row_labels) == length(rows_use))
                rowlabs <- setNames(as.character(row_labels), as.character(rows_use))

              } else if (!is.null(label_col) && is.character(label_col) && length(label_col) == 1L &&
                         !is.null(md) && label_col %in% names(md)) {

                rowlabs <- setNames(as.character(md[[label_col]][rows_use]), as.character(rows_use))
                if (step_gap > 0) y_axis_title <- label_col

              } else {
                if (!is.null(md)) {
                  if (!is.null(rownames(md)) && length(rownames(md))) {
                    rowlabs <- setNames(rownames(md)[rows_use], as.character(rows_use))
                  } else if ("RecordingID" %in% names(md)) {
                    rowlabs <- setNames(as.character(md$RecordingID[rows_use]), as.character(rows_use))
                  } else {
                    rowlabs <- setNames(paste0("row_", rows_use), as.character(rows_use))
                  }
                } else {
                  rowlabs <- setNames(paste0("row_", rows_use), as.character(rows_use))
                }
              }

              # attach label per spike row
              df$label <- unname(rowlabs[as.character(df$unit)])
              df$label[is.na(df$label) | !nzchar(df$label)] <- "NA"

            } else {
              rowlabs <- setNames(unique(df$unit), unique(df$unit)) # channel names as labels
            }

            # ---- Map rows to y positions (Case A: order/gap by label; Case B unchanged)
            if (caseA) {
              rows_use <- sort(unique(df$unit))
              lab_by_row <- unname(rowlabs[as.character(rows_use)])
              lab_by_row[is.na(lab_by_row) | !nzchar(lab_by_row)] <- "NA"

              # label levels in desired order (numeric if possible)
              lev <- unique(lab_by_row)
              lev_num <- suppressWarnings(as.numeric(lev))
              if (all(!is.na(lev_num))) {
                lev <- lev[order(lev_num, decreasing = isTRUE(step_order_decreasing))]
              } else {
                lev <- lev[order(lev, decreasing = isTRUE(step_order_decreasing))]
              }

              lab_rank <- match(lab_by_row, lev)

              step_key <- if (!is.null(md) && "Step" %in% names(md)) md$Step[rows_use] else rows_use
              run_key  <- if (!is.null(md) && "RunUID" %in% names(md)) as.character(md$RunUID[rows_use]) else as.character(rows_use)

              rows_ord <- rows_use[order(lab_rank, step_key, run_key, rows_use)]
              lab_ord  <- unname(rowlabs[as.character(rows_ord)])
              lab_ord[is.na(lab_ord) | !nzchar(lab_ord)] <- "NA"

              grp_num <- match(lab_ord, lev)

              y_pos <- seq_along(rows_ord)
              if (step_gap > 0) y_pos <- y_pos + (grp_num - 1) * step_gap

              y_map <- setNames(y_pos, as.character(rows_ord))
              df$y  <- unname(y_map[as.character(df$unit)])

              order_keys <- rows_ord

              # y-scale: per-row labels if no gaps; block midpoints if gaps
              if (step_gap > 0) {
                y_for_rows <- unname(y_map[as.character(order_keys)])
                grp_min <- tapply(y_for_rows, lab_ord, min)
                grp_max <- tapply(y_for_rows, lab_ord, max)
                y_breaks <- (grp_min + grp_max) / 2
                y_breaks <- unname(y_breaks[lev])

                y_labels <- lev
              } else {
                y_breaks <- unname(y_map[as.character(order_keys)])
                y_labels <- unname(rowlabs[as.character(order_keys)])
              }

            } else {
              order_keys <- unique(df$unit)
              y_map <- setNames(seq_along(order_keys), order_keys)
              df$y  <- unname(y_map[as.character(df$unit)])
              y_breaks <- seq_along(order_keys)
              y_labels <- unname(rowlabs)
            }

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
                linewidth = line_size, na.rm = TRUE
              ) +
              scale_y_continuous(
                breaks = y_breaks,
                labels = y_labels,
                expand = expansion(mult = c(0.02, 0.02))
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
              theme(
                panel.grid.major.y = element_blank(),
                panel.grid.minor   = element_blank()
              )

            p_out <- p_raster

            # attach tidy data
            attr(p_out, "raster_data") <- df

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
                  attr(p_out, "raster_data") <- df
                } else {
                  p_out <- p_raster
                  attr(p_out, "raster_data") <- df
                }
              }
            }

            p_out
          })
