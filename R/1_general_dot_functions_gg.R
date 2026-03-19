# ---- Internal helper shared by ggEventRaster() and ggPSTHHeatmap() -----------------

#' Auto-layout rows for EPhys* plots (shared helper)
#'
#' Internal helper that implements the shared “auto-layout” logic used by
#' \code{\link{ggEventRaster}} and \code{\link{ggPSTHHeatmap}}:
#' either one channel across many Metadata rows (Case A), or many channels within
#' a single Metadata row (Case B). It also implements label grouping + gaps via
#' \code{label_col} and \code{step_gap}, mirroring ggEventRaster.
#'
#' @param X An \code{EPhysContainer} (e.g. \code{EPhysEvents}, \code{EPhysContinuous}).
#' @param row_labels Optional explicit labels for Case A (length = nrow(Metadata(X))).
#' @param label_col Optional Metadata column name used as row label in Case A.
#' @param step_gap Non-negative numeric; if > 0 inserts gaps between label groups (Case A).
#' @param step_order_decreasing Logical; reverse label-group order (Case A).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{caseA}, \code{caseB}
#'   \item \code{rowlabs}: named character vector mapping unit id -> label
#'   \item \code{order_keys}: unit ids in plotting order
#'   \item \code{y_map}: named numeric vector mapping unit id -> y position
#'   \item \code{y_breaks}, \code{y_labels}, \code{y_axis_title}
#' }
#'
#' @keywords internal
#' @noRd
#' @importFrom EPhysData Metadata Channels
.ephys_autolayout <- function(X,
                              row_labels = NULL,
                              label_col = NULL,
                              step_gap = 0,
                              step_order_decreasing = FALSE) {

  stopifnot(is.numeric(step_gap), length(step_gap) == 1L, is.finite(step_gap), step_gap >= 0)

  md  <- EPhysData::Metadata(X)
  chs <- EPhysData::Channels(X)
  y_tick_breaks <- NULL

  n_rows <- if (!is.null(md)) nrow(md) else NA_integer_

  caseA <- !is.null(chs) && length(chs) == 1L              # one channel across rows
  caseB <- !is.null(md)  && is.finite(n_rows) && n_rows == 1L  # one row with many channels

  if (!(caseA || caseB)) {
    stop(
      "Object must satisfy one of:\n",
      "  • length(Channels(X)) == 1, or\n",
      "  • nrow(Metadata(X)) == 1."
    )
  }

  # ---- Labels
  if (caseA) {
    rows_use <- seq_len(n_rows)

    if (!is.null(row_labels)) {
      stopifnot(length(row_labels) == length(rows_use))
      rowlabs <- setNames(as.character(row_labels), as.character(rows_use))

    } else if (!is.null(label_col) &&
               is.character(label_col) && length(label_col) == 1L &&
               !is.null(md) && label_col %in% names(md)) {

      rowlabs <- setNames(as.character(md[[label_col]][rows_use]), as.character(rows_use))

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

    y_axis_title <- "Trace (Metadata row)"
    if (step_gap > 0 && is.null(row_labels) && !is.null(label_col) && label_col %in% names(md)) {
      y_axis_title <- label_col
    }

  } else {
    # Case B: channel names as labels
    ch_use <- as.character(chs)
    rowlabs <- setNames(ch_use, ch_use)
    y_axis_title <- "Channel"
  }

  # ---- Order + y positions (gaps by label in Case A)
  if (caseA) {
    rows_use <- as.integer(names(rowlabs))
    lab_by_row <- unname(rowlabs[as.character(rows_use)])
    lab_by_row[is.na(lab_by_row) | !nzchar(lab_by_row)] <- "NA"

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
    order_keys <- rows_ord

    if (step_gap > 0) {
      y_for_rows <- unname(y_map[as.character(order_keys)])
      grp_min <- tapply(y_for_rows, lab_ord, min)
      grp_max <- tapply(y_for_rows, lab_ord, max)
      y_breaks <- (grp_min + grp_max) / 2
      y_breaks <- unname(y_breaks[lev])
      y_labels <- lev
      if (length(lev) >= 2L) {
        y_tick_breaks <- (unname(grp_max[lev[-length(lev)]]) + unname(grp_min[lev[-1L]])) / 2
      }
    } else {
      y_breaks <- unname(y_map[as.character(order_keys)])
      y_labels <- unname(rowlabs[as.character(order_keys)])
    }

  } else {
    order_keys <- names(rowlabs)
    y_map <- setNames(seq_along(order_keys), order_keys)
    y_breaks <- unname(y_map)
    y_labels <- unname(rowlabs)
  }

  list(
    caseA = caseA, caseB = caseB,
    rowlabs = rowlabs,
    order_keys = order_keys,
    y_map = y_map,
    y_breaks = y_breaks,
    y_labels = y_labels,
    y_tick_breaks = y_tick_breaks,
    y_axis_title = y_axis_title
  )
}


#' Parse metadata sort specification
#'
#' Internal helper that parses a character \code{sort_by} specification of the
#' form used by plotting helpers such as \code{\link{ggEventRaster}} and
#' \code{\link{ggPSTHHeatmap.primitive}}.
#'
#' Entries in \code{sort_by} may optionally be prefixed with \code{"-"} to request
#' descending order for that metadata column. The helper strips these prefixes,
#' validates the resulting column names against \code{meta_names}, optionally
#' removes excluded columns, and returns the parsed column names together with a
#' logical vector indicating descending order.
#'
#' @param sort_by \code{NULL} or character vector naming metadata columns to sort
#'   by. A leading \code{"-"} indicates reverse sorting for that column.
#' @param meta_names Character vector of valid metadata column names, typically
#'   \code{names(meta)}.
#' @param exclude Optional character vector of metadata columns to remove from the
#'   parsed result after stripping any leading \code{"-"}. This is useful, for
#'   example, when a faceting column should not also participate in row ordering.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{cols}}{Character vector of parsed metadata column names.}
#'   \item{\code{decreasing}}{Logical vector of the same length as \code{cols};
#'   \code{TRUE} indicates descending sort order.}
#' }
#'
#' @seealso \code{\link{ggPSTHHeatmap.primitive}}, \code{\link{ggEventRaster}}
#' @keywords internal
.parse_sort_by_spec <- function(sort_by, meta_names, exclude = NULL) {
  if (is.null(sort_by) || !length(sort_by)) {
    return(list(
      cols = character(0),
      decreasing = logical(0)
    ))
  }

  stopifnot(is.character(sort_by))

  decreasing <- startsWith(sort_by, "-")
  cols <- sub("^-", "", sort_by)

  bad <- !nzchar(cols)
  if (any(bad)) {
    stop("Invalid sort_by entry/entries: ", paste(sort_by[bad], collapse = ", "))
  }

  if (!is.null(exclude) && length(exclude)) {
    keep <- !(cols %in% exclude)
    cols <- cols[keep]
    decreasing <- decreasing[keep]
  }

  miss <- setdiff(cols, meta_names)
  if (length(miss)) {
    stop("sort_by columns not in meta: ", paste(miss, collapse = ", "))
  }

  list(
    cols = cols,
    decreasing = decreasing
  )
}


#' Order metadata rows from parsed sort specification
#'
#' Internal helper that computes a row permutation for a metadata
#' \code{\link[base:data.frame]{data.frame}} from a parsed sort specification as
#' returned by \code{\link{.parse_sort_by_spec}}.
#'
#' Rows can be ordered globally or within groups defined by \code{split_by}. This
#' is useful for plot helpers that need a stable metadata-driven row order, for
#' example \code{\link{ggPSTHHeatmap.primitive}} and
#' \code{\link{ggEventRaster}}.
#'
#' @param meta A \code{data.frame} whose rows are to be reordered.
#' @param sort_spec A named list as returned by \code{\link{.parse_sort_by_spec}},
#'   with components \code{cols} and \code{decreasing}.
#' @param split_by \code{NULL} or character scalar naming a column in \code{meta}.
#'   If supplied, rows are ordered within each group of \code{meta[[split_by]]},
#'   and the resulting group-wise permutations are concatenated. For factor
#'   columns, group order follows factor levels; otherwise it follows first
#'   appearance in \code{meta}. Missing values are placed last.
#'
#' @return Integer vector of row indices giving the requested order.
#'
#' @seealso \code{\link{.parse_sort_by_spec}}, \code{\link{ggPSTHHeatmap.primitive}},
#'   \code{\link{ggEventRaster}}
#' @keywords internal
.order_meta_rows <- function(meta, sort_spec, split_by = NULL) {
  stopifnot(is.data.frame(meta))
  stopifnot(is.list(sort_spec), all(c("cols", "decreasing") %in% names(sort_spec)))
  stopifnot(length(sort_spec$cols) == length(sort_spec$decreasing))

  if (length(sort_spec$cols)) {
    miss <- setdiff(sort_spec$cols, names(meta))
    if (length(miss)) {
      stop("sort_spec columns not in meta: ", paste(miss, collapse = ", "))
    }
  }

  if (!is.null(split_by)) {
    stopifnot(is.character(split_by), length(split_by) == 1L)
    if (!split_by %in% names(meta)) {
      stop("split_by column not in meta: ", split_by)
    }
  }

  order_local <- function(m) {
    if (!length(sort_spec$cols)) {
      return(seq_len(nrow(m)))
    }

    ord_args <- Map(function(col, dec) {
      xt <- xtfrm(m[[col]])
      if (isTRUE(dec)) -xt else xt
    }, sort_spec$cols, sort_spec$decreasing)

    do.call(order, c(ord_args, list(na.last = TRUE)))
  }

  if (is.null(split_by)) {
    return(order_local(meta))
  }

  fac <- meta[[split_by]]

  if (is.factor(fac)) {
    split_levels <- levels(fac)
    split_levels <- split_levels[split_levels %in% as.character(fac)]

    ord <- unlist(lapply(split_levels, function(lvl) {
      ii <- which(as.character(fac) == lvl)
      ii[order_local(meta[ii, , drop = FALSE])]
    }), use.names = FALSE)
  } else {
    split_vals <- unique(fac[!is.na(fac)])

    ord <- unlist(lapply(split_vals, function(val) {
      ii <- which(!is.na(fac) & fac == val)
      ii[order_local(meta[ii, , drop = FALSE])]
    }), use.names = FALSE)
  }

  if (anyNA(fac)) {
    ii_na <- which(is.na(fac))
    ord <- c(ord, ii_na[order_local(meta[ii_na, , drop = FALSE])])
  }

  ord
}
