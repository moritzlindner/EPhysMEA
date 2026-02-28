#' Draw a primitive PSTH heatmap from binned counts/rates
#'
#' Creates a ggplot2 heatmap where rows correspond to units and columns correspond
#' to consecutive time bins. The function assumes `bin_mtx` and `meta` are already
#' row-aligned (same unit order), and focuses on plotting with optional unit
#' sorting and per-unit (row-wise) normalization.
#'
#' If `sort_by` is provided, the y-axis variable is derived from the corresponding
#' metadata column(s): a single column is used directly; multiple columns are
#' combined via `interaction()`. If resulting labels are not unique, they are
#' made unique (to avoid overplotting) while preserving order.
#'
#' @param bin_mtx Matrix-like object (matrix/data.frame) of size
#'   `n_units x n_bins`. Rows are units; columns are consecutive time bins
#'   in their existing order.
#' @param meta A `data.frame` with `n_units` rows containing metadata aligned
#'   to `bin_mtx` by row order (i.e., `meta[i, ]` describes `bin_mtx[i, ]`).
#' @param normalize Logical; if `TRUE`, applies per-row min-max scaling to [0, 1]
#'   before plotting. For constant rows (max == min) the normalized row is set to 0.
#' @param sort_by Optional character vector specifying one or more metadata columns
#'   (in `meta`) to sort units before plotting. Use `"-colname"` to sort a column
#'   in decreasing order. Sorting is applied left-to-right in `sort_by` (primary,
#'   secondary, ...).
#' @param bin_dur Single positive numeric scalar giving the bin duration in seconds.
#'   Used to construct the x-axis time coordinate as `time_s = (bin_index - 1) * bin_dur`.
#' @param base_size Base font size passed to `ggpubr::theme_pubr()`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The y-axis tick labels are hidden by default (to accommodate large numbers
#' of units). The first row of `bin_mtx` is plotted at the top of the heatmap.
#'
#' @examples
#' \dontrun{
#' p1 <- ggPSTHHeatmap.primitive(bin_mtx, meta, bin_dur = 0.02, base_size = 8)
#' p2 <- ggPSTHHeatmap.primitive(bin_mtx, meta, bin_dur = 0.02,
#'                              normalize = TRUE,
#'                              sort_by = c("unit_order_cond", "-dist_to_cond"),
#'                              base_size = 8)
#' print(p2)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_raster labs theme element_blank scale_fill_viridis_c
#' @importFrom ggpubr theme_pubr
#' @export
ggPSTHHeatmap.primitive <- function(bin_mtx, meta, normalize = FALSE,
                                    sort_by = NULL, bin_dur, base_size = 8) {
  bin_mtx <- as.matrix(bin_mtx)

  if (!is.data.frame(meta)) meta <- as.data.frame(meta)
  if (nrow(meta) != nrow(bin_mtx)) {
    stop("meta must have the same number of rows as bin_mtx (row-order matched).")
  }
  if (!is.numeric(bin_dur) || length(bin_dur) != 1L || !is.finite(bin_dur) || bin_dur <= 0) {
    stop("bin_dur must be a single, positive numeric value (seconds).")
  }

  # ---- optional sorting by one or more metadata columns
  if (!is.null(sort_by)) {
    if (!is.character(sort_by) || length(sort_by) < 1L) stop("sort_by must be a character vector.")

    keys <- vector("list", length(sort_by))
    for (k in seq_along(sort_by)) {
      nm  <- sort_by[k]
      dec <- startsWith(nm, "-")
      col <- sub("^-", "", nm)

      if (!col %in% names(meta)) stop("sort_by column not found in meta: ", col)

      x <- meta[[col]]
      key <- if (is.numeric(x)) x else xtfrm(x)
      if (dec) key <- -key
      keys[[k]] <- key
    }

    ord <- do.call(order, c(keys, list(na.last = TRUE)))
    meta <- meta[ord, , drop = FALSE]
    bin_mtx <- bin_mtx[ord, , drop = FALSE]
  }

  # ---- by-row 0–1 normalization
  if (isTRUE(normalize)) {
    rmin <- apply(bin_mtx, 1, min, na.rm = TRUE)
    rmax <- apply(bin_mtx, 1, max, na.rm = TRUE)
    denom <- rmax - rmin
    denom[!is.finite(denom) | denom == 0] <- NA_real_
    bin_mtx <- (bin_mtx - rmin) / denom
    bin_mtx[is.na(bin_mtx)] <- 0
  }

  # ---- y variable derived from sort_by (or row index if absent)
  nU <- nrow(bin_mtx)
  sort_cols <- if (!is.null(sort_by)) sub("^-", "", sort_by) else character(0)

  if (length(sort_cols) == 1L) {
    y_lab <- as.character(meta[[sort_cols]])
  } else if (length(sort_cols) > 1L) {
    y_lab <- as.character(interaction(meta[, sort_cols, drop = FALSE], drop = TRUE))
  } else {
    y_lab <- as.character(seq_len(nU))
  }

  # ensure y labels are unique to prevent row overplotting
  if (anyDuplicated(y_lab)) y_lab <- make.unique(y_lab, sep = " | ")

  # ---- long format by position
  nB <- ncol(bin_mtx)
  df <- data.frame(
    unit      = rep(seq_len(nU), each = nB),
    bin_index = rep(seq_len(nB), times = nU),
    value     = as.vector(t(bin_mtx))
  )
  df$time_s <- (df$bin_index - 1L) * bin_dur

  # map unit index -> label, preserve “first row at top”
  df$unit <- factor(y_lab[df$unit], levels = rev(y_lab))

  heat <- ggplot2::ggplot(df, ggplot2::aes(x = time_s, y = unit, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::labs(x = "Time (s)", y = NULL) +
    ggpubr::theme_pubr(base_size = base_size) #+
    # ggplot2::theme(
    #   axis.text.y  = ggplot2::element_blank(),
    #   axis.ticks.y = ggplot2::element_blank()
    #)

  if (isTRUE(normalize)) {
    heat <- heat + ggplot2::scale_fill_viridis_c(name = "Normalized Rate", limits = c(0, 1))
  } else {
    heat <- heat + ggplot2::scale_fill_viridis_c(name = "Rate (1/s)")
  }

  heat
}
