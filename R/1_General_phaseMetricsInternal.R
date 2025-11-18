#' Map spike times to instantaneous stimulus phase (degrees)
#'
#' Linearly interpolates a phase trace to spike times.
#' Values at the edges are extended (rule = 2).
#'
#' @param spikes Numeric vector of spike times.
#' @param t Numeric vector of sample times (strictly increasing).
#' @param phase_trace Numeric vector of phase values (degrees) at times \code{t}.
#'
#' @return A numeric vector of phases (degrees), one per spike in \code{spikes}.
#' @details Uses \code{approx(..., rule = 2, ties = "ordered")}.
#' @keywords internal
spike_phase_from_stimulus <- function(spikes, t, phase_trace) {
  approx(x = t, y = phase_trace, xout = spikes, rule = 2, ties = "ordered")$y
}


#' Rayleigh strength within a peak-centered window (per segment)
#'
#' For each segment, builds a circular phase histogram (bin width \code{bin_width_deg})
#' from spike phases (degrees), finds the peak bin, and computes the Rayleigh
#' strength statistic
#' \deqn{Z = n R^2, \quad R = \left|\frac{1}{n}\sum_i e^{i \theta_i}\right|}
#' using only spikes whose phase lies within \eqn{\pm} \code{window_half_deg}
#' degrees of the histogram peak (circularly).
#'
#' @param df A \code{data.frame} containing at least \code{phase_col} and \code{seg_col}.
#' @param phase_col Column name with spike phase in degrees (default \code{"p"}).
#' @param seg_col Column name with segment labels (default \code{"seg"}).
#' @param bin_width_deg Histogram bin width in degrees (default \code{15}).
#' @param window_half_deg Half-width of the peak-centered window in degrees (default \code{90}).
#'
#' @return A named numeric vector of length equal to the number of segments.
#'   Names are the segment labels; values are Rayleigh strength (\eqn{Z}) within
#'   the peak window. Segments with no spikes (or no spikes in-window) yield
#'   \code{NA_real_}.
#'
#' @details
#' The Rayleigh strength increases both with the sharpness of phase locking
#' and with the number of spikes. It is the test statistic used in the
#' Rayleigh test for non-uniformity of circular data.
#'
#' @seealso \link[circular]{rayleigh.test}
#' @keywords internal
rayleigh_peak_by_segment <- function(df,
                                     phase_col = "p",
                                     seg_col   = "seg",
                                     bin_width_deg   = 15,
                                     window_half_deg = 45) {
  stopifnot(is.data.frame(df),
            phase_col %in% names(df),
            seg_col   %in% names(df))
  p_all   <- df[[phase_col]]
  seg_all <- df[[seg_col]]

  # wrap degrees to [0, 360)
  wrap360 <- function(x) (x %% 360 + 360) %% 360
  circ_abs_diff <- function(a, b) abs(((a - b + 180) %% 360) - 180)

  seg_levels <- sort(unique(seg_all))
  out <- setNames(rep(NA_real_, length(seg_levels)), as.character(seg_levels))

  breaks <- seq(0, 360, by = bin_width_deg)

  for (s in seg_levels) {
    idx <- which(seg_all == s)
    if (!length(idx)) next

    p_deg <- wrap360(p_all[idx])
    if (!length(p_deg)) { out[[as.character(s)]] <- 0; next }

    # Histogram to find the peak (use left-closed intervals [a,b) with right=FALSE)
    h <- hist(p_deg, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
    counts <- h$counts
    if (all(counts == 0)) { out[[as.character(s)]] <- 0; next }

    peak_center <- h$mids[which.max(counts)]

    # spikes inside the ±window around the peak (circularly)
    in_win <- circ_abs_diff(p_deg, peak_center) <= (window_half_deg + 1e-12)
    if (!any(in_win)) { out[[as.character(s)]] <- 0; next }

    th <- p_deg[in_win] * pi / 180
    n  <- length(th)
    R  <- sqrt((sum(cos(th)))^2 + (sum(sin(th)))^2) / n
    Z  <- n * R * R
    out[[as.character(s)]] <- as.numeric(Z)
  }

  out
}

#' Phase-opponent modulation depth per segment (trial-flattened)
#'
#' Computes simple, phase-based response amplitudes per segment by pooling
#' spikes across trials and operating on spike phases (degrees).
#'
#' Supported \code{methods} (any subset; column per method in return):
#' \itemize{
#'   \item \code{"peak_trough"}: Build a circular histogram with
#'         width \code{bin_width_deg}; return \code{max(counts) - min(counts)}.
#'   \item \code{"opponent180"}: Find peak bin center; compute counts in a
#'         \eqn{\pm} \code{window_half_deg} window around the peak and around the
#'         \code{peak+180°} opponent; return \code{peak - opponent}.
#'   \item \code{"opponent90"}: As above, but opponent is the better of the two
#'         quadrature centers (\code{peak+90°} or \code{peak-90°}); return
#'         \code{peak - min(opp+90, opp-90)}.
#'   \item \code{"opponent_best"}: \code{max(opponent180, opponent90)}.
#'         Attribute \code{"opponent_best_winner"} (named vector) stores
#'         \code{180} or \code{90} per segment for traceability.
#'   \item \code{"isi50"}: Order phases on the circle; compute circular ISIs
#'         (arc lengths). Let \eqn{m = \lfloor n/2 \rfloor}. Slide a window of
#'         \eqn{m} consecutive ISIs to find the shortest arc covering half the
#'         spikes (\eqn{L_\mathrm{min}} degrees). Return
#'         \code{180 - L_min} (positive = dense half fits in < 180°).
#' }
#'
#' Notes:
#' \itemize{
#'   \item Trials must be flattened before calling (pass all spikes for a segment).
#'   \item Outputs are \emph{differences in counts} (or degrees for \code{isi50}).
#'         If you prefer rates, divide by window width or total spikes upstream.
#'   \item Segments with too few spikes yield \code{NA}.
#' }
#'
#' @param df A \code{data.frame} with a phase column (degrees) and a segment label.
#' @param phase_col Name of the phase column in degrees (default \code{"p"}).
#' @param seg_col Name of the segment column (default \code{"seg"}).
#' @param bin_width_deg Histogram bin width (default \code{15}).
#' @param window_half_deg Half-width for opponent windows (default \code{45}).
#' @param methods Character vector; any of
#'   \code{c("peak_trough","opponent180","opponent90","opponent_best","isi50")}.
#'
#' @return A numeric matrix with rows = segments (in ascending numeric order) and
#'   columns = requested \code{methods}. For \code{"opponent_best"}, an attribute
#'   \code{"opponent_best_winner"} (named integer vector of 180/90) is attached.
#' @keywords internal
phasemetrics_opponent_modulation_depth_by_segment <- function(
    df,
    phase_col       = "p",
    seg_col         = "seg",
    bin_width_deg   = 15,
    window_half_deg = 45,
    methods = c("peak_trough","opponent180","opponent90","opponent_best","isi50")
) {
  stopifnot(is.data.frame(df),
            phase_col %in% names(df),
            seg_col   %in% names(df))
  methods <- unique(methods)
  allowed <- c("peak_trough","opponent180","opponent90","opponent_best","isi50")
  if (!all(methods %in% allowed)) {
    stop("Unknown method(s): ", paste(setdiff(methods, allowed), collapse = ", "))
  }

  # --- helpers ---------------------------------------------------------------
  wrap360 <- function(x) ((x %% 360) + 360) %% 360
  circ_abs_diff <- function(a, b) abs(((a - b + 180) %% 360) - 180)

  seg_levels <- sort(unique(as.numeric(df[[seg_col]])))
  out <- matrix(NA_real_, nrow = length(seg_levels), ncol = length(methods),
                dimnames = list(as.character(seg_levels), methods))

  # histogram edges & mids for peak finding
  edges <- seq(0, 360, by = bin_width_deg)
  if (tail(edges, 1) < 360) edges <- c(edges, 360)
  mids  <- (head(edges, -1) + tail(edges, -1)) / 2

  # store winner angle for "opponent_best"
  winner <- setNames(rep(NA_integer_, length(seg_levels)), as.character(seg_levels))

  # --- per segment -----------------------------------------------------------
  for (i in seq_along(seg_levels)) {
    s <- seg_levels[i]
    p_all <- df[[phase_col]][df[[seg_col]] == s]
    if (!length(p_all)) next
    p_deg <- wrap360(p_all)
    nspk  <- length(p_deg)

    # require >= 1 spike for peak_trough/opponent; >=2 for isi50
    # histogram for peak finding
    h <- hist(p_deg, breaks = edges, plot = FALSE, include.lowest = TRUE, right = FALSE)
    counts <- h$counts
    if (!all(counts == 0)) {
      peak_center <- mids[which.max(counts)]
    } else {
      peak_center <- NA_real_
    }

    # compute window counts at given centers
    win_count <- function(center_deg) {
      if (!is.finite(center_deg)) return(NA_real_)
      sum(circ_abs_diff(p_deg, center_deg) <= (window_half_deg + 1e-12))
    }

    # --- peak_trough
    if ("peak_trough" %in% methods) {
      out[i, "peak_trough"] <- if (any(counts > 0)) max(counts) - min(counts) else NA_real_
    }

    # --- opponent180 / opponent90 / opponent_best
    opp180_val <- NA_real_
    opp90_val  <- NA_real_
    if (is.finite(peak_center)) {
      peak_ct <- win_count(peak_center)
      # 180°
      opp180_ct <- win_count(peak_center + 180)
      if (!is.na(peak_ct) && !is.na(opp180_ct)) {
        opp180_val <- peak_ct - opp180_ct
      }
      # 90° (take the more opponent of +90 / -90, i.e., smaller count → larger difference)
      opp90_ct1 <- win_count(peak_center + 90)
      opp90_ct2 <- win_count(peak_center - 90)
      opp90_ct  <- suppressWarnings(min(opp90_ct1, opp90_ct2, na.rm = TRUE))
      if (!is.na(peak_ct) && is.finite(opp90_ct)) {
        opp90_val <- peak_ct - opp90_ct
      }
    }
    if ("opponent180" %in% methods) out[i, "opponent180"] <- opp180_val
    if ("opponent90"  %in% methods) out[i, "opponent90"]  <- opp90_val
    if ("opponent_best" %in% methods) {
      best <- suppressWarnings(max(opp180_val, opp90_val, na.rm = TRUE))
      out[i, "opponent_best"] <- if (is.finite(best)) best else NA_real_
      if (is.finite(best)) {
        winner[i] <- if (!is.na(opp180_val) && best == opp180_val) 180L else 90L
      }
    }

    # --- isi50 (densest half-arc vs 180°)
    if ("isi50" %in% methods) {
      if (nspk >= 2L) {
        ps <- sort(p_deg)
        # circular ISIs (degrees)
        diffs <- c(diff(ps), (ps[1] + 360) - ps[length(ps)])
        m <- floor(nspk / 2)  # number of consecutive ISIs to sum
        if (m >= 1L) {
          # sliding sum of m consecutive ISIs on the ring
          diffs2 <- c(diffs, diffs[1:(m-1)])  # extend for wrap
          cs <- cumsum(diffs2)
          # window sums starting at each index 1..nspk
          win_sums <- cs[seq(m, length.out = nspk)] - c(0, cs[seq_len(nspk-1)])
          Lmin <- min(win_sums, na.rm = TRUE)
          out[i, "isi50"] <- 180 - Lmin
        } else {
          out[i, "isi50"] <- NA_real_
        }
      } else {
        out[i, "isi50"] <- NA_real_
      }
    }
  }

  if ("opponent_best" %in% methods) {
    attr(out, "opponent_best_winner") <- winner
  }
  out
}
