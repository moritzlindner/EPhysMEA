#' Preprocess stimulus and compute phase by level crossings
#'
#' Upsamples a stimulus trace by linear interpolation and computes
#' instantaneous phase using level crossings. Phase advances linearly
#' from \code{start_deg} to \code{start_deg + 360} between consecutive
#' crossings of the stimulus relative to a chosen level.
#'
#' @param t Numeric vector of time points (strictly increasing).
#' @param stim Numeric vector of stimulus values (same length as \code{t}).
#' @param upsample Integer scalar. Factor by which to upsample the time
#'   series using linear interpolation (default \code{1}, no upsampling).
#' @param level Numeric scalar giving the crossing level. If \code{NULL}
#'   (default), the median of the stimulus is used.
#' @param direction Character string indicating which crossings to use:
#'   \code{"up"} (default), \code{"down"}, or \code{"both"}.
#' @param start_deg Numeric scalar. Phase angle (in degrees) assigned at
#'   the crossing points (default \code{0}).
#' @param extrapolate Logical; if \code{TRUE} (default), phase is also
#'   extrapolated before the first and after the last detected crossing
#'   using the first and last period lengths.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{t}{Upsampled time vector.}
#'   \item{stim}{Upsampled stimulus values.}
#'   \item{phase_deg}{Phase in degrees, wrapped to [0, 360).}
#' }
#'
#' @details
#' The function detects stimulus level crossings by linear interpolation
#' and assigns phase values that increase linearly from one crossing to
#' the next. Each full cycle spans 360 degrees. Upsampling improves
#' precision of crossing detection if the original sampling rate is low.
#'
#' @examples
#' # Example: sine stimulus at coarse sampling
#' t <- seq(0, 1, length.out = 50)
#' stim <- sin(2 * pi * 5 * t)
#' df_phase <- phasemetrics_stimulus_preprocess(t, stim, upsample = 10,
#'                                       level = 0, direction = "up")
#' head(df_phase)
#'
#' @keywords internal
phasemetrics_stimulus_preprocess <- function(
    t, stim,
    upsample = 1L,              # e.g. 10 = 10× upsampling
    level = NULL,               # crossing level; NULL = median(stim)
    direction = c("up","down","both"),
    start_deg = 0,
    extrapolate = TRUE
) {
  stopifnot(length(t) == length(stim), is.numeric(t), is.numeric(stim))
  if (is.unsorted(t, strictly = TRUE)) stop("'t' must be strictly increasing.")
  direction <- match.arg(direction)

  # --- upsample with linear interpolation ---
  upsample <- as.integer(upsample)
  if (upsample < 1L) upsample <- 1L
  if (upsample == 1L) {
    t_up <- to_num(t)
    stim_up <- as.numeric(stim)
  } else {
    n_new <- length(t) * upsample
    t_up <- seq(from = t[1], to = t[length(t)], length.out = n_new)
    stim_up <- approx(x = t, y = stim, xout = t_up, method = "linear",
                      rule = 2, ties = "ordered")$y
  }

  # --- detect level crossings (linear interpolation) ---
  eps_zero <- .Machine$double.eps^0.5
  if (is.null(level)) level <- stats::median(stim_up, na.rm = TRUE)
  x <- stim_up - level
  x[abs(x) < eps_zero] <- 0

  up_idx   <- which(x[-length(x)] <  0 & x[-1L] >= 0)
  down_idx <- which(x[-length(x)] >  0 & x[-1L] <= 0)
  cand <- switch(direction,
                 up   = up_idx,
                 down = down_idx,
                 both = sort(c(up_idx, down_idx)))

  tc <- numeric(0)
  if (length(cand)) {
    t0 <- t_up[cand]; t1 <- t_up[cand+1L]
    x0 <- x[cand]; x1 <- x[cand+1L]
    a  <- -x0 / (x1 - x0)
    tc <- sort(t0 + a * (t1 - t0))
  }
  nZ <- length(tc)

  # --- assign phases ---
  phase_deg <- rep(NA_real_, length(t_up))
  if (nZ >= 2L) {
    k <- findInterval(t_up, tc, left.open = FALSE, rightmost.closed = FALSE)
    inside <- (k >= 1 & k < nZ)
    if (any(inside)) {
      ti <- t_up[inside]; kI <- k[inside]
      tL <- tc[kI]; tR <- tc[kI+1L]
      frac <- (ti - tL) / (tR - tL)
      phase_deg[inside] <- start_deg + 360 * (kI - 1) + 360 * frac
    }
    if (extrapolate) {
      left <- which(!inside & t_up < tc[1])
      if (length(left)) {
        T1 <- tc[2] - tc[1]; fracL <- (t_up[left] - tc[1]) / T1
        phase_deg[left] <- start_deg - 360 + 360 * (1 + fracL)
      }
      right <- which(!inside & t_up >= tc[nZ])
      if (length(right)) {
        Tn <- tc[nZ] - tc[nZ-1]; fracR <- (t_up[right] - tc[nZ]) / Tn
        phase_deg[right] <- start_deg + 360 * (nZ - 1) + 360 * fracR
      }
    }
  }

  # wrap to [0,360)
  phase_deg <- ((phase_deg %% 360) + 360) %% 360

  data.frame(
    t = t_up,
    stim = stim_up,
    phase_deg = phase_deg
  )
}
