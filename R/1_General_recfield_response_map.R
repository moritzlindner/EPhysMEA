#' Receptive-field response map (mean Response_Amplitude per pixel)
#'
#' Builds a 2-D receptive-field map for a single channel by averaging
#' \code{Response_Amplitude} (from \code{\link[=step_stim_resp_metrics]{step_stim_resp_metrics}})
#' over all frames in which each pixel was stimulated.
#'
#' @description
#' The stimulus sequence is provided as:
#' \itemize{
#'   \item \code{StimSeq[,,i]} — the spatial mask (binary/weight) for frame \eqn{i}.
#'   \item \code{change_times[i]} — the onset time (s) of frame \eqn{i}.
#' }
#' Each frame \eqn{i} is assumed to be shown \emph{until} \code{change_times[i+1]}
#' (the next onset). If \code{length(change_times) == dim(StimSeq)[3]}, the last
#' frame's offset is inferred as \code{change_times[n] + median(diff(change_times))}.
#'
#' For each frame, a scalar \emph{Response_Amplitude} is computed and assigned to
#' all pixels where \code{StimSeq[,,i] != 0}. If a pixel is stimulated multiple
#' times, values are averaged. Pixels never stimulated are returned as \code{NA}.
#'
#' @param spike_times Numeric vector of spike times (seconds), aligned so
#'   \code{t = 0} corresponds to the start of \code{StimSeq}. Can be unsorted.
#' @param StimSeq Numeric/integer/logical 3D array \code{[nrow, ncol, nframes]}.
#'   Each slice \code{StimSeq[,,i]} is the mask of pixels stimulated in frame \eqn{i}
#'   (nonzero values are treated as ON/weights).
#' @param change_times Numeric vector (seconds) of frame onset times.
#'   Must have length \code{nframes} or \code{nframes + 1}.
#' @param baserate_win Optional numeric scalar (s). If supplied, the baseline
#'   rate is computed once over \code{[change_times[1] - baserate_win, change_times[1])}
#'   and used for all frames. If \code{NULL} (default), the baseline from
#'   \code{step_stim_resp_metrics()} is used per frame.
#' @param ... Additional arguments passed to
#'   \code{\link[=step_stim_resp_metrics]{step_stim_resp_metrics}} (e.g., \code{trans_win},
#'   \code{binwidth}, \code{smooth}, \code{latency_*}).
#'
#' @return A numeric matrix \code{[nrow, ncol]} with the \emph{mean}
#'   \code{Response_Amplitude} per pixel. Pixels never stimulated are \code{NA}.
#' @importFrom stats median
#' @examples
#' \donttest{
#' # 4x4 grid, 8 frames, one pixel per frame, 0.2 s per frame
#' StimSeq <- array(0L, dim = c(4,4,8))
#' for (i in 1:8) StimSeq[((i-1) %% 4) + 1, ((i-1) %/% 4) + 1, i] <- 1L
#' change_times <- seq(0, by = 0.2, length.out = 8)  # onsets; last offset inferred
#' spikes <- c(0.05, 0.07, 0.22, 0.28, 0.43, 0.51, 0.71)
#' resp_map <- recfield_response_map(spikes, StimSeq, change_times,
#'                                   trans_win = 1.0, binwidth = 0.02)
#' image(resp_map, main = "Mean Response_Amplitude")
#' }
#'
#' @export
recfield_response_map <- function(
    spike_times,
    StimSeq,
    change_times,
    baserate_win = NULL,
    ...
) {
  # ---- Validate StimSeq ----
  if (length(dim(StimSeq)) != 3L) stop("StimSeq must be a 3D array [rows × cols × frames].")
  dims <- dim(StimSeq)
  nrow <- dims[1]; ncol <- dims[2]; nframes <- dims[3]

  # ---- Validate change_times and build offsets ----
  if (!is.numeric(change_times) || !length(change_times)) {
    stop("'change_times' must be a numeric vector of onsets (seconds).")
  }
  if (length(change_times) < nframes) {
    stop("length(change_times) must be nframes or nframes+1 (got ", length(change_times), ").")
  }
  step_onsets <- as.numeric(change_times[seq_len(nframes)])

  if (length(change_times) >= nframes + 1L) {
    step_offsets <- as.numeric(change_times[seq_len(nframes) + 1L])
  } else {
    diffs <- diff(step_onsets)
    if (!length(diffs) || !all(is.finite(diffs))) {
      stop("Cannot infer last frame duration from 'change_times'. Provide nframes+1 onsets.")
    }
    frame_dur <- median(diffs)
    step_offsets <- step_onsets + frame_dur
  }


  # Normalize and sort spikes once
  x <- sort(as.numeric(spike_times))
  x <- x[is.finite(x)]

  # ---- Optional global baseline rate (before first onset) ----
  global_base_rate <- NULL
  if (!is.null(baserate_win)) {
    stopifnot(is.numeric(baserate_win), length(baserate_win) == 1L, is.finite(baserate_win), baserate_win > 0)
    t0 <- step_onsets[1]
    base_times <- x[x >= (t0 - baserate_win) & x < t0]
    global_base_rate <- length(base_times) / baserate_win  # spikes/s
  }

  # ---- Preallocate accumulators (sum & count for internal mean) ----
  sum_map   <- matrix(0,  nrow = nrow, ncol = ncol)
  count_map <- matrix(0L, nrow = nrow, ncol = ncol)

  # ---- Main loop ----
  for (i in seq_len(nframes)) {
    m <- step_stim_resp_metrics(
      x,
      step_onset  = step_onsets[i],
      step_offset = step_offsets[i],
      latency_alpha = NA,
      shinomoto_win = NA,
      binwidth = NA,
      ...
    )
    rate_trans <- m$Response_Amplitude
    rate_base  <- if (is.null(global_base_rate)) m$Baseline_Rate else global_base_rate
    amp <- rate_trans - rate_base

    if (!is.finite(amp)) amp <- 0

    mask <- StimSeq[,,i]
    if (!is.null(dim(mask))) {
      nz <- (mask != 0)
      if (any(nz, na.rm = TRUE)) {
        # Weighted accumulation if mask has weights; mean is taken over exposures (not weight-normalized)
        sum_map[nz]   <- sum_map[nz]   + amp * mask[nz]
        count_map[nz] <- count_map[nz] + 1L
      }
    }
  }

  # ---- Compute mean per pixel; never-stimulated pixels -> NA ----
  response_map <- matrix(NA_real_, nrow = nrow, ncol = ncol)
  idx <- count_map > 0L
  # If masks can carry weights >1, you might prefer dividing by count_map or by total mask weight.
  # The MATLAB-style accumulation corresponds to dividing by count_map (exposures), as below:
  response_map[idx] <- sum_map[idx] / count_map[idx]

  response_map
}


#' Radially-flattened r50 and sector anisotropy from a 2D receptive-field map
#'
#' Computes a sub-pixel peak location and amplitude from a 2D receptive-field
#' (RF) map, estimates a radially-flattened \eqn{r_{50}} (radius at 50% of peak),
#' and quantifies anisotropy via sector-wise \eqn{r_{50}} (12 sectors, 30° each).
#' The pipeline optionally applies 2D wavelet denoising (MODWT) and Gaussian
#' smoothing before peak finding and profiling. A simple quality gate rejects
#' spurious single-pixel peaks using a peak signal-to-noise ratio (SNR) computed
#' from the outer ring of the map.
#'
#' @param rf Numeric matrix (2D receptive-field map). Values need not be
#'   normalized. Pixels may be finite everywhere; the function itself does
#'   not introduce \code{NA}s during smoothing.
#' @param pixel_width Numeric scalar: physical size of one pixel (e.g., in µm).
#'   Pixels are assumed square; \code{r50} is returned in the same units.
#' @param smooth_sd Numeric Gaussian sigma (in \emph{pixels}) for an additional
#'   isotropic blur (via \pkg{imager}). Set to \code{0} to disable the blur
#'   (default \code{0.7}).
#' @param SNR_min Numeric threshold for the peak quality gate. The peak must have
#'   \emph{both} (i) SNR \eqn{\ge} \code{SNR_min} relative to an outer-ring
#'   baseline and (ii) an absolute amplitude \eqn{\ge} \code{SNR_min} in the
#'   native units of \code{rf}. Default \code{2}.
#'
#' @details
#' \strong{Pre-processing.} The map is first denoised with
#' \code{\link[waveslim]{denoise.modwt.2d}} using \code{wf = "la8"} and
#' \code{J = 3} (shift-invariant). If \code{smooth_sd > 0}, a light Gaussian
#' blur is applied using \code{\link[imager]{isoblur}} with the given sigma in pixels.
#'
#' \strong{Sub-pixel peak.} A quadratic (3×3) fit around the integer-peak
#' provides sub-pixel row/column coordinates and the corresponding (smoothed)
#' peak amplitude.
#'
#' \strong{Radial flattening and \eqn{r_{50}}.} Values are converted to physical
#' coordinates using \code{pixel_width}. Around the sub-pixel peak, the map is
#' averaged in concentric annuli (bin width = \code{pixel_width}/2). The radial
#' profile is lightly smoothed, enforced to be non-increasing, and the 50% peak
#' crossing is found by linear interpolation to yield \code{r50}.
#'
#' \strong{Sector anisotropy.} The map is split into 12 angular sectors (30°).
#' Each sector’s local radial profile yields a sector-specific \eqn{r_{50}}.
#' The anisotropy index \code{symmetry_ratio} is defined as
#' \code{max(r50_sector) / r50_sector at +90°}. If the orthogonal sector has no
#' valid \eqn{r_{50}}, the ratio is \code{NA}.
#'
#' \strong{Quality gate (SNR).} To avoid returning \code{r50} for spurious
#' single-pixel blips, the peak must satisfy an SNR criterion computed from an
#' outer annulus: the baseline is \code{median(outer_ring)} and noise is
#' \code{1.4826 * mad(outer_ring)}. The outer ring is defined as the outer 30\%
#' of available radii. If either the SNR falls below \code{SNR_min} or the peak
#' amplitude is below \code{SNR_min}, the function returns \code{NA} for
#' \code{peak_row}, \code{peak_col}, \code{peak_amp}, \code{r50}, and
#' \code{symmetry_ratio}.
#'
#' @return A list with:
#' \describe{
#'   \item{peak_row}{Sub-pixel peak row coordinate (1-based, image indexing).}
#'   \item{peak_col}{Sub-pixel peak column coordinate (1-based).}
#'   \item{peak_amp}{Peak amplitude after pre-processing (denoise + optional blur).}
#'   \item{r50}{Radially-flattened 50\% radius (same units as \code{pixel_width}).}
#'   \item{symmetry_ratio}{\code{max(sector r50) / r50 at +90°}; \code{NA} if undefined.}
#'   \item{peak_at_margin}{Logical; \code{TRUE} if the (rounded) peak lies on the border.}
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item The Gaussian blur uses \code{imager::isoblur} with implicit
#'         reflect-like boundary handling; set \code{smooth_sd = 0} to skip it.
#'   \item The wavelet settings (\code{wf = "la8"}, \code{J = 3}) are reasonable
#'         defaults for 32–128 px RF maps; adjust \code{J} for larger fields.
#'   \item \code{r50} is undefined (\code{NA}) if the radial profile never
#'         crosses 50\% of the peak.
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("imager", quietly = TRUE) &&
#'     requireNamespace("waveslim", quietly = TRUE)) {
#'   set.seed(1)
#'   rf <- matrix(0, 48, 48)
#'   xx <- matrix(rep(seq_len(48)-0.5, each=48), 48, 48)
#'   yy <- matrix(rep(seq_len(48)-0.5, times=48), 48, 48)
#'   # Anisotropic Gaussian blob (no NAs needed)
#'   rf <- exp(-((xx-24)^2/(2*3^2) + (yy-30)^2/(2*6^2)))
#'   out <- recfield_r50_radial(rf, pixel_width = 10, smooth_sd = 0.7, SNR_min = 2)
#'   str(out)
#' }
#' }
#'
#' @seealso \code{\link[waveslim]{denoise.modwt.2d}},
#'   \code{\link[imager]{isoblur}}
#'
#' @importFrom imager as.cimg isoblur
#' @importFrom waveslim denoise.modwt.2d
#' @importFrom stats filter lm median mad coef
#' @export
recfield_r50_radial <- function(rf, pixel_width, smooth_sd = 0.7, SNR_min = 2) {
  stopifnot(is.matrix(rf), length(pixel_width) == 1, is.finite(pixel_width), pixel_width > 0)
  nr <- nrow(rf); nc <- ncol(rf)
  dr <- pixel_width / 2  # internal radial bin width

  # --- helpers ---
  gaussian_blur_imager_reflect <- function(M, sd_pix) {
    if (sd_pix <= 0) return(M)
    im <- as.cimg(M)
    sm <- isoblur(im, sd_pix)  # isotropic Gaussian blur
    as.array(sm)[,,1,1]
  }

  # Quadratic (3x3) sub-pixel peak on smoothed map
  subpixel_peak <- function(M, r0, c0) {
    if (r0 <= 1 || r0 >= nrow(M) || c0 <= 1 || c0 >= ncol(M))
      return(c(r0, c0, M[r0, c0]))
    z <- M[(r0 - 1):(r0 + 1), (c0 - 1):(c0 + 1)]
    yy <- rep(-1:1, each = 3)
    xx <- rep(-1:1, times = 3)
    A <- cbind(xx ^ 2, yy ^ 2, xx * yy, xx, yy, 1)
    cf <-
      try(coef(lm(as.vector(z) ~ A - 1)), silent = TRUE)
    if (inherits(cf, "try-error"))
      return(c(r0, c0, M[r0, c0]))
    a <-
      cf[1]
    b <- cf[2]
    cxy <- cf[3]
    d <- cf[4]
    e <- cf[5]
    f0 <- cf[6]
    D <- (2 * a) * (2 * b) - cxy ^ 2
    if (!is.finite(D) || abs(D) < 1e-9)
      return(c(r0, c0, M[r0, c0]))
    xh <- (cxy * e - 2 * b * d) / D
    xh <- max(-1, min(1, xh))
    yh <- (cxy * d - 2 * a * e) / D
    yh <- max(-1, min(1, yh))
    amp <- a * xh ^ 2 + b * yh ^ 2 + cxy * xh * yh + d * xh + e * yh + f0
    c(r0 + yh, c0 + xh, amp)
  }
  # r50 from a (r, v) profile; v should be peak-normalized
  r50_from_profile <- function(r, v) {
    ok <- is.finite(r) & is.finite(v); r <- r[ok]; v <- v[ok]
    if (!length(r)) return(NA_real_)
    o <- order(r); r <- r[o]; v <- v[o]
    if (length(v) >= 3) {
      v <- filter(v, rep(1/3, 3), sides = 2)
      v[1] <- v[2]; v[length(v)] <- v[length(v)-1]
      v <- as.numeric(v)
    }
    v <- rev(cummax(rev(v)))  # enforce non-increasing outward
    idx <- which(v <= 0.5)[1]
    if (is.na(idx)) return(NA_real_)
    if (idx == 1)  return(r[1])
    r0 <- r[idx-1]; r1 <- r[idx]; v0 <- v[idx-1]; v1 <- v[idx]
    r0 + (0.5 - v0) * (r1 - r0) / (v1 - v0)
  }

  # --- 1) smooth + subpixel peak ---
  rf_d<-denoise.modwt.2d(rf, wf = "la8", J = 3,
                   method = "universal", rule = "soft")
  rf_s <- if (smooth_sd > 0) rf_s <- gaussian_blur_imager_reflect(rf_d, smooth_sd) else rf_d
df <- data.frame(
  Row = rep(seq_len(nr), each = nr),
  Col = rep(seq_len(nr), times = nr),
  Response = as.vector(rf_s)
)

  ij <- which(rf_s == max(rf_s, na.rm = TRUE), arr.ind = TRUE)[1, ]
  sp <- subpixel_peak(rf_s, ij[1], ij[2])
  r_pk <- sp[1]; c_pk <- sp[2]; amp_pk <- sp[3]

  # guard against degenerate maps
  if (!is.finite(amp_pk) || amp_pk == 0) {
    return(list(
      peak_row = as.numeric(r_pk),
      peak_col = as.numeric(c_pk),
      peak_amp = as.numeric(amp_pk),
      r50 = NA_real_,
      symmetry_ratio = NA_real_,
      peak_at_margin = (round(r_pk) %in% c(1, nr)) || (round(c_pk) %in% c(1, nc))
    ))
  }

  # --- coordinates (centers) in physical units (square pixels) ---
  x <- (seq_len(nc) - 0.5) * pixel_width
  y <- (seq_len(nr) - 0.5) * pixel_width
  X <- matrix(rep(x, each = nr), nrow = nr)
  Y <- matrix(rep(y, times = nc), nrow = nr)
  x_pk <- (c_pk - 0.5) * pixel_width
  y_pk <- (r_pk - 0.5) * pixel_width

  # vectors
  V  <- as.vector(rf_s)
  R  <- sqrt((as.vector(X) - x_pk)^2 + (as.vector(Y) - y_pk)^2)
  TH <- (atan2(as.vector(Y) - y_pk, as.vector(X) - x_pk) * 180 / pi) %% 360

  # normalize to peak
  Vn <- V / amp_pk

  # --- 2) radially-flattened r50 (use available pixels only) ---
  rmax  <- max(R[is.finite(Vn)], na.rm = TRUE)
  edges <- seq(0, rmax + dr, by = dr)
  bin   <- cut(R, edges, right = FALSE, include.lowest = TRUE)
  prof  <- tapply(Vn, bin, function(z) mean(z, na.rm = TRUE))
  r_cent <- head(edges, -1) + dr/2
  r50   <- r50_from_profile(r_cent, as.numeric(prof))

  # --- 3) sector r50s (12 bins of 30°) and symmetry ratio ---
  nsec <- 12L; width <- 360 / nsec  # 30°
  r50_sec <- rep(NA_real_, nsec)
  for (s in 0:(nsec-1)) {
    a0 <- s * width; a1 <- a0 + width
    sel <- (TH >= a0 & TH < a1) | (s == nsec-1 & TH == 360)  # include 360 in last bin
    if (!any(sel, na.rm = TRUE)) next
    Rb <- R[sel]; Vb <- Vn[sel]
    bin_s  <- cut(Rb, edges, right = FALSE, include.lowest = TRUE)
    prof_s <- tapply(Vb, bin_s, function(z) mean(z, na.rm = TRUE))
    r50_sec[s+1] <- r50_from_profile(r_cent, as.numeric(prof_s))
  }
  ix <- function(i, n) ((i - 1L) %% n) + 1L
  i_max <- which.max(r50_sec)

  i_maj2   <- ix(i_max + nsec/2L, nsec)     # opposite of max (θ + 180°)
  i_orth1  <- ix(i_max + nsec/4L, nsec)     # θ + 90°
  i_orth2  <- ix(i_max - nsec/4L, nsec)     # θ - 90° (≡ θ + 270°)

  major <- mean(c(r50_sec[i_max],  r50_sec[i_maj2]),  na.rm = TRUE)  # major diameter
  minor <- mean(c(r50_sec[i_orth1], r50_sec[i_orth2]), na.rm = TRUE) # orthogonal diameter

  symmetry_ratio <- if (is.finite(major) && is.finite(minor) && minor > 0)
    major / minor else NA_real_

  # Outer ring: use outer 30% of available radii as baseline region
  outer_sel <- R >= (max(R[is.finite(V)], na.rm=TRUE) * 0.7)
  base_med  <- median(V[outer_sel], na.rm=TRUE)
  base_mad  <- mad(V[outer_sel], na.rm=TRUE, constant = 1) * 1.4826
  SNR_pk    <- (amp_pk - base_med) / (base_mad + 1e-9)
  if (!is.finite(SNR_pk) || SNR_pk < SNR_min || amp_pk < SNR_min) {
    return(list(
      peak_row = NA_real_, peak_col = NA_real_,
      peak_amp = NA_real_, r50 = NA_real_, symmetry_ratio = NA_real_,
      peak_at_margin = FALSE
    ))
  }

  list(
    peak_row = as.numeric(r_pk),
    peak_col = as.numeric(c_pk),
    peak_amp = as.numeric(amp_pk),
    r50 = r50,                              # distance units of pixel_width
    symmetry_ratio = symmetry_ratio,        # max sector r50 / orthogonal sector r50
    peak_at_margin = (round(r_pk) %in% c(1, nr)) || (round(c_pk) %in% c(1, nc))
  )
}
