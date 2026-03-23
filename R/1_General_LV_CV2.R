#' Spike-train irregularity metrics: Shinomoto LV and CV2
#'
#' Functions to quantify local irregularity of spike trains using
#' Shinomoto's LV and the CV2 measure.
#'
#' @description
#' For successive inter-spike intervals (ISIs) \eqn{\Delta_i}, the local
#' Shinomoto LV values are
#' \deqn{ LV_i = 3 \frac{(\Delta_i - \Delta_{i+1})^2}{(\Delta_i + \Delta_{i+1})^2}. }
#' The local CV2 values are
#' \deqn{ CV2_i = 2 \frac{|\Delta_i - \Delta_{i+1}|}{\Delta_i + \Delta_{i+1}} \in [0,2], }
#' and this implementation summarizes \eqn{|CV2 - 1|}, i.e. the distance from
#' Poisson-like irregularity on a scale from 0 to 1.
#'
#' @section Metrics provided:
#' \describe{
#'   \item{\code{spike_irreg_lv_shinomoto_local(times)}}{
#'     Vector of local LV values, one per adjacent ISI pair.}
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{
#'     Mean of \code{spike_irreg_lv_shinomoto_local(times)}.}
#'   \item{\code{spike_irreg_cv2_local(times)}}{
#'     Vector of local \eqn{|CV2 - 1|} values, one per adjacent ISI pair.}
#'   \item{\code{spike_irreg_cv2(times)}}{
#'     Single-value summary of \code{spike_irreg_cv2_local(times)}
#'     (currently the median).}
#' }
#'
#' @param times Numeric vector of spike times (seconds or samples).
#'   Spike times should be sorted in ascending order.
#'   \code{spike_irreg_lv_shinomoto()} sorts unsorted input internally;
#'   the other functions assume sorted input.
#'
#' @details
#' LV is relatively insensitive to slow firing-rate modulations and captures
#' local ISI irregularity. \eqn{|CV2 - 1|} summarizes deviation from
#' Poisson-like irregularity on a \eqn{[0,1]} scale.
#'
#' @return
#' \describe{
#'   \item{\code{spike_irreg_lv_shinomoto_local(times)}}{
#'     Numeric vector of local LV values; length \code{length(times) - 2}.
#'     If fewer than three spikes are available, returns \code{numeric(0)}.
#'     Non-finite values arising from undefined ratios such as 0/0 are returned
#'     as \code{NA_real_}.}
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{
#'     Single numeric value: the mean of
#'     \code{spike_irreg_lv_shinomoto_local(times)}.
#'     If fewer than three spikes are available, returns \code{NA_real_}.
#'     If all local values are non-finite and are removed during aggregation,
#'     returns \code{NA_real_}.}
#'   \item{\code{spike_irreg_cv2_local(times)}}{
#'     Numeric vector of local \eqn{|CV2 - 1|} values;
#'     length \code{length(times) - 2}.
#'     If fewer than three spikes are available, returns \code{numeric(0)}.
#'     Non-finite values are returned as \code{NA_real_}.}
#'   \item{\code{spike_irreg_cv2(times)}}{
#'     Single numeric summary of \code{spike_irreg_cv2_local(times)}
#'     (currently the median).
#'     If fewer than three spikes are available, returns \code{NA_real_}.
#'     If all local values are non-finite and are removed during aggregation,
#'     returns \code{NA_real_}.}
#' }
#'
#' @references
#' Shinomoto, S., et al. (2009). Relating neuronal firing patterns to functional
#' differentiation of cerebral cortex. \emph{PLoS Comput Biol}, 5(7):e1000433.
#' \doi{10.1371/journal.pcbi.1000433}
#'
#' @aliases spike_irreg_lv_shinomoto_local spike_irreg_lv_shinomoto spike_irreg_cv2_local spike_irreg_cv2
#'
#' @examples
#' times <- c(0, 0.10, 0.21, 0.33, 0.60)
#' spike_irreg_lv_shinomoto_local(times)
#' spike_irreg_lv_shinomoto(times)
#' spike_irreg_cv2_local(times)
#' spike_irreg_cv2(times)
#'
#' @name Spike_Irregularity
NULL

#' @rdname Spike_Irregularity
#' @return For \code{spike_irreg_lv_shinomoto_local()}: a numeric vector of local
#'   LV values, one per adjacent ISI pair (\code{length(times) - 2}).
#'   If fewer than 3 spikes are available, returns \code{numeric(0)}.
#'   Non-finite values arising from undefined ratios such as 0/0 are returned
#'   as \code{NA_real_}.
#' @export
spike_irreg_lv_shinomoto_local <- function(times) {
  isi <- diff(times)
  if (length(isi) < 2L) return(numeric(0))
  r <- isi[-length(isi)]
  s <- isi[-1L]
  out <- 3 * ((r - s)^2) / ((r + s)^2)
  out[!is.finite(out)] <- NA_real_
  out
}

#' @rdname Spike_Irregularity
#' @return For \code{spike_irreg_lv_shinomoto()}: a single numeric value, the mean
#'   LV over \code{spike_irreg_lv_shinomoto_local(times)}.
#'   If fewer than 3 spikes are available, returns \code{NA_real_}.
#'   If all local values are non-finite and are removed during aggregation,
#'   returns \code{NA_real_}.
#' @export
spike_irreg_lv_shinomoto <- function(times) {
  times <- as.numeric(times)
  if (!length(times)) return(NA_real_)
  if (is.unsorted(times)) times <- sort(times)

  lv_vec <- spike_irreg_lv_shinomoto_local(times)
  if (!length(lv_vec)) return(NA_real_)

  m <- mean(lv_vec, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

#' @rdname Spike_Irregularity
#' @return For \code{spike_irreg_cv2_local()}: a numeric vector of
#'   \eqn{|CV2 - 1|} values, one per adjacent ISI pair (\code{length(times) - 2}).
#'   If fewer than 3 spikes are available, returns \code{numeric(0)}.
#'   Non-finite values are returned as \code{NA_real_}.
#' @export
spike_irreg_cv2_local <- function(times) {
  isi <- diff(times)
  if (length(isi) < 2L) return(numeric(0))
  r <- isi[-length(isi)]
  s <- isi[-1L]
  out <- 2 * abs(r - s) / (r + s)      # CV2 in [0,2]
  out <- abs(out - 1)                  # distance from Poisson, in [0,1]
  out[!is.finite(out)] <- NA_real_
  out
}

#' @rdname Spike_Irregularity
#' @return For \code{spike_irreg_cv2()}: a single numeric summary of
#'   \code{spike_irreg_cv2_local(times)} (currently the median).
#'   If fewer than 3 spikes are available, returns \code{NA_real_}.
#'   If all local values are non-finite and are removed during aggregation,
#'   returns \code{NA_real_}.
#' @importFrom stats median
#' @export
spike_irreg_cv2 <- function(times) {
  vals <- spike_irreg_cv2_local(times)
  if (!length(vals)) return(NA_real_)

  m <- median(vals, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}
