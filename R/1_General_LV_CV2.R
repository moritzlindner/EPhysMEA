#' Spike‐train irregularity metrics: Shinomoto LV, CV2, and segment summaries
#'
#' Functions to quantify local irregularity of spike trains using
#' Shinomoto's LV and the CV2 measure, plus a helper to aggregate these
#' metrics within user‐defined time segments.
#'
#' @description
#' For successive inter‐spike intervals (ISIs) \eqn{\Delta_i}, the local
#' Shinomoto LV values are
#' \deqn{ LV_i = 3 \frac{(\Delta_i - \Delta_{i+1})^2}{(\Delta_i + \Delta_{i+1})^2}. }
#' The local CV2 values are
#' \deqn{ CV2_i = 2 \frac{|\Delta_i - \Delta_{i+1}|}{\Delta_i + \Delta_{i+1}} \in [0,2], }
#' and we often summarize \eqn{|CV2 - 1|} (distance from Poisson, in \eqn{[0,1]}).
#'
#' @section Metrics provided:
#' \describe{
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{Vector of local LV values (one per adjacent ISI pair).}
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{Mean of \code{spike_irreg_lv_shinomoto(times)}.}
#'   \item{\code{spike_irreg_cv2_local(times)}}{Vector of \eqn{|CV2-1|} values (one per adjacent ISI pair).}
#'   \item{\code{cv2(times)}}{Typical single‐value summary of \code{spike_irreg_cv2_local(times)}
#'     (see \emph{Value} for details).}
#'   \item{\code{segment_variability_by_breaks(times, breaks, metric = c("lv","cv2"), ...)}}{
#'     Aggregates local LV or \eqn{|CV2-1|} within time segments defined by \code{breaks}
#'     and returns a per‐segment summary table; optionally includes per‐segment
#'     p‐values from a KS test against the Poissonized null
#'     (Uniform\eqn{(0,1)} for \eqn{|CV2-1|}, and after mapping LV via
#'     \eqn{CV2 = \sqrt{4/3\, LV}}).}
#' }
#'
#' @param times Numeric vector of spike times (seconds or samples).
#'   Must be sorted in ascending order for \code{lv_shinomoto_*} and \code{cv2_*}.
#'   Non‐finite times are ignored by the segment helper.
#' @param breaks Numeric vector of cut points defining closed/open intervals
#'   for segmentation (see \code{\link[base]{cut}}). Must have length \eqn{\ge 2}.
#' @param metric Character scalar, one of \code{"lv"} or \code{"cv2"}, selecting
#'   which local measure to aggregate per segment.
#' @param labels Optional labels passed to \code{\link[base]{cut}}.
#' @param right,include.lowest Logicals forwarded to \code{\link[base]{cut}} to
#'   control interval closure.
#' @param p_values Logical; if \code{TRUE}, compute per‐segment KS p‐values
#'   against the Uniform\eqn{(0,1)} null (see Details above).
#' @param ... Additional arguments passed to \code{\link[base]{cut}}.
#'
#' @details
#' LV is relatively insensitive to slow firing‐rate modulations and captures
#' local ISI irregularity. \eqn{|CV2-1|} summarizes deviation from Poisson‐like
#' irregularity on a \eqn{[0,1]} scale.
#'
#' @return
#' \describe{
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{Numeric vector of local LV values;
#'     length \code{length(times) - 2}. If fewer than three spikes, returns \code{numeric(0)}.
#'     Non‐finite values arising from 0/0 are returned as \code{NA}.}
#'   \item{\code{spike_irreg_lv_shinomoto(times)}}{Single numeric value: the mean of
#'     \code{spike_irreg_lv_shinomoto(times)}. If fewer than three spikes, returns
#'     \code{numeric(0)}. If the mean is \code{NA}, returns \code{1}.}
#'   \item{\code{spike_irreg_cv2_local(times)}}{Numeric vector of \eqn{|CV2-1|} values;
#'     length \code{length(times) - 2}. If fewer than three spikes, returns
#'     \code{numeric(0)}. Non‐finite values are returned as \code{NA}.}
#'   \item{\code{cv2(times)}}{Single numeric summary of \code{spike_irreg_cv2_local(times)}
#'     (typically the median). If fewer than three spikes, returns \code{numeric(0)}.
#'     If the summary is \code{NA}, returns \code{0}.}
#'   \item{\code{segment_variability_by_breaks(...)}}{A \code{data.frame} with one
#'     row per segment and columns:
#'     \code{seg_id}, \code{seg_start}, \code{seg_end}, \code{seg_label},
#'     \code{n_pairs} (non‐NA local pairs used), \code{value} (segment aggregate),
#'     and optionally \code{p_value} when \code{p_values = TRUE}.}
#' }
#'
#' @references
#' Shinomoto, S., et al. (2009). Relating neuronal firing patterns to functional
#' differentiation of cerebral cortex. \emph{PLoS Comput Biol}, 5(7):e1000433.
#' \doi{10.1371/journal.pcbi.1000433}
#'
#' @aliases spike_irreg_lv_shinomoto spike_irreg_lv_shinomoto cv2 spike_irreg_cv2_local segment_variability_by_breaks
#' @seealso \code{\link[base]{cut}}
#' @examples
#' times <- c(0, 0.10, 0.21, 0.33, 0.60)
#' spike_irreg_lv_shinomoto(times)
#' spike_irreg_lv_shinomoto(times)
#' spike_irreg_cv2_local(times)
#' cv2(times)
#' segment_variability_by_breaks(times, breaks = c(0, 0.3, 0.7), metric = "lv")
#'
#' @name LV_CV2
NULL

#' @rdname LV_CV2
#' @export
#' @return For \code{spike_irreg_lv_shinomoto()}: a numeric vector of local LV values,
#'   one per adjacent ISI pair (\code{length(times) - 2}). If fewer than 3 spikes,
#'   returns \code{numeric(0)}. Non-finite values arising from 0/0 are returned as \code{NA}.
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

#' @rdname LV_CV2
#' @export
#' @return For \code{spike_irreg_lv_shinomoto()}: a single numeric value, the mean LV over
#'   \code{spike_irreg_lv_shinomoto(times)}. If fewer than 3 spikes, returns
#'   \code{numeric(0)}. If the mean is \code{NA}, returns \code{1}.
#' @export
spike_irreg_lv_shinomoto <- function(times) {
  times <- as.numeric(times)
  if (!length(times)) return(numeric(0))
  if (is.unsorted(times)) times <- sort(times)

  lv_vec <- spike_irreg_lv_shinomoto_local(times)
  if (!length(lv_vec)) return(numeric(0))
  m <- mean(lv_vec, na.rm = TRUE)
  if (is.na(m)) 1 else m
}
#' @rdname LV_CV2
#' @return For \code{spike_irreg_cv2_local()}: a numeric vector of \eqn{|CV2-1|} values,
#'   one per adjacent ISI pair. If fewer than 3 spikes, returns \code{numeric(0)}.
#'   Non-finite values are returned as \code{NA}.
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

#' @rdname LV_CV2
#' @return For \code{cv2()}: a single numeric summary of \code{spike_irreg_cv2_local(times)}
#'   (typically the median). If fewer than 3 spikes, returns \code{numeric(0)}.
#'   If the summary is \code{NA}, returns \code{0}.
#' @importFrom stats median
#' @export
cv2 <- function(times) {
  vals <- cv2_times(times)
  if (!length(vals)) return(numeric(0))
  m <- median(vals, na.rm = TRUE) #  m <- mean(vals, na.rm = TRUE)
  if (is.na(m)) 0 else m
}
