#' Step-response metrics (Farrow-style) + ISI-drop latency + Shinomoto
#'
#' Computes Farrow ON/OFF counts and indices, a transience index from a
#' peak-normalized PSTH, an ISI-drop latency relative to onset/offset (or auto),
#' a simple response amplitude (mean rate in the transience window), and
#' baseline statistics (rate and Shinomoto LV).
#'
#' @param x            Numeric (or units) vector of spike/event times (seconds).
#' @param step_onset   Numeric (or units) step onset time (seconds).
#' @param step_offset  Numeric (or units) step offset time (seconds).
#' @param on_win       Length of ON counting window (s), default 0.4.
#' @param off_win      Length of OFF counting window (s), default 0.4.
#' @param trans_win    Length of transience PSTH window (s), default 1.0.
#' @param binwidth     PSTH bin width (s), default 0.02 (20 ms).
#' @param smooth       Logical, Gaussian-smooth PSTH before AUC? default FALSE.
#' @param gauss_sd     Gaussian kernel SD (s) if smoothing, default 0.02.
#' @param return_psth  If TRUE, returns a small data.frame with the PSTH used.
#' @param latency_alpha       Quantile for ISI threshold (e.g., 0.05). Skip latency calculation if NA.
#' @param latency_preduration Pre-event baseline duration (s) used for robust ISI fallback.
#' @param latency_postduration Post-event horizon (s) for robust fallback; default = trans_win.
#' @param latency_event  Which event to use for latency: "auto" (default: onset if Bias>0 else offset),
#'   "onset", or "offset".
#' @param shinomoto_win Window for baseline firing rate and baseline stability (Shinomoto).
#'   Skip Shinomoto calculation if NA.
#'
#' @importFrom units set_units
#' @importFrom stats quantile filter dnorm
#'
#' @return
#' A named list with:
#' \describe{
#'   \item{Farrow_ON}{Spike count in the ON window \code{[step_onset, step_onset + on_win)}.}
#'   \item{Farrow_OFF}{Spike count in the OFF window \code{[step_offset, step_offset + off_win)}.}
#'   \item{Farrow_BiasIndex}{\code{(Farrow_ON - Farrow_OFF) / (Farrow_ON + Farrow_OFF)}; 0 if both are 0.}
#'   \item{Farrow_TransienceIndex}{Area under the peak-normalized PSTH over \code{trans_win},
#'         divided by \code{trans_win} (dimensionless).}
#'   \item{ISIdrop_Latency}{Latency (s) from the chosen event to the time of the second spike
#'         in the first post-event ISI that is “short”; short if
#'         \code{postISI < max(Q_alpha, 0.6 * median(ISI_pre))}.}
#'   \item{Response_Amplitude}{Mean spike rate (spikes/s) in the transience window
#'         \code{[trans_start, trans_start + trans_win)}, where \code{trans_start} is
#'         \code{step_onset} if \code{Farrow_BiasIndex > 0} and \code{step_offset} otherwise.}
#'   \item{Baseline_Rate}{Mean baseline rate in \code{[step_onset - on_win, step_onset)}.}
#'   \item{Shinomoto_Baseline_LV}{Local Variation from adjacent ISIs in the baseline window.}
#'   \item{psth}{(optional) A data.frame with columns \code{t_left}, \code{t_right}, \code{rate}
#'         for the PSTH used in the transience calculation.}
#' }
#'
#' @section Helpers (documented on this page):
#' \describe{
#'   \item{\code{step_response_amplitude(x, trans_win)}}{Compute mean rate in a window.}
#'   \item{\code{step_response_latency(x, event, latency_alpha, latency_preduration, latency_postduration)}}{ISI-drop latency.}
#' }
#'
#' (For Farrow metrics, see https://pmc.ncbi.nlm.nih.gov/articles/PMC3075295/#sec1.
#' For LV see 10.1016/j.biosystems.2004.09.023.)
#'
#' @name step_stim_resp_metrics
#' @export
step_stim_resp_metrics <- function(x,
                                   step_onset,
                                   step_offset,
                                   on_win = 0.4,
                                   off_win = 0.4,
                                   trans_win = 1.0,
                                   binwidth = 0.02,
                                   smooth = FALSE,
                                   gauss_sd = 0.02,
                                   return_psth = FALSE,
                                   latency_alpha = 0.05,
                                   latency_preduration = 0.5,
                                   latency_postduration = NULL,
                                   latency_event = c("auto", "onset", "offset"),
                                   shinomoto_win = 2) {
  latency_event <- match.arg(latency_event)
  if (is.null(latency_postduration))
    latency_postduration <- trans_win

  x                   <- to_num(x)
  step_onset          <- to_num(step_onset)
  step_offset         <- to_num(step_offset)
  on_win              <- to_num(on_win)
  off_win             <- to_num(off_win)
  trans_win           <- to_num(trans_win)
  binwidth            <- to_num(binwidth)
  gauss_sd            <- to_num(gauss_sd)
  latency_alpha       <- to_num(latency_alpha)
  latency_preduration <- to_num(latency_preduration)
  latency_postduration <- to_num(latency_postduration)
  shinomoto_win        <- to_num(shinomoto_win)

  stopifnot(is.numeric(x),
            is.numeric(step_onset),
            is.numeric(step_offset))
  stopifnot(
    is.finite(on_win)  && on_win  > 0,
    is.finite(off_win) && off_win > 0,
    is.finite(trans_win) && trans_win > 0,
    is.finite(latency_preduration) &&
      latency_preduration > 0,
    is.finite(latency_postduration) &&
      latency_postduration > 0
  )

  x <- x[x < step_offset + max(off_win, trans_win, latency_postduration)]
  x <- x[x > step_onset - max(on_win, latency_preduration)]
  x <- sort(x[is.finite(x)])

  # --- Farrow ON / OFF and Bias ---
  Farrow_ON  <- count_in(x, step_onset,  step_onset  + on_win)
  Farrow_OFF <- count_in(x, step_offset, step_offset + off_win)
  den <- Farrow_ON + Farrow_OFF
  Farrow_BiasIndex <-
    if (den > 0)      {
      (Farrow_ON - Farrow_OFF) / den
    }   else  {
      0
    }

  # --- Transience AUC over chosen 1 s window (same rule as before) ---
  trans_start <-
    if (Farrow_BiasIndex > 0)      {
      step_onset
    }  else    {
      step_offset
    }
  trans_end   <- trans_start + trans_win
  in_win <- x[x >= trans_start & x < trans_end]

  if (!is.na(binwidth)){
    stopifnot(is.finite(binwidth) && binwidth > 0)
    n_breaks <- floor((trans_end - trans_start) / binwidth)
    breaks   <- trans_start + seq(0, n_breaks) * binwidth
    if (tail(breaks, 1) < trans_end - 1e-12)   { breaks <- c(breaks, trans_end)}


    h <-
      hist(
        in_win,
        breaks = breaks,
        right = FALSE,
        plot = FALSE,
        include.lowest = TRUE
      )
    rate <- h$counts / diff(h$breaks) # spikes/s

    if (smooth && length(rate)) {
      k_idx  <- seq(-3 * gauss_sd, 3 * gauss_sd, by = binwidth)
      kernel <- dnorm(k_idx, sd = gauss_sd)
      kernel <- kernel / sum(kernel)
      r_sm <- filter(rate, kernel, sides = 2, circular = FALSE)
      if (anyNA(r_sm)) {
        nz <- which(!is.na(r_sm))
        if (length(nz)) {
          left_gap  <- seq_len(min(nz) - 1L)
          right_gap <- if (max(nz) < length(r_sm)) (max(nz) + 1L):length(r_sm) else integer(0)
          if (length(left_gap))  r_sm[left_gap]  <- r_sm[min(nz)]
          if (length(right_gap)) r_sm[right_gap] <- r_sm[max(nz)]
        } else {
          r_sm <- rate
        }
      }
      rate <- as.numeric(r_sm)
    }

    mx <- suppressWarnings(max(rate, na.rm = TRUE))
    if (is.finite(mx) && mx > 0) {
      rate <- rate / mx                 # peak-normalize to 1
    } else {
      rate[] <- 0
    }

    Farrow_TransienceIndex <-
      if (length(rate))    {
        sum(rate) * binwidth
      }  else  {
        0
      }

    Farrow_TransienceIndex <- (sum(rate) * binwidth) / trans_win
  } else {
    Farrow_ON <- NA
    Farrow_OFF <- NA
    Farrow_TransienceIndex <- NA
  }

  if (!is.na(latency_alpha)) {
    evt_for_latency <-
      if (latency_event == "onset") {
        step_onset
      } else if (latency_event == "offset") {
        step_offset
      } else {
        if (Farrow_BiasIndex > 0)
          step_onset
        else
          step_offset
      }
    if (!(latency_alpha > 0 && latency_alpha < 1)) {
      stop ("'latency_alpha' must be a numeric in the range of 0 to 1 or NA.")
    }
    ISIdrop_Latency <- step_response_latency(x,
                                             evt_for_latency,
                                             latency_alpha,
                                             latency_preduration,
                                             latency_postduration)
  } else {
    ISIdrop_Latency <- NA
  }
  # Response Amplitude
  Response_Amplitude <- step_response_amplitude(in_win, trans_win)

  # Resting stability
  base_times <- x[x >= (step_onset - shinomoto_win) & x < step_onset]
  if (!is.na(shinomoto_win)) {
    Baseline_Rate <- length(base_times) / shinomoto_win
    if(length(base_times) == 0){
      Shinomoto_Baseline_LV <- 1
    } else {
      Shinomoto_Baseline_LV <- spike_irreg_lv_shinomoto(base_times)
    }
    if (length(Shinomoto_Baseline_LV) == 0L)
      Shinomoto_Baseline_LV <- 1
  } else {
    Baseline_Rate <- NA
    Shinomoto_Baseline_LV <- NA
  }

  out <- list(
    Farrow_ON              = as.numeric(Farrow_ON),
    Farrow_OFF             = as.numeric(Farrow_OFF),
    Farrow_BiasIndex       = as.numeric(Farrow_BiasIndex),
    Farrow_TransienceIndex = as.numeric(Farrow_TransienceIndex),
    ISIdrop_Latency         = as.numeric(ISIdrop_Latency),
    Response_Amplitude      = as.numeric(Response_Amplitude),
    Baseline_Rate           = as.numeric(Baseline_Rate),
    Shinomoto_Baseline_LV   = as.numeric(Shinomoto_Baseline_LV)
  )

  if (isTRUE(return_psth)) {
    out$psth <- data.frame(
      t_left  = h$breaks[-length(h$breaks)],
      t_right = h$breaks[-1],
      rate    = rate,
      check.names = FALSE
    )
  }

  out
}


#' Helper: response amplitude (mean rate) in a window
#'
#' Computes the average spike rate (spikes/s) given the spikes that fall inside a window
#' and the window length.
#'
#' @param x Numeric vector of spike times already restricted to the analysis window (can be length 0).
#' @param trans_win Window length (seconds). Numeric or `units`.
#'
#' @return A numeric scalar: spikes per second; `NA_real_` if `trans_win` is invalid.
#'
#' @rdname step_stim_resp_metrics
#' @export
step_response_amplitude <- function(x, trans_win) {
  trans_win <- to_num(trans_win)

  if (!is.finite(trans_win) || trans_win <= 0)
    return(NA_real_)

  ra <- length(x) / trans_win
  if (!is.finite(ra))
    NA_real_
  else
    as.numeric(ra)
}


#' Helper: ISI-drop latency detector
#'
#' Latency (s) from `event` to the second spike of the first post-event ISI that is
#' shorter than a threshold derived from pre-event ISIs:
#' \deqn{postISI < \max(Q_\alpha, 0.6 \cdot \mathrm{median(ISI_{pre})}).}
#'
#' The post-event search is restricted to \code{[event, event + latency_postduration]}.
#'
#' @param x Numeric vector of spike/event times (s), single train.
#' @param event Numeric event time (s) to compute latency from.
#' @param latency_alpha Quantile for the low-ISI threshold (e.g., 0.05).
#' @param latency_preduration Duration (s) of the pre-event baseline window.
#' @param latency_postduration Post-event horizon (s) to search for short ISIs.
#'
#' @return Latency in seconds or `NA_real_` if no short-ISI response is detected.
#' @importFrom stats quantile median
#' @rdname step_stim_resp_metrics
#' @export
step_response_latency <- function(x,
                                  event,
                                  latency_alpha,
                                  latency_preduration,
                                  latency_postduration) {
  # Pre-event baseline ISIs
  pre <- x[x >= (event - latency_preduration) & x < event]
  if (length(pre) >= 2L) {
    preISIs <- diff(pre)
  } else if (length(pre) == 1L) {
    preISIs <- latency_preduration / 2
  } else {
    preISIs <- latency_preduration
  }
  if (!length(preISIs) || !all(is.finite(preISIs)))
    return(NA_real_)

  # Threshold from baseline
  low_q <-
    as.numeric(quantile(preISIs, probs = latency_alpha, names = FALSE))
  med   <- median(preISIs)
  ISIthr <- max(low_q, 0.60 * med)

  # Post-event spikes, restricted to the horizon
  rel <- x[x >= event & x < (event + latency_postduration)] - event

  # Need at least two spikes within the window to form a post-event ISI
  if (length(rel) < 2L)
    return(NA_real_)

  postISIs <- diff(rel)
  idx <- which(postISIs < ISIthr)[1]

  if (is.finite(idx))
    rel[idx + 1]
  else
    NA_real_
}



#' left Inclusive, right exclusive window counting
#' @keywords internal
#' @noMd
#' @param a,b Left () and right margins of window
count_in <- function(x, a, b) {
  if (!length(x))
    return(0L)
  sum(x >= a & x < b)
}

