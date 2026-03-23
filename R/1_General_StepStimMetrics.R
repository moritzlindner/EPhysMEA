#' Step-stimulus response metrics for `EPhysEvents`
#'
#' Compute step-response metrics for each \strong{run × channel} in an
#' \code{\link[EPhysData:EPhysEvents-class]{EPhysEvents}} object and return them
#' as a tidy wide \code{data.frame}. The function applies
#' \code{\link{step_stim_resp_metrics}} to every run/channel combination and can
#' optionally average the resulting metrics across repeated runs that share the
#' same \code{RecordingID}.
#'
#' This page documents both the high-level \code{StepStimMetrics()} method and
#' the underlying metric definitions used internally.
#'
#' @param X An \code{EPhysEvents} instance.
#' @param step_range Numeric or \code{units} length-2 vector \code{c(start, end)}
#'   defining the step interval. It is interpreted as \code{[start, end)} and
#'   must satisfy \code{start < end}. If given as \code{units}, it is converted
#'   to \code{\link[EPhysData:TimeUnits]{TimeUnits}(X)}.
#' @param average_by_recordingID Logical; if \code{TRUE}, average all numeric
#'   metric columns across rows sharing the same \code{RecordingID} and
#'   \code{Channel} using \code{mean(na.rm = TRUE)}. Default: \code{FALSE}.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @inheritDotParams step_stim_resp_metrics
#'   on_win off_win trans_win binwidth smooth gauss_sd return_psth
#'   latency_alpha latency_preduration latency_postduration latency_event
#'   shinomoto_win
#' @param ... Further arguments forwarded to \code{\link{step_stim_resp_metrics}}.
#'
#' @details
#' \strong{Iteration level:}
#' \code{StepStimMetrics()} iterates over the contents of \code{X} via
#' \code{EPhysData::lapply(X, ...)} and computes the metrics separately for each
#' run and channel. The resulting list output is flattened into a wide
#' \code{data.frame}.
#'
#' \strong{Choice of ON-like versus OFF-like response window:}
#' several metrics depend on whether a unit behaves more ON-like or OFF-like.
#' This is decided from the Farrow bias index:
#' \deqn{
#' \mathrm{Farrow\_BiasIndex} =
#' \frac{\mathrm{Farrow\_ON} - \mathrm{Farrow\_OFF}}
#'      {\mathrm{Farrow\_ON} + \mathrm{Farrow\_OFF}}
#' }
#' with the value defined as 0 when both counts are 0.
#'
#' If \code{Farrow_BiasIndex > 0}, the response is treated as ON-dominated and
#' the transience and response-amplitude windows start at \code{step_onset}.
#' Otherwise they start at \code{step_offset}.
#'
#' \strong{Computed metrics:}
#' \describe{
#'   \item{\code{Farrow_ON}}{Spike count in the ON window
#'   \code{[step_onset, step_onset + on_win)}.}
#'
#'   \item{\code{Farrow_OFF}}{Spike count in the OFF window
#'   \code{[step_offset, step_offset + off_win)}.}
#'
#'   \item{\code{Farrow_BiasIndex}}{ON-versus-OFF preference, computed as
#'   \code{(Farrow_ON - Farrow_OFF) / (Farrow_ON + Farrow_OFF)} and set to
#'   \code{0} if both counts are zero. Positive values indicate relatively
#'   stronger ON responses, negative values relatively stronger OFF responses.}
#'
#'   \item{\code{Farrow_TransienceIndex}}{A transience measure derived from a
#'   PSTH in the selected response window of length \code{trans_win}. Spikes are
#'   histogrammed using bins of width \code{binwidth}; optional Gaussian
#'   smoothing can be applied. The PSTH is then peak-normalized so its maximum
#'   equals 1, and the index is computed as the area under that normalized PSTH,
#'   divided by \code{trans_win}. Smaller values indicate a more transient
#'   response, larger values a more sustained response over the analysis window.}
#'
#'   \item{\code{ISIdrop_Latency}}{Latency from a chosen event
#'   (\code{step_onset}, \code{step_offset}, or automatically selected based on
#'   the bias index) to the \emph{second spike} of the first post-event interspike
#'   interval that is shorter than a threshold derived from pre-event baseline
#'   ISIs. The threshold is
#'   \code{max(Q_alpha, 0.6 * median(ISI_pre))}, where \code{Q_alpha} is the
#'   \code{latency_alpha} quantile of pre-event ISIs. Post-event spikes are only
#'   searched within \code{[event, event + latency_postduration)}.}
#'
#'   \item{\code{Response_Amplitude}}{Mean spike rate (spikes/s) over baseline in the selected
#'   response window \code{[trans_start, trans_start + trans_win)}, where
#'   \code{trans_start} is \code{step_onset} for ON-dominated responses and
#'   \code{step_offset} otherwise.}
#'
#'   \item{\code{Baseline_Rate}}{Mean firing rate in the baseline window
#'   immediately preceding step onset:
#'   \code{[step_onset - shinomoto_win, step_onset)}.}
#'
#'   \item{\code{Shinomoto_Baseline_LV}}{Mean local variation (LV) of adjacent
#'   ISIs in the baseline window, computed with
#'   \code{\link{spike_irreg_lv_shinomoto}}.}
#' }
#'
#' If \code{return_psth = TRUE}, the low-level function
#' \code{\link{step_stim_resp_metrics}} additionally returns the PSTH used for
#' the transience calculation as a small \code{data.frame}; when called through
#' \code{StepStimMetrics()}, this object is carried along inside the intermediate
#' list output but is typically not useful for averaging across runs.
#'
#' @return
#' If \code{average_by_recordingID = FALSE} (default), a wide \code{data.frame}
#' with one row per run × channel, containing identifiers such as
#' \code{RunUID}, \code{RecordingID}, and \code{Channel}, plus one column for
#' each metric returned by \code{step_stim_resp_metrics()}.
#'
#' If \code{average_by_recordingID = TRUE}, a wide \code{data.frame} with one row
#' per \code{RecordingID × Channel}; all numeric metric columns are averaged with
#' \code{mean(na.rm = TRUE)}.
#'
#' @section Low-level functions on this page:
#' \describe{
#'   \item{\code{step_stim_resp_metrics(x, ...)}}{Compute the metrics for a single
#'   spike train.}
#'   \item{\code{step_spikerate(x, t)}}{Helper that computes mean spike rate
#'   (spikes/s) within a window from spike count and window duration.}
#'   \item{\code{step_response_latency(x, event, ...)}}{Helper implementing the
#'   ISI-drop latency rule.}
#' }
#'
#' @references
#' Farrow K, Masland RH
#' (2011) Physiological clustering of visual channels in the mouse retina. \emph{ J Neurophysiol}. 2011;105(4):1516-30.
#'
#' Shinomoto S, Shima K, Tanji J. Differences in spiking patterns among cortical
#' neurons. \emph{Biosystems}. 2005;79(1-3):57-65.
#'
#' @examples
#' # Per-(run, channel) metrics:
#' # df <- StepStimMetrics(ev, step_range = c(1, 2))
#'
#' # Averaged by RecordingID:
#' # avg <- StepStimMetrics(
#' #   ev,
#' #   step_range = c(1, 2),
#' #   average_by_recordingID = TRUE,
#' #   latency_alpha = 0.05
#' # )
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels TimeUnits lapply
#' @importFrom units set_units
#' @seealso \code{\link{step_stim_resp_metrics}},
#'   \code{\link{step_spikerate}},
#'   \code{\link{step_response_latency}},
#'   \code{EPhysData::\link[EPhysData:Metadata]{Metadata}},
#'   \code{EPhysData::\link[EPhysData:Channels]{Channels}},
#'   \code{EPhysData::\link[EPhysData:TimeUnits]{TimeUnits}}
#' @export
setGeneric("StepStimMetrics",
           function(X, step_range, average_by_recordingID = FALSE, parallel = !interactive(), ...)
             standardGeneric("StepStimMetrics")
)

#' @rdname StepStimMetrics
#' @export
setMethod("StepStimMetrics", "EPhysEvents",
          function(X, step_range, average_by_recordingID = FALSE, parallel = !interactive(), ...) {

            md  <- Metadata(X)
            dat <- X@Data
            chn <- Channels(X)

            if (length(unique(md$Experiment))!=1) {
              stop("Only EPhysEvents containing one single type of Experiment as defined in the respective metadata column can safely be processed by this method.")
            }

            # Validate step_range (allow units or numeric)
            if ("units" %in% class(step_range)) {
              tryCatch({
                step_range<-set_units(step_range, TimeUnits(X), mode = "standard")
              }, error = function(e){
                stop("Time unit conversion for step_range failed with message: ",e, "Is provided unit compatible with 'TimeUnits(X)'?")
              })
            }
            step_rangen <- suppressWarnings(as.numeric(step_range))

            if (!is.numeric(step_rangen) || length(step_rangen) != 2L || any(!is.finite(step_rangen))) {
              stop("step_range must be a finite numeric (or units) length-2 vector.")
            }
            step_rangen <- sort(step_rangen)
            if (!(step_rangen[2] > step_rangen[1])) stop("step_range must satisfy start < end.")
            step_onset  <- step_rangen[1]
            step_offset <- step_rangen[2]


            # lapply over channels in this trace

            metrics_list <-
              lapply(
                X,
                step_stim_resp_metrics,
                step_onset  = step_onset,
                step_offset = step_offset,
                ...,
                parallel = parallel
              )

            metrics_df <-nested2df(X,metrics_list)

            if (isFALSE(average_by_recordingID)) {
              return(metrics_df)
            } else {
              num_cols <- setdiff(names(metrics_df)[sapply(metrics_df, is.numeric)],
                                  c("RecordingID", "Channel"))
              keys <- c("RecordingID", "Channel")

              out <- aggregate(metrics_df[num_cols],
                               by = metrics_df[keys],
                               FUN = function(x) mean(x, na.rm = TRUE))
              out$RunUID<-"aggregated"
              return(out)

            }
          }
)



#' @rdname StepStimMetrics
#'
#' @param x Numeric vector of spike or event times for a single spike train.
#' @param step_onset Numeric or \code{units} scalar giving step onset time.
#' @param step_offset Numeric or \code{units} scalar giving step offset time.
#' @param on_win Length of the ON counting window. Default: \code{0.4}.
#' @param off_win Length of the OFF counting window. Default: \code{0.4}.
#' @param trans_win Length of the response window used for transience and
#'   response amplitude. Default: \code{1.0}.
#' @param binwidth PSTH bin width used for the transience calculation.
#'   Default: \code{0.02}.
#' @param smooth Logical; if \code{TRUE}, smooth the PSTH with a Gaussian kernel
#'   before peak-normalization and AUC calculation.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'   Default: \code{0.02}.
#' @param return_psth Logical; if \code{TRUE}, include the PSTH used for the
#'   transience calculation in the returned list.
#' @param latency_alpha Quantile used to define the low-ISI threshold for the
#'   latency detector. Set to \code{NA} to skip latency calculation.
#' @param latency_preduration Duration of the pre-event window used to obtain
#'   baseline ISIs for the latency detector.
#' @param latency_postduration Duration of the post-event search window used by
#'   the latency detector. If \code{NULL}, defaults to \code{trans_win}.
#' @param latency_event Character string; one of \code{"auto"},
#'   \code{"onset"}, or \code{"offset"}. In \code{"auto"} mode, onset is used
#'   when \code{Farrow_BiasIndex > 0}, otherwise offset.
#' @param shinomoto_win Length of the baseline window before step onset used for
#'   \code{Baseline_Rate} and \code{Shinomoto_Baseline_LV}. Set to \code{NA} to
#'   skip both quantities.
#'
#' @return For \code{step_stim_resp_metrics()}, a named list containing the
#'   metrics described under \emph{Computed metrics} on this page, and optionally
#'   \code{psth} if \code{return_psth = TRUE}.
#'
#' @importFrom stats quantile filter dnorm median
#' @importFrom units set_units
#' @seealso \code{\link{spike_irreg_lv_shinomoto}}
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
  # Resting stability
  base_times <- x[x >= (step_onset - shinomoto_win) & x < step_onset]
  if (!is.na(shinomoto_win)) {
    Baseline_Rate <- step_spikerate(base_times, base_times)
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

  # Response Amplitude
  Response_Amplitude <- step_spikerate(in_win, trans_win) - Baseline_Rate

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


#' @rdname StepStimMetrics
#'
#' @param x Numeric vector of spike times already restricted to the relevant
#'   analysis window.
#' @param t Duration of the analysis window in seconds, used to convert spike
#'   count to mean spike rate.
#'
#' @return For \code{step_spikerate()}, a numeric scalar giving mean spike rate
#'   in spikes/s, or \code{NA_real_} if the window duration is invalid.
#' @export
step_spikerate <- function(x, t) {
  t <- to_num(t)

  if (!is.finite(t) || t <= 0)
    return(NA_real_)

  ra <- length(x) / t
  if (!is.finite(ra))
    NA_real_
  else
    as.numeric(ra)
}


#' @rdname StepStimMetrics
#'
#' @param event Event time relative to which latency is measured.
#'
#' @return For \code{step_response_latency()}, the latency in seconds to the
#'   second spike of the first post-event ISI that falls below the threshold
#'   derived from pre-event ISIs, or \code{NA_real_} if no such event is found.
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

