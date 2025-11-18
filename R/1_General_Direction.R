#' Mean response amplitude per direction window (full-window rate)
#'
#' Computes, for each stimulus window \code{[start, end)}, the mean firing
#' rate as \code{count(start <= t < end) / (end - start)}. This function
#' does **not** perform latency detection, ON/OFF splitting, or baseline
#' subtraction; it is a straight full-window rate using
#' \code{\link{step_response_amplitude}}.
#'
#' @param spikes Numeric (or units) vector of spike times in seconds.
#' @param windows_df A data.frame with at least the columns:
#'   \describe{
#'     \item{\code{start}}{Window start time (s).}
#'     \item{\code{end}}{Window end time (s) (right-open: spikes at \code{end} are excluded).}
#'     \item{\code{direction_deg}}{Direction label (degrees). Not used in the computation,
#'           kept for downstream tuning metrics.}
#'   }
#'   Optionally \code{window_id}. Time columns may be numeric or units
#'   convertible to seconds.
#'
#' @details
#' Spike counting uses half-open intervals \code{[start, end)}. Windows with
#' non-finite or non-positive duration yield \code{NA}. This function assumes
#' a helper \code{to_num()} is available to coerce inputs to numeric seconds
#' (e.g., from \pkg{units}); if not, define it or provide pure numeric inputs.
#'
#' @return
#' A numeric vector of length \code{nrow(windows_df)} with mean firing rates
#' (spikes/s) per window; entries may be \code{NA_real_} for invalid windows.
#'
#' @seealso \code{\link{step_response_amplitude}},
#'   \code{\link{direction_selectivity_index}},
#'   \code{\link{direction_preferred_direction}}
#'
#' @examples
#' # Toy spikes and two direction windows
#' spikes <- c(0.10, 0.15, 0.51, 0.70, 0.71)
#' windows_df <- data.frame(
#'   start = c(0.00, 0.50),
#'   end   = c(0.40, 0.90),
#'   direction_deg = c(0, 180)
#' )
#' direction_response_amplitude(spikes, windows_df)
#'
#' @export
direction_response_amplitude <- function(
    spikes,
    windows_df
) {
  spikes <- sort(to_num(spikes))
  stopifnot(is.numeric(spikes))

  req <- c("start","end","direction_deg")
  if (!all(req %in% names(windows_df)))
    stop("windows_df must contain columns: ", paste(req, collapse=", "))

  do_one <- function(i) {
    s   <- to_num(windows_df$start[i])
    e   <- to_num(windows_df$end[i])

    if (!is.finite(s) || !is.finite(e) || e <= s) {
      return(data.frame(
        window_id = (if ("window_id" %in% names(windows_df)) windows_df$window_id[i] else i),
        direction_deg = windows_df$direction_deg[i],
        start = s, end = e, mid = NA_real_,
        baseline_rate = NA_real_,
        Response_Amplitude = NA_real_
      ))
    }
    amp_win <- NA_real_
    width_win <- e - s
    if (width_win > 0) {
      spikes_win <- spikes[spikes >= s & spikes < e]
      amp_win <- step_response_amplitude(spikes_win, width_win)
    }
    amp_win
  }

  out <- do.call(c, lapply(seq_len(nrow(windows_df)), do_one))
  out
}

#' Direction tuning metrics: gDSI (vector strength) and preferred direction
#'
#' Computes direction selectivity as the normalized vector sum (a.k.a. gDSI)
#' and the preferred direction from per-direction responses.
#'
#' Let \eqn{R_d \ge 0} be the response magnitude for direction \eqn{\theta_d}
#' (degrees). The resultant vector is
#' \deqn{\mathbf{V}=\left(\sum_d R_d\cos\theta_d,\ \sum_d R_d\sin\theta_d\right)}
#' and the gDSI / vector strength is \eqn{\|\mathbf{V}\| / \sum_d R_d}.
#' The preferred direction is \eqn{\mathrm{atan2}(V_y,V_x)} in degrees.
#'
#' @param df   Data frame containing at least the response and angle columns.
#' @param resp Name of the response column in `df`. Values are baseline-subtracted
#'             and rectified internally via `pmax(0, ·)`. Default `"Response_Amplitude"`.
#' @param angle Name of the direction column in degrees (0–360). Default `"direction_deg"`.
#'
#' @return
#' - `direction_selectivity_index()`: a single numeric (gDSI in [0,1]) or `NA_real_`
#'   if responses are missing or sum to zero.
#' - `direction_preferred_direction()`: a single numeric preferred direction in degrees
#'   (0–360), or `NA_real_` if undefined.
#'
#' @examples
#' df <- data.frame(direction_deg = seq(0, 315, by = 45),
#'                  Response_Amplitude = c(0,1,3,5,4,2,1,0))
#' direction_preferred_direction(df)
#' direction_selectivity_index(df)  # same as vector strength
#'
#' @name direction_tuning_metrics
NULL


#' @rdname direction_tuning_metrics
#' @export
direction_selectivity_index <- function(df, resp = "Response_Amplitude", angle = "direction_deg") {
  th <- (df[[angle]] %% 360) * pi/180
  R  <- pmax(0, as.numeric(df[[resp]]))
  ok <- is.finite(th) & is.finite(R)
  th <- th[ok]; R <- R[ok]
  if (!length(R) || sum(R) == 0) return(NA_real_)
  C <- sum(R * cos(th)); S <- sum(R * sin(th))
  sqrt(C*C + S*S) / sum(R)
}

#' @rdname direction_tuning_metrics
#' @export
direction_preferred_direction <- function(df, resp = "Response_Amplitude", angle = "direction_deg") {
  th <- (df[[angle]] %% 360) * pi/180
  R  <- pmax(0, as.numeric(df[[resp]]))
  ok <- is.finite(th) & is.finite(R)
  th <- th[ok]; R <- R[ok]
  if (!length(R) || sum(R) == 0) return(NA_real_)
  pd_deg <- (atan2(sum(R * sin(th)), sum(R * cos(th))) * 180/pi) %% 360
  as.numeric(pd_deg)
}
