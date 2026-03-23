#' Mean response amplitude per direction window (full-window rate)
#'
#' Computes, for each stimulus window \code{[start, end)}, the mean firing
#' rate as \code{count(start <= t < end) / (end - start)}. This function
#' does **not** perform latency detection, ON/OFF splitting, or baseline
#' subtraction; it is a straight full-window rate using
#' \code{\link{step_spikerate}}.
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
#' @seealso \code{\link{step_spikerate}},
#'   \code{\link{motion_dsi}},
#'   \code{\link{motion_preferred_direction}}
#'
#' @examples
#' # Toy spikes and two direction windows
#' spikes <- c(0.10, 0.15, 0.51, 0.70, 0.71)
#' windows_df <- data.frame(
#'   start = c(0.00, 0.50),
#'   end   = c(0.40, 0.90),
#'   direction_deg = c(0, 180)
#' )
#' motion_direction_amplitude(spikes, windows_df)
#'
#' @export
motion_direction_amplitude <- function(
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
      amp_win <- step_spikerate(spikes_win, width_win)
    }
    amp_win
  }

  out <- do.call(c, lapply(seq_len(nrow(windows_df)), do_one))
  out
}
#' Motion tuning metrics: gDSI (vector strength), gOSI, preferred direction, and area score
#'
#' Computes motion direction selectivity as the normalized vector sum (gDSI / vector strength),
#' the preferred direction (vector-sum angle), motion orientation (axis) selectivity as the
#' doubled-angle normalized vector sum (gOSI), and an area-based tuning score
#' (\code{motion_tuning_area_score()}).
#'
#' Let \eqn{R_d \ge 0} be the response magnitude for direction \eqn{\theta_d} (degrees).
#' Responses are rectified internally via \code{pmax(0, \cdot)}; if baseline subtraction is needed,
#' supply baseline-subtracted values in \code{resp}.
#'
#' \strong{Direction selectivity (gDSI):}
#' \deqn{\mathrm{gDSI}=\frac{\left|\sum_d R_d e^{i\theta_d}\right|}{\sum_d R_d}
#'      =\frac{\sqrt{\left(\sum_d R_d\cos\theta_d\right)^2+\left(\sum_d R_d\sin\theta_d\right)^2}}{\sum_d R_d}}
#'
#' Preferred direction is \eqn{\mathrm{atan2}(V_y,V_x)} in degrees (0–360), where
#' \eqn{V_x=\sum_d R_d\cos\theta_d} and \eqn{V_y=\sum_d R_d\sin\theta_d}.
#'
#' \strong{Orientation (axis) selectivity (gOSI):}
#' \deqn{\mathrm{gOSI}=\frac{\left|\sum_d R_d e^{i2\theta_d}\right|}{\sum_d R_d}}
#' This quantity is insensitive to a 180° flip (i.e., \eqn{\theta} vs \eqn{\theta+\pi})
#' and is therefore appropriate for bidirectional (axis/orientation) tuning.
#'
#' \strong{Area-based motion tuning score:}
#' \code{motion_tuning_area_score()} treats the per-direction responses as radii in polar
#' coordinates, forms the corresponding polygon in the plane, and returns
#' \deqn{1 - A/A_{\max}}
#' where \eqn{A} is the polygon area after max-normalizing responses to 1 and \eqn{A_{\max}} is the
#' area obtained for the same sampled angles when all radii equal 1. The score is in \eqn{[0,1]},
#' with values near 0 indicating broad/untuned responses and values near 1 indicating sharply
#' tuned responses.
#'
#' @param df    Data frame containing at least the response and angle columns.
#' @param resp  Name of the response column in \code{df}. Default \code{"Response_Amplitude"}.
#' @param angle Name of the direction column in degrees (0–360). Default \code{"direction_deg"}.
#'
#' @return
#' \itemize{
#'   \item \code{motion_dsi()}: a single numeric (gDSI in [0,1]) or \code{NA_real_} if
#'         responses are missing or sum to zero.
#'   \item \code{motion_preferred_direction()}: a single numeric preferred direction in degrees
#'         (0–360), or \code{NA_real_} if undefined.
#'   \item \code{motion_osi()}: a single numeric (gOSI in [0,1]) or \code{NA_real_} if
#'         responses are missing or sum to zero.
#'   \item \code{motion_tuning_area_score()}: an area-based score in [0,1] (see Details).
#' }
#'
#' @references
#' Direction selectivity (gDSI): \doi{10.1016/j.neuron.2021.07.008}.
#'
#' Orientation (axis) selectivity (gOSI using doubled angles): \doi{10.1038/s41467-025-64321-1}.
#' (Note: that paper also defines an alternative OSI as \eqn{(R_{PO}-R_{NO})/(R_{PO}+R_{NO})};
#' this package implements the doubled-angle gOSI.)
#'
#' @examples
#' df <- data.frame(direction_deg = seq(0, 315, by = 45),
#'                  Response_Amplitude = c(0,1,3,5,4,2,1,0))
#' motion_preferred_direction(df)
#' motion_dsi(df)
#' motion_osi(df)
#' motion_tuning_area_score(df)
#'
#' @name motion_tuning_metrics
NULL



#' @rdname motion_tuning_metrics
#' @export
motion_dsi <- function(df, resp = "Response_Amplitude", angle = "direction_deg") {
  th <- (df[[angle]] %% 360) * pi/180
  R  <- pmax(0, as.numeric(df[[resp]]))
  ok <- is.finite(th) & is.finite(R)
  th <- th[ok]; R <- R[ok]
  if (!length(R) || sum(R) == 0) return(NA_real_)
  C <- sum(R * cos(th)); S <- sum(R * sin(th))
  sqrt(C*C + S*S) / sum(R)
}

#' @rdname motion_tuning_metrics
#' @export
motion_preferred_direction <- function(df, resp = "Response_Amplitude", angle = "direction_deg") {
  th <- (df[[angle]] %% 360) * pi/180
  R  <- pmax(0, as.numeric(df[[resp]]))
  ok <- is.finite(th) & is.finite(R)
  th <- th[ok]; R <- R[ok]
  if (!length(R) || sum(R) == 0) return(NA_real_)
  pd_deg <- (atan2(sum(R * sin(th)), sum(R * cos(th))) * 180/pi) %% 360
  as.numeric(pd_deg)
}

#' @rdname motion_tuning_metrics
#' @export
motion_osi <- function(df, resp = "Response_Amplitude", angle = "direction_deg") {
  th <- (df[[angle]] %% 360) * pi/180
  R  <- pmax(0, as.numeric(df[[resp]]))
  ok <- is.finite(th) & is.finite(R)
  th <- th[ok]; R <- R[ok]
  if (!length(R) || sum(R) == 0) return(NA_real_)
  C2 <- sum(R * cos(2 * th)); S2 <- sum(R * sin(2 * th))
  sqrt(C2*C2 + S2*S2) / sum(R)
}


#' @rdname motion_tuning_metrics
#' @export
motion_tuning_area_score <- function(df,
                            resp = "Response_Amplitude",
                            angle = "direction_deg") {

  ang_deg <- as.numeric(df[[angle]]) %% 360
  R <- as.numeric(df[[resp]])

  ok <- is.finite(ang_deg) & is.finite(R)
  ang_deg <- ang_deg[ok]
  R <- R[ok]

  if (!length(R)) return(NA_real_)

  R <- pmax(0, R)

  # average duplicates (by exact degree value)
  Ru <- tapply(R, ang_deg, mean, na.rm = TRUE)
  ang_u_deg <- as.numeric(names(Ru))
  R_u <- as.numeric(Ru)

  # need at least 3 vertices for a polygon
  if (length(R_u) < 3L) {
    out <- NA_real_
    return(out)
  }

  if (sum(R_u) == 0) {
    out <- NA_real_
    return(out)
  }

  mx <- max(R_u)
  if (!is.finite(mx) || mx <= 0) {
    out <- NA_real_
    return(out)
  }
  R_u <- R_u / mx

  # sort by angle to avoid self-crossing due to arbitrary ordering
  o <- order(ang_u_deg)
  ang_u_deg <- ang_u_deg[o]
  R_u <- R_u[o]

  th <- ang_u_deg * pi/180

  poly_area <- function(x, y) {
    # shoelace; assumes vertices ordered; closes polygon internally
    x2 <- c(x, x[1]); y2 <- c(y, y[1])
    abs(sum(x2[-length(x2)] * y2[-1] - x2[-1] * y2[-length(y2)])) / 2
  }

  # observed polygon (r = R_u)
  x <- R_u * cos(th)
  y <- R_u * sin(th)
  A <- poly_area(x, y)

  # maximal polygon for this angular sampling (r = 1)
  x1 <- cos(th)
  y1 <- sin(th)
  Amax <- poly_area(x1, y1)

  if (!is.finite(A) || !is.finite(Amax) || Amax <= 0) {
    out <- NA_real_
    return(out)
  }

  score <- 1 - (A / Amax)
  # numeric safety
  score <- max(0, min(1, score))

  as.numeric(score)
}
