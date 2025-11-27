#' Direction tuning for EPhysEvents (multi-run × multi-channel)
#'
#' @description
#' For each run (row of \code{Metadata}) and each channel, this method:
#' \enumerate{
#'   \item computes the mean response per direction window via
#'         \code{\link{direction_response_amplitude}} (one value per row of \code{windows_df});
#'   \item summarizes direction tuning with
#'         \code{\link{direction_selectivity_index}} (gDSI; normalized vector sum) and
#'         \code{\link{direction_preferred_direction}} (preferred direction in degrees).
#' }
#' It returns two tidy data.frames (analogous to \code{Recfield}): a long-form stack of
#' per-direction responses (\code{maps}) and a per-(Run, Channel) table of metrics. Channels without enough events/spikes to reach a resolution for DSI of at least 0.2 are ignored and results set to NA.
#'
#' @details
#' \itemize{
#'   \item \strong{EPhysEvents structure:} \code{X@Data} is a list over runs/trials; each
#'         \code{X@Data[[i]]} is a \emph{named list of channels}, and each channel entry is a
#'         numeric vector of spike times (seconds).
#'   \item \strong{Windows:} \code{windows_df} must contain \code{start}, \code{end},
#'         and \code{direction_deg}. Times may be numeric or \pkg{units} convertible to seconds.
#'         If \code{check_equal_durations=TRUE}, all window durations must be equal within \code{tol}.
#'   \item \strong{Counting convention:} spikes are counted in half-open intervals \code{[start, end)}.
#' }
#'
#' @param X      An \code{EPhysEvents} X.
#' @inheritParams direction_response_amplitude
#' @param check_equal_durations Logical; if \code{TRUE} (default) enforce identical window
#'                    durations across rows of \code{windows_df} within \code{tol}.
#' @param tol         Numeric tolerance for duration equality (seconds). Default \code{1e-9}.
#' @param ReturnMap   Logical; if \code{TRUE} (default) returns the long-form responses table
#'                    in the \code{maps} element.
#' @inheritParams EPhysData::`lapply-EPhys`
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{maps}}{(if \code{ReturnMap}) data.frame with columns
#'         \code{Run}, \code{Channel}, \code{direction_deg}, \code{Response_Amplitude},
#'         \code{start}, \code{end}.}
#'   \item{\code{metrics}}{data.frame with columns
#'         \code{Run}, \code{Channel}, \code{dsi}, \code{preferred_direction}.}
#' }
#'
#' @seealso \code{\link{direction_response_amplitude}},
#'   \code{\link{direction_selectivity_index}},
#'   \code{\link{direction_preferred_direction}},
#'   \code{\link{Recfield}}
#'
#' @examples
#' \donttest{
#' # E: EPhysEvents, dir_win: data.frame(start, end, direction_deg)
#' res <- DirectionTuning(E, windows_df = dir_win)
#' head(res$maps)
#' head(res$metrics)
#' }
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @export
setGeneric("DirectionTuning",
           function(X, windows_df,
                    check_equal_durations = TRUE,
                    tol = 1e-9,
                    ReturnMap = TRUE,
                    parallel = !interactive(),
                    error = c("stop", "warn")[1],
                    progress = interactive()) standardGeneric("DirectionTuning")
)

#' @rdname DirectionTuning
#' @export
setMethod("DirectionTuning", signature(X = "EPhysEvents"),
          function(X,
                   windows_df,
                   check_equal_durations = TRUE,
                   tol = 1e-9,
                   ReturnMap = TRUE,
                   parallel = !interactive(),
                   error = c("stop", "warn")[1],
                   progress = interactive()) {

            # --- window sanity (optional equal durations) ---
            if (isTRUE(check_equal_durations)) {
              starts <- to_num(windows_df$start)
              ends   <- to_num(windows_df$end)
              dur    <- ends - starts
              if (any(!is.finite(dur) | dur <= 0))
                stop("Non-finite or non-positive window duration in 'windows_df'.")
              if (length(dur) > 1 && any(abs(dur - dur[1]) > tol)) {
                stop(sprintf("Direction windows must have equal duration (tol=%g). Got: %s",
                             tol, paste(round(dur, 6), collapse = ", ")))
              }
            }

            # --- per-(run×channel) worker applied by lapply(EPhysEvents, ...) ---
            single_worker <- function(spikes, windows_df, ReturnMap) {
              # responses per window (half-open [start, end))
              ra <- direction_response_amplitude(spikes, windows_df)

              # long-form per-direction responses
              map_df <- data.frame(
                direction_deg      = windows_df$direction_deg,
                Response_Amplitude = as.numeric(ra),
                start              = windows_df$start,
                end                = windows_df$end,
                check.names = FALSE
              )

              # metrics from the map
              if (sum(map_df$Response_Amplitude*(map_df$end-windows_df$start))>=10){ # require a minimum total of 10 spikes within the windows yielding a resolution of 0.2
                dsi <- direction_selectivity_index(map_df,
                                                   resp  = "Response_Amplitude",
                                                   angle = "direction_deg")
                pd  <- direction_preferred_direction(map_df,
                                                     resp  = "Response_Amplitude",
                                                     angle = "direction_deg")
              } else {
                dsi <- NA
                pd <- NA
              }


              met_df <- data.frame(
                dsi                 = as.numeric(dsi),
                preferred_direction = as.numeric(pd),
                check.names = FALSE
              )

              if (isTRUE(ReturnMap)) {
                list(metrics = met_df, map = map_df)
              } else {
                list(metrics = met_df)
              }
            }

            # --- run×channel application via your EPhysEvents lapply ---
            res_list <- lapply(
              X,
              single_worker,
              windows_df,
              ReturnMap,
              parallel = parallel,
              error    = error,
              progress = progress
            )

            # --- collect metrics ---
            metrics_nested <- lapply(res_list, function(x) lapply(x, function(y) as.data.frame(y$metrics)))
            metrics_df     <- nested2df(X, metrics_nested)

            if (isTRUE(ReturnMap)) {
              maps_nested <- lapply(res_list, function(x) lapply(x, function(y) as.data.frame(y$map)))
              maps_df     <- nested2df(X, maps_nested)
              return(list(maps = maps_df,
                          metrics = metrics_df))
            } else {
              return(list(metrics = metrics_df))
            }
          })
