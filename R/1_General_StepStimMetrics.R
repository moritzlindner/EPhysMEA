#' Step-stimulus response metrics for `EPhysEvents` (per run × channel)
#'
#' Applies \code{step_stim_resp_metrics()} to every \strong{Channel} of every
#' \strong{run} (trace) in an \code{EPhysEvents} object and returns a tidy
#' \emph{wide} table (one column per metric). Optionally, metrics are averaged
#' across repeated trials with the same \code{RecordingID}.
#'
#' Iteration is performed via \code{EPhysData::lapply(X, ...)}. Any additional
#' arguments (\code{...}) are forwarded to \code{\link{step_stim_resp_metrics}},
#' including the \code{parallel} flag which that function may use internally.
#'
#' @param X An \code{EPhysEvents} instance.
#' @param step_range Numeric or \code{units} length-2 vector \code{c(start, end)}
#'   specifying the step window, interpreted as \emph{[start, end)}. Must satisfy
#'   \code{start < end}. If \code{units}, it is converted to \code{TimeUnits(X)}.
#' @param average_by_recordingID Logical; if \code{TRUE}, average metrics across
#'   runs that share the same \code{RecordingID} in \code{Metadata(X)} using
#'   \code{mean(na.rm = TRUE)}. Default \code{FALSE}.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @inheritDotParams step_stim_resp_metrics
#'   on_win off_win trans_win binwidth smooth gauss_sd return_psth
#'   latency_alpha latency_preduration latency_postduration latency_event
#' @param ... Further arguments forwarded to \code{\link{step_stim_resp_metrics}}.
#'
#' @details
#' Channel order follows \code{Channels(X)}. Row names for the per-run metric
#' matrices are set from \code{Channels(X)} before stacking; this assumes that
#' each run exposes metrics for the same set of channels in a consistent way.
#' If your runs differ in available channels, consider harmonizing beforehand.
#'
#' @return
#' \itemize{
#'   \item If \code{average_by_recordingID = FALSE} (default): a wide
#'         \code{data.frame} with columns \code{RunUID}, \code{Channel},
#'         one column per metric returned by \code{step_stim_resp_metrics()},
#'         and the \code{RecordingID} column joined from \code{Metadata(X)}.
#'   \item If \code{average_by_recordingID = TRUE}: one row per
#'         \code{RecordingID} × \code{Channel} with the same metric columns
#'         averaged across runs and an additional \code{n_trials} column.
#' }
#'
#' @seealso \code{\link{step_stim_resp_metrics}}, \code{EPhysData::Metadata},
#'   \code{EPhysData::Channels}, \code{EPhysData::TimeUnits}
#'
#' @examples
#' # Per-(run, channel) metrics:
#' # df <- StepStimMetrics(ev, step_range = c(1, 2))
#'
#' # Averaged by RecordingID with custom latency alpha:
#' # avg <- StepStimMetrics(ev, step_range = c(1, 2),
#' #                        average_by_recordingID = TRUE,
#' #                        latency_alpha = 0.05)
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels TimeUnits lapply
#' @importFrom units set_units
#' @importFrom reshape2 melt dcast
#' @seealso step_stim_resp_metrics
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
