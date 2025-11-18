#' Compute inter‐event‐intervals (IEI) from an EPhysEvents object
#'
#' Transforms each channel’s timestamp vectors in an \code{EPhysEvents} object
#' into a two‐column data.frame of \code{Time} and \code{IEI}, and returns an
#' \code{EPhysIEI} object.
#'
#' @inheritParams Bin
#' @return An \code{EPhysIEI} object with:
#' \describe{
#'   \item{Metadata}{Copied from \code{x@Metadata}.}
#'   \item{Data}{A list of length \code{nrow(x@Metadata)}, each element a named list of
#'     channels. For each channel, a data.frame with columns:
#'     \itemize{
#'       \item \code{Time}: original zero‐based timestamps,
#'       \item \code{IEI}: inter‐event intervals (\code{c(NA, diff(Time))}).
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' rec_ev <- LoadEPhysEvents("data/session1", "types.tsv")
#' rec_iei <- ComputeIEI(rec_ev)
#' str(rec_iei@Data[[1]]$ChannelA)
#' }
#' @importClassesFrom EPhysData EPhysEvents
#' @export
setGeneric("ComputeIEI", function(x) standardGeneric("ComputeIEI"))

#' @describeIn ComputeIEI Method for EPhysEvents
#' @export
setMethod("ComputeIEI", "EPhysEvents", function(x) {
  iei_list <- EPhysData::lapply(x, function(ts) {
    if (length(ts) > 0L) {
      cbind(
        Time = ts,
        IEI  = c(NA_real_, diff(ts))
      )
    } else {
      cbind(
        Time = NA_real_,
        IEI  = NA_real_
      )
    }
  })

  new(
    "EPhysIEI",
    Metadata      = x@Metadata,
    Data          = iei_list,
    ExamInfo      = x@ExamInfo,
    SubjectInfo   = x@SubjectInfo,
    Imported      = x@Imported,
    TimeTrace     = numeric(0),
    Channels      = x@Channels,
    StimulusTrace = x@StimulusTrace,
    TimeUnits     = x@TimeUnits,
    StimulusUnits = x@StimulusUnits
  )
})
