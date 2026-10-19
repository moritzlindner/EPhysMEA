#' Count events per metadata row and channel
#'
#' Counts the number of events in each run (metadata row) and channel of an
#' \code{\link[EPhysData:EPhysEvents-class]{EPhysEvents}} object.
#'
#' \code{Count()} is intentionally a simple event-counting operation. It does
#' not divide by recording duration, average across repetitions, bin the data,
#' or perform baseline subtraction. To count events in a defined time window,
#' first restrict the object with \code{EPhysData::Subset(..., TimeRange = ...)}
#' and then call \code{Count()}.
#'
#' @param X An \code{EPhysEvents} object.
#'
#' @return An integer matrix with one row per row of \code{Metadata(X)} and one
#'   column per entry of \code{Channels(X)}. Row names are the corresponding
#'   \code{RunUID} values. Empty spike trains have count 0.
#'
#' @examples
#' \dontrun{
#' # Count all events in every run and channel
#' counts <- Count(events)
#'
#' # Count events in a defined analysis window
#' on <- EPhysData::Subset(events, TimeRange = c(2.0, 2.4))
#' n_on <- Count(on)
#'
#' # For a single selected channel, calculate its mean firing rate
#' # across repetitions:
#' channel <- EPhysData::Channels(on)[1]
#' r_on <- mean(n_on[, channel]) / 0.4
#' }
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @export
setGeneric("Count", function(X) standardGeneric("Count"))

#' @describeIn Count Count events for each metadata row and channel of an EPhysEvents object.
#' @export
setMethod("Count", "EPhysEvents", function(X) {
  md <- Metadata(X)
  channels <- Channels(X)

  counts <- matrix(
    0L,
    nrow = nrow(md),
    ncol = length(channels),
    dimnames = list(as.character(md$RunUID), channels)
  )

  if (!nrow(md) || !length(channels)) {
    return(counts)
  }

  for (i in seq_len(nrow(md))) {
    run <- X@Data[[i]]

    counts[i, ] <- vapply(seq_along(channels), function(j) {
      ch <- channels[[j]]

      if (!is.null(names(run)) && ch %in% names(run)) {
        return(as.integer(length(run[[ch]])))
      }

      if (j <= length(run)) {
        return(as.integer(length(run[[j]])))
      }

      0L
    }, integer(1L))
  }

  counts
})
