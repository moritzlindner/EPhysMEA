#' Assign Event Timestamps to Table of Content (Metadata)
#'
#' For each window defined in a Table of Content (ToC) table, extracts the timestamps
#' from each channel in a named list of time-series data that fall within that window,
#' applies optional start/end offsets, and then subtracts the offset start so that times
#' are zero-based.
#'
#' @param ToC A data frame with at least three columns:
#'   \describe{
#'     \item{Start}{Numeric start time of each window.}
#'     \item{Stop}{Numeric stop time of each window.}
#'     \item{RunUID}{Unique identifier for each window (used as output names).}
#'   }
#' @param Data A named list of numeric vectors, where each element represents the
#'   timestamps for one channel.
#'
#' @return A named list of length \code{nrow(ToC)}. Each element corresponds to one
#'   row of \code{ToC} (named by its \code{RunUID}) and is itself a named list
#'   of channels. Within each sublist, each channel entry is a numeric vector of
#'   **zero-based** timestamps (after subtracting \code{start_offset}) that lie in
#'   \code{[Start + start_offset, Stop + end_offset]}. Channels with zero matching
#'   timestamps are omitted.
#'
#' @details
#' The function loops over each row of \code{ToC}, applies the specified offsets
#' to the \code{Start} and \code{Stop} times, subsets each channel’s timestamp
#' vector to that interval, subtracts the adjusted start to zero-base the times,
#' and drops any channels without hits.
#'
#' @examples
#' \dontrun{
#' ToC <- data.frame(
#'   RunUID = c("A1", "B2"),
#'   Start    = c(0.0, 5.0),
#'   Stop     = c(4.9, 9.9),
#'   stringsAsFactors = FALSE
#' )
#' Data <- list(
#'   ch1 = c(-0.1, 0.5, 2.0, 6.5),
#'   ch2 = c(1.5, 5.5, 8.0)
#' )
#' # apply a 0.1s start offset and a 0.2s end offset
#' result <- Assign_Events_To_ToC(ToC, Data)
#' # result[["A1"]]$ch1  # times in [0.1, 5.1] minus 0.1
#' }
#'
#' @export
Assign_Events_To_ToC <- function(ToC,
                                 Data) {
  idx <- seq_len(nrow(ToC))
  res <- lapply(idx, function(i) {
    row <- ToC[i, ]
    adj_start <- row$Start
    adj_stop  <- row$Stop
    ev  <- lapply(Data, function(ts) {
      sub <- ts[ts >= adj_start & ts <= adj_stop] - adj_start
      if (length(sub)) sub else NULL
    })
    lapply(ev, function(x) if (is.null(x)) numeric(0) else x)
  })
  names(res) <- as.character(ToC$RunUID)
  res
}
