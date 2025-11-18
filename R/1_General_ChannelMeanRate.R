## Generic (created only if not already defined)
setGeneric("ChannelMeanRate", function(X, ...)
  standardGeneric("ChannelMeanRate"))

#' Mean spike rate per channel (EPhysEvents)
#'
#' Counts events per run and channel (using the length of each timestamp vector),
#' divides by recording duration from \code{Metadata$Duration} (or  \code{Metadata$Diff} for backward compatibility), and averages the
#' resulting per-run rates across runs (unweighted mean).
#'
#' @inheritParams Bin
#' @param na.rm Logical; if \code{TRUE} (default), skip NA rates.
#'              (NA rates arise if a channel is missing in a run or \code{Diff} is invalid.)
#' @return A named numeric vector: mean spikes/second per channel.
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData lapply
#' @export
setMethod("ChannelMeanRate", "EPhysEvents", function(X, na.rm = TRUE) {
  md   <- X@Metadata

  if (!"Diff" %in% names(md)) {
    if (!"Duration" %in% names(md)){
      stop("Metadata must contain a numeric column named 'Diff' (legacy) or 'Duration' (recording duration).")
    } else {
      dur <- as.numeric(md$Duration)
    }
  } else {
    dur <- as.numeric(md$Diff)
  }

  counts<-simplify2array(lapply(lapply(X, length),unlist))
  rates <- sweep(counts, 2, dur, "/")
  rates[rates==Inf]<-NA_integer_
  out<-apply(rates, 1, mean, na.rm = na.rm)
  out
})
