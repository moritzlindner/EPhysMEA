#' Bin spikes into an EPhysContinuous object
#'
#' Bins per‐trial spike times from an \code{EPhysEvents} object into a
#' left‐inclusive, right‐exclusive 3D numeric count array
#' (\code{time × trial × channel}), storing counts and metadata.
#'
#' @param X An \code{EPhysEvents} object.
#' @param binWidth A \pkg{units} object specifying the width of each time bin;
#'                 its unit will be set to \code{TimeUnits(X)}.
#' @param StimulusBinFunction Function to bin stimulus trace by (if present).
#' @inheritParams EPhysData::`lapply-EPhys`
#'
#' @return An \code{EPhysContinuous} object with:
#' \describe{
#'   \item{Data}{A 3D numeric array \code{[time × trial × channel]} of spike counts.}
#'   \item{TimeTrace}{Numeric vector of bin starts (in \code{TimeUnits(X)}).}
#'   \item{Channels}{Character vector of channel names.}
#'   \item{StimulusTrace}{Empty numeric vector.}
#'   \item{TimeUnits}{Copied from \code{TimeUnits(X)}.}
#'   \item{Metadata}{Copied from \code{Metadata(X)}.}
#' }
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels ChannelMetadata TimeTrace StimulusTrace TimeUnits StimulusUnits
#' @importFrom units as_units
setGeneric(
  name = "Bin",
  def = function(X,
                 binWidth,
                 StimulusBinFunction = mean, parallel = !interactive())
  {
    standardGeneric("Bin")
  }
)
#' @export
setMethod("Bin",
          signature(X = "EPhysEvents"),
          function(X, binWidth, StimulusBinFunction = mean, parallel = !interactive()) {

            ## --- validate bin width and StimulusBinFunction ---
            set_units(binWidth, TimeUnits(X), mode = "standard")
            bw <- as.numeric(binWidth)

            stimbinfx <- StimulusBinFunction
            tryCatch({
              stopifnot(is.function(stimbinfx))
              outfx <- stimbinfx(1:10)
              stopifnot(is.numeric(outfx), length(outfx) == 1L)
            }, error = function(e) {
              stop("StimulusBinFunction must return a single numeric for a numeric StimulusTrace subset.")
            })

            ## --- metadata and basic checks ---
            md       <- Metadata(X)
            data     <- X@Data
            channels <- Channels(X)

            Experiment <- unique(md$Experiment)
            if (length(Experiment) != 1L) {
              stop("Expected unique 'Experiment', found: ", paste(Experiment, collapse = ", "))
            }
            durs <- md$Diff
            if (!all(durs == durs[1])) stop("Run durations (Diff) must be identical across all runs")
            duration <- durs[1]

            ## --- define bins ---
            breaks <- seq(0, duration, by = bw)
            if (tail(breaks, 1) < duration) breaks <- c(breaks, duration)
            nBins   <- length(breaks) - 1L
            starts  <- breaks[-length(breaks)]
            bin_names <- paste0("b", seq_len(nBins))   # stable column names for nested2df

            nRuns <- length(md$RunUID)
            nCh   <- length(channels)

            ## --- per (run × channel) counting via lapply(EPhysEvents, ...) ---
            # returns named numeric vector (length = nBins) per channel; zeros if channel missing/empty
            res_nested <- lapply(
              X  = X,
              FUN = function(ts) {
                if (is.null(ts) || !length(ts)) {
                  return(stats::setNames(rep(0L, nBins), bin_names))
                }
                ts <- ts[is.finite(ts)]
                ts <- ts[ts >= 0 & ts < duration]
                cnt <- hist(ts, breaks = breaks, plot = FALSE, right = FALSE)$counts
              },
              , parallel = parallel # (use defaults: error="stop", progress = interactive())
            )

            ## --- reshape nested output with your helper ---
            # nested2df(X, res_nested) must produce: RunUID, RecordingID, Channel, <bin columns>
            ## --- reshape nested output without nested2df: build 3D array [bins × runs × channels] ---
            cnts <- array(
              0,
              dim = c(nBins, nRuns, nCh)
            )

            for (i in seq_len(nRuns)) {
              run_res <- res_nested[[i]]

              # Build [bins × channels] for this run.
              M <- vapply(seq_len(nCh), function(j) {
                v <- run_res[[j]]
                if (is.null(v)) {
                  rep(0, nBins)
                } else {
                  v <- as.numeric(v)
                  v
                }
              }, numeric(nBins))

              # Place into array slice for this run
              cnts[, i,] <- M
            }

            ## --- Re-binning stimulus with LOCF for empty bins ---
            stim_out       <- StimulusRebin(X, breaks, bin_fun = stimbinfx)
            stim_units_out <- if (length(stim_out)) StimulusUnits(X) else ""

            ## --- construct EPhysContinuous ---
            newEPhysContinuous(
              Data             = cnts,
              TimeTrace        = starts,
              Metadata         = md,
              Channels         = channels,
              Channel_Metadata = X@Channel_Metadata,
              StimulusTrace    = stim_out,
              TimeUnits        = TimeUnits(X),
              StimulusUnits    = stim_units_out,
              ExamInfo         = X@ExamInfo,
              SubjectInfo      = X@SubjectInfo,
              Imported         = X@Imported
            )
          })


#' Re-bin a stimulus trace onto target bin edges (LOCF for empty bins)
#'
#' Aggregates a regularly or irregularly sampled stimulus trace onto the
#' left-closed, right-open bins defined by `breaks`. If a bin contains no
#' samples, the last observed value before the bin end is carried forward
#' (LOCF). If no prior sample exists for an early bin, the first available
#' value is used.
#'
#' This is a small internal worker used by \code{StimulusRebin()}.
#'
#' @param orig_time Numeric vector of original stimulus sample times.
#' @param orig_vals Numeric vector of stimulus values (same length as `orig_time`).
#' @param breaks    Numeric vector of bin edges (strictly increasing).
#' @param bin_fun   Function to aggregate samples within a bin (must return a scalar).
#'
#' @return A numeric vector of length `length(breaks) - 1L` with one value per bin.
#' @keywords internal
#' @noRd
.rebin_stimulus_locf <- function(orig_time, orig_vals, breaks, bin_fun = mean) {
  stopifnot(is.numeric(orig_time), is.numeric(orig_vals), length(orig_time) == length(orig_vals))
  stopifnot(is.numeric(breaks), length(breaks) >= 2L)

  # Ensure increasing order
  o <- order(orig_time)
  t <- orig_time[o]
  v <- orig_vals[o]

  n_bins <- length(breaks) - 1L
  edges  <- breaks

  vapply(seq_len(n_bins), function(k) {
    left  <- edges[k]
    right <- edges[k + 1]
    idx   <- which(t >= left & t < right)
    if (length(idx) > 0L) {
      val <- bin_fun(v[idx])
      if (length(val) != 1L || !is.finite(val)) stop("`bin_fun` must return a single finite value.")
      val
    } else {
      pos <- max(which(t < right), na.rm = TRUE)
      if (is.finite(pos)) v[pos] else v[1]
    }
  }, numeric(1))
}

#' Stimulus re-binning for EPhysContainer
#'
#' Re-bins the stimulus trace stored in an \code{EPhysContainer} (or subclass)
#' onto the bin edges provided in \code{breaks} using a per-bin aggregation
#' function (default: \code{mean}). Bins are treated as left-closed, right-open.
#'
#' If either \code{TimeTrace} or \code{StimulusTrace} is missing (length 0),
#' or their lengths differ, this method returns \code{numeric(0)}.
#'
#' @param X       An \code{EPhysContainer}.
#' @param breaks  Numeric vector of bin edges (strictly increasing).
#' @param bin_fun Function to aggregate samples within a bin (must return a single value).
#'
#' @return A numeric vector of length \code{length(breaks) - 1} with one value per bin,
#'   or \code{numeric(0)} if re-binning cannot be performed.
#' @examples
#' \dontrun{
#'   stim_binned <- StimulusRebin(X, breaks, mean)
#' }
#' @export
setGeneric("StimulusRebin", function(X, breaks, bin_fun = mean) standardGeneric("StimulusRebin"))

#' @rdname StimulusRebin
setMethod("StimulusRebin", "EPhysContainer",
          function(X, breaks, bin_fun = mean) {
            # Quick capability checks
            breaks<-to_num(breaks)
            has_tt   <- length(TimeTrace(X)) > 0L
            has_stim <- length(StimulusTrace(X)) > 0L
            if (!(has_tt && has_stim)) return(numeric(0))
            if (length(StimulusTrace(X)) != length(TimeTrace(X))) return(numeric(0))

            # Validate bin_fun quickly
            tryCatch({
              test <- bin_fun(c(1, 2, 3))
              if (length(test) != 1L || !is.finite(test)) stop()
            }, error = function(e) {
              stop("`bin_fun` must return a single finite value for a numeric vector.")
            })

            .rebin_stimulus_locf(
              orig_time = to_num(TimeTrace(X)),
              orig_vals = StimulusTrace(X),
              breaks    = breaks,
              bin_fun   = bin_fun
            )
          })

