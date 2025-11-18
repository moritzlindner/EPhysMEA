#' Add stimulus trace to EPhysContinuous by binning stimulus samples
#'
#' Given a continuous stimulus vector sampled at \code{sample_rate} (with units),
#' aggregates the stimulus into the time bins defined by \code{X@TimeTrace} and
#' stores the result (with units, if present) in the \code{StimulusTrace} slot,
#' updating \code{StimulusUnits} accordingly.
#'
#' @param X An \code{EPhysEvents} or \code{EPhysContinuous} object.
#' @param stim_values  Numeric or \pkg{units} vector of stimulus samples.
#' @param sample_rate  A \pkg{units} scalar specifying the sampling rate;
#'                     it will be converted to \code{X@TimeUnits} before use.
#' @param FUN          Aggregation function applied to samples in each bin
#'                     (default \code{mean}).
#' @param ...          Additional arguments passed to \code{FUN}.
#'
#' @return The same \code{EPhysContinuous} object, with:
#' \describe{
#'   \item{StimulusTrace}{Numeric or \pkg{units} vector of length \code{length(X@TimeTrace)}.}
#'   \item{StimulusUnits}{Updated to the unit of \code{stim_values} if it was a \pkg{units} object.}
#' }
#'
#' @importFrom units set_units deparse_unit as_units
#' @importClassesFrom EPhysData EPhysContainer
#' @importFrom EPhysData StimulusTrace TimeTrace StimulusUnits TimeUnits StimulusUnits<- StimulusTrace<-
#' @name AddStimulus
setGeneric(
  name = "AddStimulus",
  def = function(X, stim_values, sample_rate, FUN = mean, ...)
  {
    standardGeneric("AddStimulus")
  }
)
#' @export
setMethod("AddStimulus",
          signature(X = "EPhysContainer"),
          function(X, stim_values, sample_rate, FUN = mean, ...) {
            # convert sample_rate into object's time units
            if(!("units" %in% class(sample_rate))){
              warning("Units of sample_rate will be assumed to be the same as Time Trace.")
            }

            if (length(TimeTrace(X)) < 2L) stop("TimeTrace must have at least two points.")

            sample_rate <- 1/set_units(sample_rate, TimeUnits(X), mode = "standard")
            sr_num <- as.numeric(sample_rate)

            # handle stimulus units
            if (inherits(stim_values, "units")) {
              stim_unit   <- deparse_unit(stim_values)
              stim_numeric <- as.numeric(stim_values)
            } else {
              stim_unit   <- NULL
              stim_numeric <- stim_values
            }

            # generate sample times
            n <- length(stim_numeric)
            sample_times <- seq(0, by = 1/sr_num, length.out = n)

            # define bins from TimeTrace
            centers <- drop_units(TimeTrace(X))
            if (length(centers) < 2) {
              stop("TimeTrace must have at least two points to define bins.")
            }
            bw      <- centers[2] - centers[1]
            half_bw <- bw / 2
            edges   <- c(centers - half_bw, centers[length(centers)] + half_bw)

            # bin stimulus
            binned_vals <- vapply(seq_along(centers), function(i) {
              idx <- which(sample_times >= edges[i] & sample_times < edges[i+1])
              if (length(idx) > 0) FUN(stim_numeric[idx], ...) else NA_real_
            }, numeric(1))

            # re‐apply units if needed
            if (!is.null(stim_unit)) {
                StimulusUnits(X) <- stim_unit
            }

            StimulusTrace(X) <- binned_vals
            validObject(X)
            X
          })
