#' Phase-based response metrics for `EPhysEvents`
#'
#' Compute per-segment phase-locking or phase-modulation metrics from an
#' \code{\link[EPhysData:EPhysEvents-class]{EPhysEvents}} object and return them
#' as an \code{\link[EPhysData:EPhysContinuous-class]{EPhysContinuous}} object.
#'
#' For each run and channel, spike times are assigned a stimulus phase derived
#' from \code{StimulusTrace} using
#' \code{\link{phasemetrics_stimulus_preprocess}} together with
#' \code{\link{spike_phase_from_stimulus}}. Spikes are then grouped into
#' user-defined temporal segments, and one phase metric is computed separately
#' for each segment. The resulting output is an \code{EPhysContinuous} whose time
#' axis corresponds to the numeric segment labels rather than physical time.
#'
#' This page documents both the high-level \code{PhaseMetrics()} method and the
#' metric definitions used internally.
#'
#' @param X An \code{EPhysEvents} instance with non-empty
#'   \code{\link[EPhysData:TimeTrace]{TimeTrace}} and
#'   \code{\link[EPhysData:StimulusTrace]{StimulusTrace}}.
#' @param method Character string specifying which phase metric to compute. One of:
#'   \describe{
#'     \item{\code{"rayleighZ"}}{Rayleigh phase-locking strength
#'     \eqn{Z = n R^2}, computed within a peak-centered phase window for each
#'     segment. Higher values indicate stronger concentration of spike phases.}
#'
#'     \item{\code{"peak_trough"}}{Difference between the largest and smallest
#'     bin count of the circular phase histogram within a segment. Larger values
#'     indicate stronger phase modulation.}
#'
#'     \item{\code{"opponent180"}}{Spike count in a peak-centered phase window
#'     minus the spike count in the corresponding opponent window centered
#'     180\eqn{^\circ} away.}
#'
#'     \item{\code{"opponent90"}}{Spike count in a peak-centered phase window
#'     minus the smaller of the two orthogonal comparison windows centered at
#'     peak \eqn{\pm 90^\circ}.}
#'
#'     \item{\code{"opponent_best"}}{The larger of \code{opponent180} and
#'     \code{opponent90}.}
#'
#'     \item{\code{"isi50"}}{A compactness metric defined as
#'     \eqn{180^\circ} minus the shortest circular arc containing half of all
#'     spike phases within the segment. Larger values indicate tighter phase
#'     clustering.}
#'   }
#' @param seg Defines the segmentation of time into output bins:
#'   \describe{
#'     \item{\code{NULL}}{Default. Segments are defined by
#'     \code{floor(TimeTrace(X))}.}
#'     \item{function}{A function mapping spike times to numeric segment labels.}
#'     \item{numeric vector}{A numeric vector aligned to
#'     \code{TimeTrace(X)} giving one segment label per time sample; spike-wise
#'     segment labels are then obtained by step-function interpolation.}
#'   }
#' @param bin_width_deg Width of the circular phase-histogram bins in degrees.
#'   Used by the histogram-based metrics. Default: \code{15}.
#' @param window_half_deg Half-width in degrees of the peak-centered and opponent
#'   comparison windows used by the windowed phase metrics. Default: \code{45}.
#' @param flatten_trials Logical; if \code{TRUE} (default), first call
#'   \code{\link[EPhysData:FlattenTrials]{FlattenTrials}} so each output trial
#'   corresponds to one \code{RecordingID}. If \code{FALSE}, metrics are
#'   computed for each input run separately.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @param ... Further arguments forwarded to
#'   \code{\link{phasemetrics_stimulus_preprocess}}, for example phase
#'   preprocessing options such as \code{upsample}, \code{level}, or
#'   \code{direction}.
#'
#' @details
#' \strong{Workflow:}
#' \enumerate{
#'   \item Optionally flatten repeated trials by \code{RecordingID} using
#'   \code{\link[EPhysData:FlattenTrials]{FlattenTrials}}.
#'   \item Derive a continuous phase trace from the stimulus with
#'   \code{\link{phasemetrics_stimulus_preprocess}}.
#'   \item Map each spike to a stimulus phase using
#'   \code{\link{spike_phase_from_stimulus}}.
#'   \item Assign each spike to a segment.
#'   \item Compute the selected phase metric separately for each segment and
#'   each channel.
#'   \item Store the results in an \code{EPhysContinuous} object with dimensions
#'   \code{segment × trial × channel}.
#' }
#'
#' \strong{Segmentation:}
#' Segment labels define the output time axis. They need not correspond to
#' physical seconds, although that is the default when \code{seg = NULL}. The
#' output \code{TimeTrace} is simply the sorted unique numeric segment labels.
#'
#' \strong{Output semantics:}
#' The returned object is called \code{EPhysContinuous} because it stores one
#' numeric value per segment, run, and channel. However, those values are not
#' samples of a continuous analog signal; they are summary statistics computed
#' from all spikes falling into each segment.
#'
#' \strong{Metric implementation:}
#' Exact peak finding, circular windowing, and tie handling are delegated to the
#' underlying helper functions
#' \code{\link{rayleigh_peak_by_segment}} and
#' \code{\link{phasemetrics_opponent_modulation_depth_by_segment}}.
#'
#' @return
#' An \code{\link[EPhysData:EPhysContinuous-class]{EPhysContinuous}} object with
#' \code{Data} shaped \code{[nSegments × nTrialsOut × nChannels]}. Each entry is
#' the selected phase metric for one segment, one output trial, and one channel.
#'
#' The returned object has:
#' \describe{
#'   \item{\code{TimeTrace}}{The numeric segment labels.}
#'   \item{\code{StimulusTrace}}{An empty numeric vector, because the output no
#'   longer represents a stimulus trace.}
#'   \item{\code{Metadata}}{The metadata of the processed object, after optional
#'   trial flattening.}
#'   \item{\code{Channels}}{The same channel vector as the processed input.}
#' }
#'
#' @section Low-level functions used on this page:
#' \describe{
#'   \item{\code{\link{phasemetrics_stimulus_preprocess}}}{Convert a stimulus
#'   trace into a phase trace.}
#'   \item{\code{\link{spike_phase_from_stimulus}}}{Assign stimulus phases to
#'   individual spike times.}
#'   \item{\code{\link{rayleigh_peak_by_segment}}}{Compute per-segment Rayleigh
#'   peak strength.}
#'   \item{\code{\link{phasemetrics_opponent_modulation_depth_by_segment}}}{Compute
#'   histogram-based phase modulation metrics by segment.}
#' }
#'
#' @examples
#' # Compute per-second phase metrics after flattening repeated trials:
#' # pm <- PhaseMetrics(ev, method = "rayleighZ")
#'
#' # Compute opponent modulation on custom segments:
#' # pm <- PhaseMetrics(
#' #   ev,
#' #   method = "opponent_best",
#' #   seg = function(t) floor(t / 0.5),
#' #   window_half_deg = 30
#' # )
#'
#' @importClassesFrom EPhysData EPhysEvents EPhysContinuous
#' @importFrom EPhysData Metadata Channels TimeTrace TimeUnits StimulusTrace lapply HasStimulus
#' @seealso
#'   \code{\link{phasemetrics_stimulus_preprocess}},
#'   \code{\link{spike_phase_from_stimulus}},
#'   \code{\link{rayleigh_peak_by_segment}},
#'   \code{\link{phasemetrics_opponent_modulation_depth_by_segment}},
#'   \code{EPhysData::\link[EPhysData:FlattenTrials]{FlattenTrials}},
#'   \code{EPhysData::\link[EPhysData:TimeTrace]{TimeTrace}},
#'   \code{EPhysData::\link[EPhysData:StimulusTrace]{StimulusTrace}}
#' @export
setGeneric("PhaseMetrics", function(X,
                                    method = c("rayleighZ","peak_trough","opponent180","opponent90","opponent_best","isi50"),
                                    seg = NULL,
                                    bin_width_deg = 15,
                                    window_half_deg = 45,
                                    flatten_trials = TRUE,
                                    parallel = !interactive(),
                                    progress = interactive(),
                                    lapply_error = c("stop","warn")[1],
                                    ...) {
  standardGeneric("PhaseMetrics")
})

#' @rdname PhaseMetrics
#' @export
setMethod("PhaseMetrics", signature(X = "EPhysEvents"),
          function(X,
                   method = c("rayleighZ","peak_trough","opponent180","opponent90","opponent_best","isi50"),
                   seg = NULL,
                   bin_width_deg = 15,
                   window_half_deg = 45,
                   flatten_trials = TRUE,
                   parallel = !interactive(),
                   progress = interactive(),
                   lapply_error = c("stop","warn")[1],
                   ...) {

            method <- match.arg(method)

            # --- optionally flatten trials by RecordingID
            x_use <- if (isTRUE(flatten_trials)) FlattenTrials(X) else X

            # --- basic checks
            if (length(TimeTrace(x_use)) == 0L)
              stop("TimeTrace is empty; required to derive phases/segments.")
            if (!(HasStimulus(x_use)))
              stop("StimulusTrace is empty; cannot derive phases. Provide a StimulusTrace.")

            # --- build phase trace from StimulusTrace
            phase_df <- phasemetrics_stimulus_preprocess(
              t    = to_num(TimeTrace(x_use)),
              stim = StimulusTrace(x_use),
              ...
            )
            if (!is.data.frame(phase_df) || !all(c("t","phase_deg") %in% names(phase_df)))
              stop("phasemetrics_stimulus_preprocess() must return a data.frame with columns 't' and 'phase_deg'.")

            # --- segment labels (define global segment grid)
            tt <- TimeTrace(x_use)
            seg_vec <- NULL; seg_fun <- NULL; tu <- TimeUnits(x_use)

            if (is.null(seg)) {
              seg_vec <- floor(tt)
            } else if (is.function(seg)) {
              seg_fun <- seg; seg_vec <- seg(tt); tu <- ""
            } else if (is.numeric(seg)) {
              if (length(seg) != length(tt))
                stop("If 'seg' is numeric, it must have the same length as TimeTrace.")
              if ("units" %in% class(seg)) tu <- deparse_unit(seg)
              seg_vec <- as.numeric(seg)
            } else {
              stop("'seg' must be NULL, a function, or a numeric vector aligned to TimeTrace.")
            }

            seg_levels <- sort(unique(as.numeric(seg_vec)))
            nSeg <- length(seg_levels)

            # --- dims & metadata
            md  <- Metadata(x_use)
            ch  <- Channels(x_use)
            nTr <- nrow(md)
            nCh <- length(ch)

            # --- helpers ------------------------------------------------------
            t_phase  <- phase_df$t
            ph_trace <- phase_df$phase_deg

            map_spikes_to_seg <- function(spk) {
              if (length(spk) == 0L) return(numeric(0))
              if (!is.null(seg_fun)) {
                as.numeric(seg_fun(spk))
              } else {
                approx(
                  x = tt,
                  y = as.numeric(seg_vec),
                  xout = spk,
                  method = "constant",
                  # <- step function, not linear
                  f      = 0,
                  rule = 2,
                  ties = "ordered"
                )$y
              }
            }

            compute_by_segment <- function(df) {
              if (method == "rayleighZ") {
                # must return a named numeric with names = segment labels
                v <- rayleigh_peak_by_segment(
                  df,
                  phase_col       = "p",
                  seg_col         = "seg",
                  bin_width_deg   = bin_width_deg,
                  window_half_deg = window_half_deg
                )
                # ensure names present and numeric
                if (is.null(names(v)) && length(v))
                  names(v) <- as.character(sort(unique(df$seg)))
                v
              } else {
                m <- phasemetrics_opponent_modulation_depth_by_segment(
                  df,
                  phase_col       = "p",
                  seg_col         = "seg",
                  bin_width_deg   = bin_width_deg,
                  window_half_deg = window_half_deg,
                  methods         = method
                )
                v <- as.numeric(m[, method, drop = TRUE])
                names(v) <- rownames(m)
                v
              }
            }

            # per-channel worker for lapply(EPhysEvents, ...)
            per_channel_fun <- function(spk) {
              if (is.null(spk) || !length(spk)) return(numeric(0))
              p_deg  <- spike_phase_from_stimulus(spk, t_phase, ph_trace)
              seg_spk <- map_spikes_to_seg(spk)
              df <- data.frame(p = p_deg, seg = seg_spk)
              compute_by_segment(df)  # named numeric by segment
            }

            results <- lapply(
              X        = x_use,
              FUN      = per_channel_fun,
              parallel = parallel,
              error    = match.arg(lapply_error, c("stop","warn")),
              progress = progress
            )
            # --- allocate output and pack results into [seg × trial × channel]
            trial_labels <- as.character(md$RunUID)
            names(results)<-NULL
            out <- array(NA_real_, dim = c(nSeg, nTr, nCh),
                         dimnames = list(
                           time    = as.character(seg_levels),
                           trial   = trial_labels,
                           channel = ch
                         ))

            for (i in seq_len(nTr)) {
              row_res <- results[[i]]
              if (!is.list(row_res)) {
                stop("PhaseMetrics: Aggregation into matrix failed: no valid data for run ",
                     md$RunUID[i], ".")
              }

              for (cj in seq_len(nCh)) {
                v <- row_res[[cj]]

                if (!is.numeric(v)) {
                  stop("PhaseMetrics: expected numeric segment vector for channel ", cj, ".")
                }

                vals <- rep(0, nSeg)
                if (length(v) > 0) {
                  lab  <- as.numeric(names(v))
                  idx  <- match(lab, seg_levels)  # position of each label in global grid
                  keep <- !is.na(idx)
                  if (any(!keep)) {
                    warning("PhaseMetrics: some segment labels for channel ", cj,
                            " (run ", md$RunUID[i], ") do not match seg_levels and are dropped.")
                  }
                  vals[idx[keep]] <- v[keep]
                }

                out[, i, cj] <- vals
              }
            }


            # --- build EPhysContinuous
            newEPhysContinuous(
              Data             = out,
              TimeTrace        = as.numeric(seg_levels),
              Channels         = ch,
              Channel_Metadata = x_use@Channel_Metadata,
              Metadata         = md,
              StimulusTrace    = numeric(0),
              TimeUnits        = tu,
              StimulusUnits    = x_use@StimulusUnits,
              ExamInfo         = x_use@ExamInfo,
              SubjectInfo      = x_use@SubjectInfo,
              Imported         = x_use@Imported
            )
          })
