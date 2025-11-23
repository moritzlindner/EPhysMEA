#' Convert EPhysEvents to per-segment phase metrics (continuous form)
#'
#' Computes a per-segment phase-based amplitude metric from an \code{EPhysEvents}
#' object, using phases derived from \code{StimulusTrace} via
#' \code{phasemetrics_stimulus_preprocess()}, and returns an \code{EPhysContinuous}
#' whose time axis is the set of segment labels.
#'
#' @inheritParams Bin
#' @param method Character; one of:
#'   \itemize{
#'     \item \code{"rayleighZ"} — Rayleigh strength (\eqn{Z = nR^2}) inside a peak-centered window.
#'     \item \code{"peak_trough"} — max−min of circular phase histogram counts.
#'     \item \code{"opponent180"} — windowed count at peak minus count at peak+180°.
#'     \item \code{"opponent90"}  — windowed count at peak minus the smaller of peak±90° windows.
#'     \item \code{"opponent_best"} — \code{max(opponent180, opponent90)}.
#'     \item \code{"isi50"} — 180° minus the shortest circular arc covering half the spikes.
#'   }
#' @param seg Either \code{NULL} (default; uses \code{floor(TimeTrace)}),
#'   a \emph{function} mapping spike times to segment labels, or a \emph{numeric}
#'   vector the same length as \code{TimeTrace} giving a label per time sample
#'   (spike labels obtained by interpolation).
#' @param bin_width_deg Histogram bin width (default \code{15}°).
#' @param window_half_deg Half-width of the peak/opponent windows (default \code{45}°).
#' @param flatten_trials Logical; if \code{TRUE} (default) first call
#'   \code{FlattenTrials(X)} so each output trial corresponds to one \code{RecordingID}.
#'   If \code{FALSE}, metrics are computed per input trial.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @param ... Passed to \code{phasemetrics_stimulus_preprocess()} (e.g., \code{upsample}, \code{level}, \code{direction}).
#'
#' @return An \code{EPhysContinuous} with \code{Data} shaped
#'   \code{[nSegments × nTrialsOut × nChannels]} containing the chosen metric per segment.
#'   \code{TimeTrace} equals the numeric segment labels; \code{StimulusTrace} is length 0.
#' @keywords internal
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels TimeTrace StimulusTrace lapply HasStimulus
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
