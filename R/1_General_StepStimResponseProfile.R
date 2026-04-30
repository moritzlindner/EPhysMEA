#' Step-stimulus response profile metrics for `EPhysEvents`
#'
#' \strong{StepStimResponseProfile:} Compute response-profile metrics for step responses in an
#' \code{\link[EPhysData:EPhysEvents-class]{EPhysEvents}} object and optionally
#' summarize them across repeated runs that share the same \code{RecordingID}.
#'
#' The method proceeds in four stages:
#' \enumerate{
#'   \item For each run × channel, compute baseline-subtracted ON/OFF AUC
#'   components together with derived run-level \code{PolarityScore},
#'   \code{PolarityIndex}, and \code{BiasIndex}.
#'   \item Summarize these run-level polarity and bias metrics per
#'   \code{RecordingID × Channel}. If multiple runs are available, repeatability
#'   statistics are computed across runs and used to derive the unit-level
#'   calls. If only one run is available, repeatability testing is skipped and
#'   the observed run-level polarity and bias calls are used directly.
#'   \item Compute downstream per-run metrics whose definition depends on the
#'   summarized polarity and bias results.
#'   \item If \code{average_by_recordingID = TRUE}, summarize these downstream
#'   metrics per \code{RecordingID × Channel}.
#' }
#'
#' Polarity and bias summarize different aspects of the response:
#' \describe{
#'   \item{\code{PolarityScore} / \code{PolarityIndex}}{Whether positive
#'   modulation or ON-step suppression dominates. Positive values indicate that
#'   positive modulation dominates; negative values indicate dominant
#'   ON-step suppression.}
#'   \item{\code{BiasIndex}}{Whether the response modulation is predominantly
#'   associated with stimulus onset or stimulus offset, irrespective of polarity.
#'   Positive values indicate onset-dominated modulation; negative values
#'   indicate offset-dominated modulation.}
#' }
#'
#' The downstream metrics currently comprise:
#' \describe{
#'   \item{\code{Shinomoto_Baseline_LV}}{Shinomoto local variation (LV) of the
#'   spike train in the baseline window \code{[step_onset - baseline_win,
#'   step_onset)}. The LV formula follows Shinomoto et al.; its application to
#'   the pre-step baseline window is a derived use in this method.}
#'   \item{\code{Response_Amplitude}}{Non-negative response amplitude in spikes/s,
#'   chosen from the ON, OFF, or suppressive ON component according to the
#'   summarized polarity and bias calls.}
#'   \item{\code{ISI_Latency}}{ISI-based response latency measured relative to
#'   step onset or step offset according to the summarized polarity and bias
#'   calls. Negative-polarity units are treated as suppressive and use the
#'   latency-to-suppression rule.}
#'   \item{\code{Farrow_TransienceIndex}}{Derived peak-normalized transience
#'   index computed on the response component selected by the summarized
#'   polarity and bias calls. For positive responses, the positive
#'   baseline-subtracted component is used; for negative responses, the
#'   suppressive component is rectified and used analogously.}
#'   \item{\code{Liang_OnOff_PreferenceIndex}}{Liang-style ON/OFF preference
#'   index computed from baseline-subtracted ON- and OFF-aligned response traces
#'   of length \code{step_offset - step_onset}. As in Liang et al., positive
#'   values are set to zero for suppressive responses before the index is
#'   computed. In this implementation, suppressive cases are identified from the
#'   summarized \code{ProfilePolarity}, and an additional post-hoc sign
#'   inversion is applied so that ON-step suppression yields a negative,
#'   OFF-labeled value.}
#' }
#'
#' @param X An \code{EPhysEvents} object.
#' @param step_range Numeric or \code{units} length-2 vector \code{c(start, end)}
#'   defining the step interval. It is interpreted as \code{[start, end)} and
#'   must satisfy \code{start < end}. If given as \code{units}, it is converted
#'   to \code{\link[EPhysData:TimeUnits]{TimeUnits}(X)}.
#' @param average_by_recordingID Logical; if \code{TRUE}, average all numeric
#'   metric columns across rows sharing the same \code{RecordingID} and
#'   \code{Channel} using \code{mean(na.rm = TRUE)}. In addition,
#'   \code{Suppressed_By_Contrast} is summarized as the proportion of runs
#'   classified as suppressive and returned as
#'   \code{Suppressed_By_Contrast_Prop}. Default: \code{FALSE}.
#' @param decision_alpha Significance threshold used to convert the run-level
#'   repeatability tests for \code{PolarityIndex} and \code{BiasIndex} into the
#'   unit-level calls \code{ProfilePolarity} and \code{ProfileBias}. If only one
#'   run is available for a \code{RecordingID × Channel} group, repeatability
#'   testing is skipped and this threshold is not used for that group. Default:
#'   \code{0.05}.
#' @param latency_alpha Quantile used to define the baseline ISI threshold for
#'   the latency detector. Default: \code{0.05}.
#' @param p_adjust Multiple-testing adjustment method passed to
#'   \code{\link[stats]{p.adjust}} for the group-level polarity and bias
#'   p-values. If all processed \code{RecordingID × Channel} groups have only
#'   one run, no repeatability p-values are generated and this argument is not
#'   used. Default: \code{"none"}.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @inheritDotParams step_stim_response_auc_components
#'   on_win off_win baseline_win binwidth smooth gauss_sd
#' @param ... Further arguments forwarded to
#'   \code{\link{step_stim_response_auc_components}}.
#'
#' @return
#' If \code{average_by_recordingID = FALSE} (default), a wide \code{data.frame}
#' with one row per run × channel and identifier columns plus the following
#' metrics:
#' \code{Baseline_Rate}, \code{Dominant_Component}, \code{PolarityIndex},
#' \code{PolarityScore}, \code{BiasIndex}, \code{Shinomoto_Baseline_LV},
#' \code{Response_Amplitude}, \code{ISI_Latency},
#' \code{Farrow_TransienceIndex}, and
#' \code{Liang_OnOff_PreferenceIndex}.
#'
#' If \code{average_by_recordingID = TRUE}, a wide \code{data.frame} with one row
#' per \code{RecordingID × Channel}, containing the same metrics summarized
#' across runs. Numeric metrics are averaged with \code{mean(na.rm = TRUE)}.
#' \code{Dominant_Component} is summarized as the unique most frequent
#' non-missing category across runs and returned as \code{NA} if no unique mode
#' exists.
#'
#' @section Low-level functions on this page:
#' \describe{
#'   \item{\code{step_stim_response_auc_components(x, ...)}}{Compute the
#'   baseline-subtracted ON/OFF AUC components for a single spike train.}
#'   \item{\code{step_stim_response_polarity(auc_metrics)}}{Derive
#'   \code{PolarityScore} and \code{PolarityIndex} from the AUC components.}
#'   \item{\code{step_stim_response_bias(auc_metrics)}}{Derive
#'   \code{BiasIndex} from the AUC components.}
#'   \item{\code{step_stim_response_polarity_bias(x, ...)}}{Convenience wrapper
#'   that combines the per-run AUC-component, polarity, and bias calculations.}
#'   \item{\code{step_stim_response_baseline_lv(x, ...)}}{Compute the baseline
#'   Shinomoto LV.}
#'   \item{\code{step_stim_response_amplitude(...)}}{Compute the downstream
#'   response amplitude from the summarized polarity and bias calls.}
#'   \item{\code{step_stim_response_isi_latency(x, ...)}}{Compute the downstream
#'   ISI-based response latency from the summarized polarity and bias calls.}
#'   \item{\code{step_stim_response_transience(x, ...)}}{Compute the downstream
#'   transience index from the summarized polarity and bias calls.}
#'   \item{\code{step_stim_response_liang_index(x, ...)}}{Compute the Liang-style
#'   ON/OFF preference index from the summarized polarity call.}
#'   \item{\code{step_stim_response_dependent_metrics(x, ...)}}{Convenience
#'   wrapper that combines the downstream per-run calculations.}
#' }
#'
#' #' @references
#' Farrow K, Masland RH.
#' Physiological clustering of visual channels in the mouse retina.
#' \emph{J Neurophysiol}. 2011;105(4):1516-1530.
#'
#' Shinomoto S, Shima K, Tanji J.
#' Differences in spiking patterns among cortical neurons.
#' \emph{Biosystems}. 2005;79(1-3):57-65.
#'
#' Liang L, Fratzl A, Goldey GJ, Ramesh RN, Sugden AU, Morgan JL, Chen C,
#' Andermann ML. A Fine-Scale Functional Logic to Convergence from Retina to
#' Thalamus. \emph{Cell}. 2018;173(6):1343-1355.e24.
#' doi:10.1016/j.cell.2018.04.041
#'
#' @importClassesFrom EPhysData EPhysEvents
#' @importFrom EPhysData Metadata Channels TimeUnits lapply
#' @importFrom units set_units
#' @importFrom stats binom.test p.adjust mad median quantile
#' @export
setGeneric(
  "StepStimResponseProfile",
  function(X,
           step_range,
           average_by_recordingID = FALSE,
           decision_alpha = 0.05,
           latency_alpha = 0.05,
           p_adjust = "none",
           parallel = !interactive(),
           ...)
    standardGeneric("StepStimResponseProfile")
)

#' @rdname StepStimResponseProfile
#' @export
setMethod(
  "StepStimResponseProfile",
  "EPhysEvents",
  function(X,
           step_range,
           average_by_recordingID = FALSE,
           decision_alpha = 0.05,
           latency_alpha = 0.05,
           p_adjust = "none",
           parallel = !interactive(),
           ...) {

    p_adjust <- match.arg(p_adjust, p.adjust.methods)
    p_value_args_nondefault <- !identical(decision_alpha, 0.05) ||
      !identical(p_adjust, "none")

    md <- Metadata(X)

    if (length(unique(md$Experiment)) != 1L) {
      stop(
        "Only EPhysEvents containing one single type of Experiment as defined ",
        "in the respective metadata column can safely be processed by this method."
      )
    }

    dots <- list(...)
    on_win <- if ("on_win" %in% names(dots)) dots$on_win else 0.4
    off_win <- if ("off_win" %in% names(dots)) dots$off_win else on_win
    baseline_win <- if ("baseline_win" %in% names(dots)) dots$baseline_win else 2
    binwidth <- if ("binwidth" %in% names(dots)) dots$binwidth else 0.02
    smooth <- if ("smooth" %in% names(dots)) dots$smooth else FALSE
    gauss_sd <- if ("gauss_sd" %in% names(dots)) dots$gauss_sd else 0.02

    on_win <- to_num(on_win)
    off_win <- to_num(off_win)
    baseline_win <- to_num(baseline_win)
    binwidth <- to_num(binwidth)
    gauss_sd <- to_num(gauss_sd)

    if ("units" %in% class(step_range)) {
      tryCatch({
        step_range <- set_units(step_range, TimeUnits(X), mode = "standard")
      }, error = function(e) {
        stop(
          "Time unit conversion for step_range failed with message: ",
          e,
          " Is provided unit compatible with 'TimeUnits(X)'?"
        )
      })
    }

    step_rangen <- suppressWarnings(as.numeric(step_range))
    if (!is.numeric(step_rangen) ||
        length(step_rangen) != 2L ||
        any(!is.finite(step_rangen))) {
      stop("step_range must be a finite numeric (or units) length-2 vector.")
    }

    step_rangen <- sort(step_rangen)
    if (!(step_rangen[2] > step_rangen[1])) {
      stop("step_range must satisfy start < end.")
    }

    step_onset  <- step_rangen[1]
    step_offset <- step_rangen[2]

    ## Stage 1: per-run AUC, polarity, and bias
    metrics_list <- EPhysData::lapply(
      X,
      step_stim_response_polarity_bias,
      step_onset  = step_onset,
      step_offset = step_offset,
      ...,
      parallel = parallel
    )
    metrics_df <- nested2df(X, metrics_list)

    keys <- c("RecordingID", "Channel")

    wanted_metrics <- c(
      "Baseline_Rate",
      "Dominant_Component",
      "PolarityIndex",
      "Polarity_SignTest_p",
      "PolarityScore",
      "BiasIndex",
      "Bias_SignTest_p",
      "Shinomoto_Baseline_LV",
      "Response_Amplitude",
      "ISI_Latency",
      "Farrow_TransienceIndex",
      "Liang_OnOff_PreferenceIndex"
    )

    id_cols <- c("RunUID", "RecordingID", "Channel")

    ## Stage 2: summarize per-run polarity/bias and compute repeatability p values
    mean_or_na <- function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    }

    mode_or_na <- function(x) {
      x <- x[!is.na(x) & nzchar(as.character(x))]
      if (!length(x)) return(NA_character_)
      tab <- table(x)
      winners <- names(tab)[tab == max(tab)]
      if (length(winners) == 1L) winners else NA_character_
    }

    num_cols <- setdiff(
      names(metrics_df)[vapply(metrics_df, is.numeric, logical(1))],
      "RecordingID"
    )

    out_num <- aggregate(
      metrics_df[num_cols],
      by = metrics_df[keys],
      FUN = mean_or_na
    )

    split_idx <- interaction(metrics_df$RecordingID, metrics_df$Channel, drop = TRUE)

    out_rep <- do.call(
      rbind,
      lapply(split(metrics_df, split_idx), function(df) {

        ## Polarity repeatability from PolarityIndex
        p <- df$PolarityIndex
        p <- p[is.finite(p)]
        p_nz <- p[p != 0]

        n_runs <- nrow(df)
        n_nonzero <- length(p_nz)
        n_pos <- sum(p_nz > 0)
        n_neg <- sum(p_nz < 0)

        Polarity_Call <- if (!n_nonzero) {
          NA_character_
        } else if (n_pos > n_neg) {
          "positive"
        } else if (n_neg > n_pos) {
          "negative"
        } else {
          "ambiguous"
        }
        single_run <- n_runs <= 1L

        Polarity_SignTest_p <- if (!single_run && n_nonzero > 0L) {
          binom.test(
            x = max(n_pos, n_neg),
            n = n_nonzero,
            p = 0.5,
            alternative = "two.sided"
          )$p.value
        } else {
          NA_real_
        }

        ## Bias repeatability from BiasIndex
        b <- df$BiasIndex
        b <- b[is.finite(b)]
        b_nz <- b[b != 0]

        n_bias <- length(b)
        n_bias_nonzero <- length(b_nz)
        n_onset <- sum(b_nz > 0)
        n_offset <- sum(b_nz < 0)

        Bias_Call <- if (!n_bias_nonzero) {
          NA_character_
        } else if (n_onset > n_offset) {
          "onset"
        } else if (n_offset > n_onset) {
          "offset"
        } else {
          "balanced"
        }

        Bias_SignTest_p <- if (!single_run && n_bias_nonzero > 0L) {
          binom.test(
            x = max(n_onset, n_offset),
            n = n_bias_nonzero,
            p = 0.5,
            alternative = "two.sided"
          )$p.value
        } else {
          NA_real_
        }

        data.frame(
          RecordingID = df$RecordingID[1],
          Channel = df$Channel[1],
          RunUID = if (single_run) as.character(df$RunUID[1]) else "aggregated",
          N_Runs = n_runs,
          Polarity_Call = Polarity_Call,
          N_Runs_PolarityNonzero = n_nonzero,
          N_Runs_PolarityPositive = n_pos,
          N_Runs_PolarityNegative = n_neg,
          PolarityIndex_Median = if (length(p)) median(p, na.rm = TRUE) else NA_real_,
          PolarityIndex_MAD = if (length(p)) {
            mad(p, center = median(p, na.rm = TRUE), constant = 1, na.rm = TRUE)
          } else {
            NA_real_
          },
          Polarity_SignTest_p = Polarity_SignTest_p,

          Bias_Call = Bias_Call,
          N_Runs_BiasFinite = n_bias,
          N_Runs_BiasNonzero = n_bias_nonzero,
          N_Runs_OnsetDominated = n_onset,
          N_Runs_OffsetDominated = n_offset,
          BiasIndex_Median = if (length(b)) median(b, na.rm = TRUE) else NA_real_,
          BiasIndex_MAD = if (length(b)) {
            mad(b, center = median(b, na.rm = TRUE), constant = 1, na.rm = TRUE)
          } else {
            NA_real_
          },
          Bias_SignTest_p = Bias_SignTest_p,

          stringsAsFactors = FALSE
        )
      })
    )

    has_single_run_groups <- any(out_rep$N_Runs <= 1L, na.rm = TRUE)
    has_multi_run_groups <- any(out_rep$N_Runs > 1L, na.rm = TRUE)

    if (has_single_run_groups && !has_multi_run_groups && p_value_args_nondefault) {
      warning(
        "StepStimResponseProfile(): all RecordingID × Channel groups have only one run; ",
        "repeatability p-value testing is skipped, so non-default ",
        "'decision_alpha' and/or 'p_adjust' values are ignored.",
        call. = FALSE
      )
    } else if (has_single_run_groups && p_value_args_nondefault) {
      warning(
        "StepStimResponseProfile(): some RecordingID × Channel groups have only one run; ",
        "repeatability p-value testing is skipped for those groups, so non-default ",
        "'decision_alpha' and/or 'p_adjust' values are ignored for them.",
        call. = FALSE
      )
    }

    if (p_adjust != "none" && has_multi_run_groups) {
      if ("Polarity_SignTest_p" %in% names(out_rep)) {
        ok <- !is.na(out_rep$Polarity_SignTest_p)
        out_rep$Polarity_SignTest_p_adj <- NA_real_
        out_rep$Polarity_SignTest_p_adj[ok] <- p.adjust(
          out_rep$Polarity_SignTest_p[ok],
          method = p_adjust
        )
      }
      if ("Bias_SignTest_p" %in% names(out_rep)) {
        ok <- !is.na(out_rep$Bias_SignTest_p)
        out_rep$Bias_SignTest_p_adj <- NA_real_
        out_rep$Bias_SignTest_p_adj[ok] <- p.adjust(
          out_rep$Bias_SignTest_p[ok],
          method = p_adjust
        )
      }
    }

    out_rep$ProfilePolarity <- ifelse(
      out_rep$N_Runs <= 1L &
        out_rep$Polarity_Call %in% c("positive", "negative"),
      out_rep$Polarity_Call,
      ifelse(
        !is.na(out_rep$Polarity_SignTest_p) &
          out_rep$Polarity_SignTest_p < decision_alpha &
          out_rep$Polarity_Call %in% c("positive", "negative"),
        out_rep$Polarity_Call,
        "ambiguous"
      )
    )

    out_rep$ProfileBias <- ifelse(
      out_rep$N_Runs <= 1L &
        out_rep$Bias_Call %in% c("onset", "offset"),
      out_rep$Bias_Call,
      ifelse(
        !is.na(out_rep$Bias_SignTest_p) &
          out_rep$Bias_SignTest_p < decision_alpha &
          out_rep$Bias_Call %in% c("onset", "offset"),
        out_rep$Bias_Call,
        "ambiguous"
      )
    )

    ## Stage 3: per-run downstream metrics conditional on summarized polarity/bias
    raw_nested <- EPhysData::lapply(
      X,
      function(z) z,
      parallel = parallel
    )
    x_flat <- unlist(raw_nested, recursive = FALSE, use.names = FALSE)

    if (length(x_flat) != nrow(metrics_df)) {
      stop(
        "Internal error: flattened spike-train list length does not match ",
        "the number of per-run rows."
      )
    }

    key_run <- paste(metrics_df$RecordingID, metrics_df$Channel, sep = "\r")
    key_rep <- paste(out_rep$RecordingID, out_rep$Channel, sep = "\r")
    rep_match <- out_rep[match(key_run, key_rep), , drop = FALSE]

    dep_list <- Map(
      function(x, profile_polarity, profile_bias, baseline_rate,
               on_auc_pos, on_auc_neg, off_auc_pos) {
        step_stim_response_dependent_metrics(
          x = x,
          step_onset = step_onset,
          step_offset = step_offset,
          baseline_win = baseline_win,
          on_win = on_win,
          off_win = off_win,
          binwidth = binwidth,
          smooth = smooth,
          gauss_sd = gauss_sd,
          baseline_rate = baseline_rate,
          on_auc_pos = on_auc_pos,
          on_auc_neg = on_auc_neg,
          off_auc_pos = off_auc_pos,
          profile_polarity = profile_polarity,
          profile_bias = profile_bias,
          latency_alpha = latency_alpha
        )
      },
      x_flat,
      rep_match$ProfilePolarity,
      rep_match$ProfileBias,
      metrics_df$Baseline_Rate,
      metrics_df$On_AUC_Pos,
      metrics_df$On_AUC_Neg,
      metrics_df$Off_AUC_Pos
    )

    dep_df <- do.call(
      rbind,
      lapply(dep_list, function(z) as.data.frame(z, stringsAsFactors = FALSE))
    )
    rownames(dep_df) <- NULL

    rep_small <- rep_match[, c(
      "Polarity_Call", "Polarity_SignTest_p",
      "Bias_Call", "Bias_SignTest_p",
      "ProfilePolarity", "ProfileBias"
    ), drop = FALSE]

    per_run_out <- cbind(metrics_df, dep_df)

    if (isFALSE(average_by_recordingID)) {
      keep <- intersect(c(id_cols, wanted_metrics), names(per_run_out))
      return(per_run_out[, keep, drop = FALSE])
    }

    ## Stage 4: summarize downstream per-run metrics and restrict output columns
    num_keep <- intersect(
      c(
        "Baseline_Rate",
        "PolarityIndex",
        "PolarityScore",
        "Polarity_SignTest_p",
        "BiasIndex",
        "Bias_SignTest_p",
        "Shinomoto_Baseline_LV",
        "Response_Amplitude",
        "ISI_Latency",
        "Farrow_TransienceIndex",
        "Liang_OnOff_PreferenceIndex"
      ),
      names(per_run_out)
    )

    out_num_final <- aggregate(
      per_run_out[num_keep],
      by = per_run_out[keys],
      FUN = mean_or_na
    )

    out_dom <- aggregate(
      per_run_out["Dominant_Component"],
      by = per_run_out[keys],
      FUN = mode_or_na
    )

    key_num <- paste(out_num_final$RecordingID, out_num_final$Channel, sep = "\r")
    key_dom <- paste(out_dom$RecordingID, out_dom$Channel, sep = "\r")
    out_dom <- out_dom[match(key_num, key_dom), , drop = FALSE]

    out <- cbind(
      out_num_final,
      out_dom["Dominant_Component"]
    )

    key_out <- paste(out$RecordingID, out$Channel, sep = "\r")
    key_rep <- paste(out_rep$RecordingID, out_rep$Channel, sep = "\r")
    out$RunUID <- out_rep$RunUID[match(key_out, key_rep)]

    keep <- intersect(c(id_cols, wanted_metrics), names(out))
    out[, keep, drop = FALSE]
  }
)

#' Low-level ON/OFF AUC components for a single step response
#'
#' \strong{step_stim_response_auc_components:} Low-level function to compute baseline-subtracted ON/OFF AUC components for a single spike train.
#'
#' @param x Numeric vector of spike or event times for a single spike train.
#' @param step_onset Numeric or \code{units} scalar giving step onset time.
#' @param step_offset Numeric or \code{units} scalar giving step offset time.
#' @param on_win Length of the ON analysis window. Default: \code{0.4}.
#' @param off_win Length of the OFF/post-step analysis window.
#'   Default: \code{on_win}.
#' @param baseline_win Length of the pre-step baseline window immediately
#'   preceding \code{step_onset}. Default: \code{2}.
#' @param binwidth Bin width used for the rate traces and AUC calculations.
#'   Default: \code{0.02}.
#' @param smooth Logical; if \code{TRUE}, smooth the binned traces with a
#'   Gaussian kernel before AUC calculation. Default: \code{FALSE}.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'   Default: \code{0.02}.
#'
#' @return For step_stim_response_auc_components: A named list containing:
#' \describe{
#'   \item{\code{Baseline_Rate}}{Mean baseline spike rate in spikes/s.}
#'   \item{\code{Mean_Step_Rate}}{Mean spike rate over the full step interval
#'     \code{[step_onset, step_offset)}.}
#'   \item{\code{On_AUC_Signed}}{Signed baseline-subtracted AUC in the ON window.}
#'   \item{\code{Off_AUC_Signed}}{Signed baseline-subtracted AUC in the OFF window.}
#'   \item{\code{On_AUC_Pos}}{Positive baseline-subtracted AUC in the ON window.}
#'   \item{\code{On_AUC_Neg}}{Negative baseline-subtracted AUC in the ON window.}
#'   \item{\code{Off_AUC_Pos}}{Positive baseline-subtracted AUC in the OFF window.}
#'   \item{\code{Off_AUC_Neg}}{Negative baseline-subtracted AUC in the OFF window.}
#'   \item{\code{Positive_AUC}}{\code{max(On_AUC_Pos, Off_AUC_Pos)}.}
#'   \item{\code{Negative_AUC}}{\code{On_AUC_Neg}.}
#'   \item{\code{Dominant_Component}}{Component with the largest contribution
#'     among \code{On_Pos}, \code{Off_Pos}, and \code{On_Neg}.}
#' }
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_auc_components <- function(x,
                                              step_onset,
                                              step_offset,
                                              on_win = 0.4,
                                              off_win = on_win,
                                              baseline_win = 2,
                                              binwidth = 0.02,
                                              smooth = FALSE,
                                              gauss_sd = 0.02) {

  x <- sort(to_num(x)[is.finite(to_num(x))])

  step_onset   <- to_num(step_onset)
  step_offset  <- to_num(step_offset)
  on_win       <- to_num(on_win)
  off_win      <- to_num(off_win)
  baseline_win <- to_num(baseline_win)
  binwidth     <- to_num(binwidth)
  gauss_sd     <- to_num(gauss_sd)

  stopifnot(
    is.numeric(x),
    is.numeric(step_onset),
    is.numeric(step_offset),
    is.finite(step_offset - step_onset),
    (step_offset - step_onset) > 0,
    is.finite(on_win), on_win > 0,
    is.finite(off_win), off_win > 0,
    is.finite(baseline_win), baseline_win > 0,
    is.finite(binwidth), binwidth > 0
  )

  step_win <- step_offset - step_onset

  x <- x[x >= (step_onset - baseline_win) & x < (step_offset + off_win)]

  Baseline_Rate <- count_in(x, step_onset - baseline_win, step_onset) / baseline_win

  Mean_Step_Rate <- window_mean_rate(
    spk = x,
    t0 = step_onset,
    win = step_win,
    binwidth = binwidth,
    gauss_sd = gauss_sd
  )

  on_obj <- make_rate_trace(
    spk = x,
    t0 = step_onset,
    win = on_win,
    binwidth = binwidth,
    smooth = smooth,
    gauss_sd = gauss_sd
  )
  off_obj <- make_rate_trace(
    spk = x,
    t0 = step_offset,
    win = off_win,
    binwidth = binwidth,
    smooth = smooth,
    gauss_sd = gauss_sd
  )

  on_dt <- diff(on_obj$breaks)
  off_dt <- diff(off_obj$breaks)

  on_trace  <- on_obj$rate  - Baseline_Rate
  off_trace <- off_obj$rate - Baseline_Rate

  On_AUC_Signed <- sum(on_trace * on_dt, na.rm = TRUE)
  Off_AUC_Signed <- sum(off_trace * off_dt, na.rm = TRUE)

  On_AUC_Pos <- sum(pmax(on_trace, 0) * on_dt, na.rm = TRUE)
  On_AUC_Neg <- sum(-pmin(on_trace, 0) * on_dt, na.rm = TRUE)

  Off_AUC_Pos <- sum(pmax(off_trace, 0) * off_dt, na.rm = TRUE)
  Off_AUC_Neg <- sum(-pmin(off_trace, 0) * off_dt, na.rm = TRUE)

  Positive_AUC <- max(On_AUC_Pos, Off_AUC_Pos)
  Negative_AUC <- On_AUC_Neg

  dom_vals <- c(On_Pos = On_AUC_Pos, Off_Pos = Off_AUC_Pos, On_Neg = On_AUC_Neg)
  Dominant_Component <- if (all(!is.finite(dom_vals)) || max(dom_vals, na.rm = TRUE) <= 0) {
    NA_character_
  } else {
    names(dom_vals)[which.max(dom_vals)]
  }

  list(
    Baseline_Rate = as.numeric(Baseline_Rate),
    Mean_Step_Rate = as.numeric(Mean_Step_Rate),
    On_AUC_Signed = as.numeric(On_AUC_Signed),
    Off_AUC_Signed = as.numeric(Off_AUC_Signed),
    On_AUC_Pos = as.numeric(On_AUC_Pos),
    On_AUC_Neg = as.numeric(On_AUC_Neg),
    Off_AUC_Pos = as.numeric(Off_AUC_Pos),
    Off_AUC_Neg = as.numeric(Off_AUC_Neg),
    Positive_AUC = as.numeric(Positive_AUC),
    Negative_AUC = as.numeric(Negative_AUC),
    Dominant_Component = as.character(Dominant_Component)
  )
}

#' Derive polarity metrics from step-response AUC components
#'
#' \strong{step_stim_response_polarity:} Low-level function to compute polarity metrics from the output of
#' \code{\link{step_stim_response_auc_components}}.
#'
#' The polarity score is defined as
#' \deqn{
#' \mathrm{PolarityScore} =
#' \mathrm{Positive\_AUC} - \mathrm{Negative\_AUC},
#' }
#' where \code{Positive_AUC = max(On_AUC_Pos, Off_AUC_Pos)} and
#' \code{Negative_AUC = On_AUC_Neg}. The normalized polarity index is
#' \deqn{
#' \mathrm{PolarityIndex} =
#' \frac{\mathrm{Positive\_AUC} - \mathrm{Negative\_AUC}}
#'      {\mathrm{Positive\_AUC} + \mathrm{Negative\_AUC}}.
#' }
#'
#' @param auc_metrics Named list produced by
#'   \code{\link{step_stim_response_auc_components}}.
#'
#' @return For step_stim_response_polarity: A named list with elements \code{PolarityScore} and
#'   \code{PolarityIndex}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_polarity <- function(auc_metrics) {
  eps <- sqrt(.Machine$double.eps)

  if (is.null(auc_metrics$Positive_AUC) || is.null(auc_metrics$Negative_AUC)) {
    stop(
      "'auc_metrics' must be a single output list from ",
      "'step_stim_response_auc_components()' containing ",
      "'Positive_AUC' and 'Negative_AUC'."
    )
  }

  Positive_AUC <- as.numeric(auc_metrics$Positive_AUC)
  Negative_AUC <- as.numeric(auc_metrics$Negative_AUC)

  PolarityScore <- Positive_AUC - Negative_AUC
  PolarityIndex <- if ((Positive_AUC + Negative_AUC) > eps) {
    (Positive_AUC - Negative_AUC) / (Positive_AUC + Negative_AUC)
  } else {
    0
  }

  list(
    PolarityScore = as.numeric(PolarityScore),
    PolarityIndex = as.numeric(PolarityIndex)
  )
}

#' Derive onset-versus-offset bias from step-response AUC components
#'
#' \strong{step_stim_response_bias:} Low-level function to compute onset-versus-offset bias from the output of
#' \code{\link{step_stim_response_auc_components}}.
#'
#' The bias index is defined as
#' \deqn{
#' \mathrm{BiasIndex} =
#' \frac{(\mathrm{On\_AUC\_Pos} + \mathrm{On\_AUC\_Neg}) -
#'       (\mathrm{Off\_AUC\_Pos} + \mathrm{Off\_AUC\_Neg})}
#'      {(\mathrm{On\_AUC\_Pos} + \mathrm{On\_AUC\_Neg}) +
#'       (\mathrm{Off\_AUC\_Pos} + \mathrm{Off\_AUC\_Neg})}.
#' }
#' Positive values indicate onset-dominated modulation; negative values indicate
#' offset-dominated modulation.
#'
#' @param auc_metrics Named list produced by
#'   \code{\link{step_stim_response_auc_components}}.
#'
#' @return For step_stim_response_bias: A named list with element \code{BiasIndex}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_bias <- function(auc_metrics) {
  eps <- sqrt(.Machine$double.eps)

  On_Modulation_AUC  <- as.numeric(auc_metrics$On_AUC_Pos)  + as.numeric(auc_metrics$On_AUC_Neg)
  Off_Modulation_AUC <- as.numeric(auc_metrics$Off_AUC_Pos) + as.numeric(auc_metrics$Off_AUC_Neg)

  den <- On_Modulation_AUC + Off_Modulation_AUC

  BiasIndex <- if (den > eps) {
    (On_Modulation_AUC - Off_Modulation_AUC) / den
  } else {
    0
  }

  list(
    BiasIndex = as.numeric(BiasIndex)
  )
}

#' Shinomoto LV
#'
#' \strong{step_stim_response_baseline_lv:} Low-level function to compute
#' Shinomoto local variation (LV) in the baseline window before step onset
#'
#' @param x Numeric vector of spike times for a single spike train.
#' @param step_onset Numeric scalar giving step onset time.
#' @param baseline_win Length of the baseline window immediately preceding
#'   \code{step_onset}.
#'
#' @return For \code{step_stim_response_baseline_lv()}, a named list with element
#'   \code{Shinomoto_Baseline_LV}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_baseline_lv <- function(x,
                                           step_onset,
                                           baseline_win = 2) {
  x <- sort(to_num(x)[is.finite(to_num(x))])
  step_onset <- to_num(step_onset)
  baseline_win <- to_num(baseline_win)

  base_times <- x[x >= (step_onset - baseline_win) & x < step_onset]

  if (length(base_times) == 0L) {
    lv <- 1
  } else {
    lv <- spike_irreg_lv_shinomoto(base_times)
  }

  if (length(lv) == 0L || is.na(lv)) {
    lv <- 1
  }

  list(
    Shinomoto_Baseline_LV = as.numeric(lv)
  )
}

#' Response Amplitude
#'
#' \strong{step_stim_response_amplitude:} low-level function to compute response amplitude from summarized polarity and bias calls
#'
#' @param on_auc_pos Positive ON AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param on_auc_neg Negative ON AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param off_auc_pos Positive OFF AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param on_win Length of the ON analysis window.
#' @param off_win Length of the OFF analysis window.
#' @param profile_polarity Unit-level summarized polarity call.
#' @param profile_bias Unit-level summarized bias call.
#'
#' @return For \code{step_stim_response_amplitude()}, a named list with element
#'   \code{Response_Amplitude}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_amplitude <- function(on_auc_pos,
                                         on_auc_neg,
                                         off_auc_pos,
                                         on_win,
                                         off_win,
                                         profile_polarity = c("positive", "negative", "ambiguous"),
                                         profile_bias = c("onset", "offset", "ambiguous")) {

  on_auc_pos <- to_num(on_auc_pos)
  on_auc_neg <- to_num(on_auc_neg)
  off_auc_pos <- to_num(off_auc_pos)
  on_win <- to_num(on_win)
  off_win <- to_num(off_win)

  on_amp <- on_auc_pos / on_win
  off_amp <- off_auc_pos / off_win
  neg_amp <- on_auc_neg / on_win

  amp <- if (identical(profile_polarity, "negative")) {
    neg_amp
  } else if (identical(profile_polarity, "positive")) {
    if (identical(profile_bias, "onset")) {
      on_amp
    } else if (identical(profile_bias, "offset")) {
      off_amp
    } else {
      vals <- c(on_amp, off_amp)
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }
  } else {
    NA_real_
  }

  list(
    Response_Amplitude = as.numeric(amp)
  )
}

#' ISI-based response latency
#'
#' \strong{step_stim_response_isi_latency:} low-level function to compute ISI-based response latency from summarized polarity and bias calls
#'
#' @param x Numeric vector of spike times for a single spike train.
#' @param step_onset Numeric scalar giving step onset time.
#' @param step_offset Numeric scalar giving step offset time.
#' @param baseline_win Length of the pre-event baseline window. Internally, the
#'   latency baseline window is restricted to be strictly smaller than
#'   \code{step_offset - step_onset}.
#' @param on_win Length of the ON analysis window.
#' @param off_win Length of the OFF analysis window.
#' @param profile_polarity Unit-level summarized polarity call.
#' @param profile_bias Unit-level summarized bias call.
#' @param latency_alpha Quantile used to define the ISI threshold.
#'
#' @return For \code{step_stim_response_isi_latency()}, a named list with element
#'   \code{ISI_Latency}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_isi_latency <- function(x,
                                           step_onset,
                                           step_offset,
                                           baseline_win = 2,
                                           on_win = 0.4,
                                           off_win = on_win,
                                           profile_polarity = c("positive", "negative", "ambiguous"),
                                           profile_bias = c("onset", "offset", "ambiguous"),
                                           latency_alpha = 0.05) {

  x <- sort(to_num(x)[is.finite(to_num(x))])
  step_onset <- to_num(step_onset)
  step_offset <- to_num(step_offset)
  baseline_win <- to_num(baseline_win)
  on_win <- to_num(on_win)
  off_win <- to_num(off_win)
  latency_alpha <- to_num(latency_alpha)

  eps <- sqrt(.Machine$double.eps)
  step_win <- step_offset - step_onset
  latency_baseline_win <- min(baseline_win, max(step_win - eps, eps))

  latency_one <- function(event, postduration, direction = c("decrease", "increase")) {
    direction <- match.arg(direction)

    pre <- x[x >= (event - latency_baseline_win) & x < event]
    if (length(pre) >= 2L) {
      preISIs <- diff(pre)
    } else if (length(pre) == 1L) {
      preISIs <- latency_baseline_win / 2
    } else {
      preISIs <- latency_baseline_win
    }

    if (!length(preISIs) || !all(is.finite(preISIs))) {
      return(NA_real_)
    }

    low_q  <- as.numeric(quantile(preISIs, probs = latency_alpha, names = FALSE))
    high_q <- as.numeric(quantile(preISIs, probs = 1 - latency_alpha, names = FALSE))
    med    <- median(preISIs)

    ISIthr <- if (direction == "decrease") {
      max(low_q, 0.60 * med)
    } else {
      min(high_q, med / 0.60)
    }

    rel <- x[x >= event & x < (event + postduration)] - event
    if (length(rel) < 2L) {
      return(NA_real_)
    }

    postISIs <- diff(rel)
    idx <- if (direction == "decrease") {
      which(postISIs < ISIthr)[1]
    } else {
      which(postISIs > ISIthr)[1]
    }

    if (is.finite(idx)) rel[idx + 1] else NA_real_
  }

  lat <- if (identical(profile_polarity, "negative")) {
    latency_one(step_onset, on_win, "increase")
  } else if (identical(profile_polarity, "positive")) {
    if (identical(profile_bias, "onset")) {
      latency_one(step_onset, on_win, "decrease")
    } else if (identical(profile_bias, "offset")) {
      latency_one(step_offset, off_win, "decrease")
    } else {
      vals <- c(
        latency_one(step_onset, on_win, "decrease"),
        latency_one(step_offset, off_win, "decrease")
      )
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }
  } else {
    NA_real_
  }

  list(
    ISI_Latency = as.numeric(lat)
  )
}

#' Farrow-style transience index
#'
#' \strong{step_stim_response_transience:} low-level function to compute transience from summarized polarity and bias calls
#'
#' @param x Numeric vector of spike times for a single spike train.
#' @param step_onset Numeric scalar giving step onset time.
#' @param step_offset Numeric scalar giving step offset time.
#' @param baseline_rate Scalar baseline rate in spikes/s.
#' @param on_win Length of the ON analysis window.
#' @param off_win Length of the OFF analysis window.
#' @param profile_polarity Unit-level summarized polarity call.
#' @param profile_bias Unit-level summarized bias call.
#' @param binwidth Bin width for the response trace.
#' @param smooth Logical; whether to smooth the trace.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'
#' @return For \code{step_stim_response_transience()}, a named list with element
#'   \code{Farrow_TransienceIndex}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_transience <- function(x,
                                          step_onset,
                                          step_offset,
                                          baseline_rate,
                                          on_win = 0.4,
                                          off_win = on_win,
                                          profile_polarity = c("positive", "negative", "ambiguous"),
                                          profile_bias = c("onset", "offset", "ambiguous"),
                                          binwidth = 0.02,
                                          smooth = FALSE,
                                          gauss_sd = 0.02) {

  x <- sort(to_num(x)[is.finite(to_num(x))])
  step_onset <- to_num(step_onset)
  step_offset <- to_num(step_offset)
  baseline_rate <- to_num(baseline_rate)
  on_win <- to_num(on_win)
  off_win <- to_num(off_win)
  binwidth <- to_num(binwidth)
  gauss_sd <- to_num(gauss_sd)

  transience_one <- function(t0, win, mode = c("positive", "negative")) {
    mode <- match.arg(mode)

    tr <- make_rate_trace(
      spk = x,
      t0 = t0,
      win = win,
      binwidth = binwidth,
      smooth = smooth,
      gauss_sd = gauss_sd
    )

    trace <- tr$rate - baseline_rate
    mod <- if (mode == "positive") pmax(trace, 0) else -pmin(trace, 0)

    mx <- suppressWarnings(max(mod, na.rm = TRUE))
    if (!is.finite(mx) || mx <= 0) {
      return(0)
    }

    mod_norm <- mod / mx
    sum(mod_norm * diff(tr$breaks), na.rm = TRUE) / win
  }

  ti <- if (identical(profile_polarity, "negative")) {
    transience_one(step_onset, on_win, "negative")
  } else if (identical(profile_polarity, "positive")) {
    if (identical(profile_bias, "onset")) {
      transience_one(step_onset, on_win, "positive")
    } else if (identical(profile_bias, "offset")) {
      transience_one(step_offset, off_win, "positive")
    } else {
      vals <- c(
        transience_one(step_onset, on_win, "positive"),
        transience_one(step_offset, off_win, "positive")
      )
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }
  } else {
    NA_real_
  }

  list(
    Farrow_TransienceIndex = as.numeric(ti)
  )
}

#' Liang-style ON/OFF preference index
#'
#' \strong{step_stim_response_liang_index:} low-level function to compute Liang-style ON/OFF preference from summarized polarity
#'
#' @param x Numeric vector of spike times for a single spike train.
#' @param step_onset Numeric scalar giving step onset time.
#' @param step_offset Numeric scalar giving step offset time.
#' @param baseline_rate Scalar baseline rate in spikes/s.
#' @param profile_polarity Unit-level summarized polarity call.
#' @param binwidth Bin width for the response trace.
#' @param smooth Logical; whether to smooth the trace.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'
#' @return For \code{step_stim_response_liang_index()}, a named list with element
#'   \code{Liang_OnOff_PreferenceIndex}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_liang_index <- function(x,
                                           step_onset,
                                           step_offset,
                                           baseline_rate,
                                           profile_polarity = c("positive", "negative", "ambiguous"),
                                           binwidth = 0.02,
                                           smooth = FALSE,
                                           gauss_sd = 0.02) {

  x <- sort(to_num(x)[is.finite(to_num(x))])
  step_onset <- to_num(step_onset)
  step_offset <- to_num(step_offset)
  baseline_rate <- to_num(baseline_rate)
  binwidth <- to_num(binwidth)
  gauss_sd <- to_num(gauss_sd)

  if (identical(profile_polarity, "ambiguous")) {
    return(list(Liang_OnOff_PreferenceIndex = NA_real_))
  }

  eps <- sqrt(.Machine$double.eps)
  step_win <- step_offset - step_onset

  on_trace <- make_rate_trace(
    spk = x,
    t0 = step_onset,
    win = step_win,
    binwidth = binwidth,
    smooth = smooth,
    gauss_sd = gauss_sd
  )$rate - baseline_rate

  off_trace <- make_rate_trace(
    spk = x,
    t0 = step_offset,
    win = step_win,
    binwidth = binwidth,
    smooth = smooth,
    gauss_sd = gauss_sd
  )$rate - baseline_rate

  if (identical(profile_polarity, "negative")) {
    on_trace[on_trace > 0] <- 0
    off_trace[off_trace > 0] <- 0
  }

  w <- abs(on_trace) + abs(off_trace)
  den_pref <- sum(w^2, na.rm = TRUE)

  idx <- if (is.finite(den_pref) && den_pref > eps) {
    sum((on_trace - off_trace) * w, na.rm = TRUE) / den_pref
  } else {
    NA_real_
  }

  if (identical(profile_polarity, "negative") && is.finite(idx)) {
    idx <- -abs(idx)
  }

  list(
    Liang_OnOff_PreferenceIndex = as.numeric(idx)
  )
}



#' Low-level response-profile metrics for a single spike train
#'
#' \strong{step_stim_response_polarity_bias}: Wrapper that computes AUC components, polarity metrics, and bias
#' metrics for one spike train.
#'
#' This helper is restricted to per-run calculations. It does not compute any
#' downstream metrics whose definition depends on run-averaged
#' \code{PolarityIndex}, \code{BiasIndex}, or their corresponding significance
#' results.
#'
#' @param x Numeric vector of spike or event times for a single spike train.
#' @param step_onset Numeric or \code{units} scalar giving step onset time.
#' @param step_offset Numeric or \code{units} scalar giving step offset time.
#' @param on_win Length of the ON analysis window. Default: \code{0.4}.
#' @param off_win Length of the OFF/post-step analysis window.
#'   Default: \code{on_win}.
#' @param baseline_win Length of the pre-step baseline window immediately
#'   preceding \code{step_onset}. Default: \code{2}.
#' @param binwidth Bin width used for the rate traces and AUC calculations.
#'   Default: \code{0.02}.
#' @param smooth Logical; if \code{TRUE}, smooth the binned traces with a
#'   Gaussian kernel before AUC calculation. Default: \code{FALSE}.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'   Default: \code{0.02}.
#'
#' @return for step_stim_response_polarity_bias: A named list containing the outputs of
#'   \code{step_stim_response_auc_components()},
#'   \code{step_stim_response_polarity()}, and
#'   \code{step_stim_response_bias()}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_polarity_bias <- function(x,
                                               step_onset,
                                               step_offset,
                                               on_win = 0.4,
                                               off_win = on_win,
                                               baseline_win = 2,
                                               binwidth = 0.02,
                                               smooth = FALSE,
                                               gauss_sd = 0.02) {

  auc_metrics <- step_stim_response_auc_components(
    x = x,
    step_onset = step_onset,
    step_offset = step_offset,
    on_win = on_win,
    off_win = off_win,
    baseline_win = baseline_win,
    binwidth = binwidth,
    smooth = smooth,
    gauss_sd = gauss_sd
  )

  c(
    auc_metrics,
    step_stim_response_polarity(auc_metrics),
    step_stim_response_bias(auc_metrics)
  )
}



#' Wrapper for conditional metrics
#'
#' \strong{step_stim_response_dependent_metrics:} Low-level wrapper function to compute downstream per-run metrics conditional on summarized polarity and bias
#'
#' @param x Numeric vector of spike times for a single spike train.
#' @param step_onset Numeric scalar giving step onset time.
#' @param step_offset Numeric scalar giving step offset time.
#' @param baseline_win Length of the baseline window immediately preceding
#'   \code{step_onset}.
#' @param on_win Length of the ON analysis window.
#' @param off_win Length of the OFF analysis window.
#' @param binwidth Bin width for the response trace.
#' @param smooth Logical; whether to smooth the trace.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#' @param baseline_rate Scalar baseline rate in spikes/s.
#' @param on_auc_pos Positive ON AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param on_auc_neg Negative ON AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param off_auc_pos Positive OFF AUC from
#'   \code{step_stim_response_auc_components()}.
#' @param profile_polarity Unit-level summarized polarity call.
#' @param profile_bias Unit-level summarized bias call.
#' @param latency_alpha Quantile used to define the latency threshold.
#'
#' @return For \code{step_stim_response_dependent_metrics()}, a named list with
#'   elements \code{Shinomoto_Baseline_LV}, \code{Response_Amplitude},
#'   \code{ISI_Latency}, \code{Farrow_TransienceIndex}, and
#'   \code{Liang_OnOff_PreferenceIndex}.
#'
#' @rdname StepStimResponseProfile
#' @export
step_stim_response_dependent_metrics <- function(x,
                                                 step_onset,
                                                 step_offset,
                                                 baseline_win = 2,
                                                 on_win = 0.4,
                                                 off_win = on_win,
                                                 binwidth = 0.02,
                                                 smooth = FALSE,
                                                 gauss_sd = 0.02,
                                                 baseline_rate,
                                                 on_auc_pos,
                                                 on_auc_neg,
                                                 off_auc_pos,
                                                 profile_polarity = c("positive", "negative", "ambiguous"),
                                                 profile_bias = c("onset", "offset", "ambiguous"),
                                                 latency_alpha = 0.05) {

  c(
    step_stim_response_baseline_lv(
      x = x,
      step_onset = step_onset,
      baseline_win = baseline_win
    ),
    step_stim_response_amplitude(
      on_auc_pos = on_auc_pos,
      on_auc_neg = on_auc_neg,
      off_auc_pos = off_auc_pos,
      on_win = on_win,
      off_win = off_win,
      profile_polarity = profile_polarity,
      profile_bias = profile_bias
    ),
    step_stim_response_isi_latency(
      x = x,
      step_onset = step_onset,
      step_offset = step_offset,
      baseline_win = baseline_win,
      on_win = on_win,
      off_win = off_win,
      profile_polarity = profile_polarity,
      profile_bias = profile_bias,
      latency_alpha = latency_alpha
    ),
    step_stim_response_transience(
      x = x,
      step_onset = step_onset,
      step_offset = step_offset,
      baseline_rate = baseline_rate,
      on_win = on_win,
      off_win = off_win,
      profile_polarity = profile_polarity,
      profile_bias = profile_bias,
      binwidth = binwidth,
      smooth = smooth,
      gauss_sd = gauss_sd
    ),
    step_stim_response_liang_index(
      x = x,
      step_onset = step_onset,
      step_offset = step_offset,
      baseline_rate = baseline_rate,
      profile_polarity = profile_polarity,
      binwidth = binwidth,
      smooth = smooth,
      gauss_sd = gauss_sd
    )
  )
}




#' Build a binned response-rate trace for a fixed analysis window
#'
#' Internal helper used by \code{\link{step_stim_resp_metrics}}.
#'
#' @param spk Numeric vector of spike times.
#' @param t0 Window start time.
#' @param win Window length.
#' @param binwidth Bin width for the histogram-based rate trace. If \code{NA},
#'   the helper returns \code{list(breaks = NULL, rate = NULL)}.
#' @param smooth Logical; if \code{TRUE}, smooth the rate trace with a Gaussian
#'   kernel.
#' @param gauss_sd Standard deviation of the Gaussian smoothing kernel.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{breaks}}{Histogram breakpoints.}
#'   \item{\code{rate}}{Binned spike rate in spikes/s.}
#' }
#'
#' @keywords internal
make_rate_trace <- function(spk, t0, win, binwidth, smooth = FALSE, gauss_sd = 0.02) {
  if (is.na(binwidth)) {
    return(list(breaks = NULL, rate = NULL))
  }
  stopifnot(is.finite(win), win > 0, is.finite(binwidth), binwidth > 0)

  t1 <- t0 + win
  spk <- spk[is.finite(spk) & spk >= t0 & spk < t1]

  breaks <- seq(t0, t1, by = binwidth)
  if (length(breaks) < 2L || tail(breaks, 1) < (t1 - 1e-12)) {
    breaks <- c(breaks, t1)
  }
  breaks <- unique(breaks)

  h <- hist(
    spk,
    breaks = breaks,
    right = FALSE,
    plot = FALSE,
    include.lowest = TRUE
  )
  rate <- h$counts / diff(h$breaks)

  if (smooth && length(rate)) {
    k_idx  <- seq(-3 * gauss_sd, 3 * gauss_sd, by = binwidth)
    kernel <- dnorm(k_idx, sd = gauss_sd)
    kernel <- kernel / sum(kernel)

    r_sm <- filter(rate, kernel, sides = 2, circular = FALSE)
    if (anyNA(r_sm)) {
      nz <- which(!is.na(r_sm))
      if (length(nz)) {
        left_gap  <- seq_len(min(nz) - 1L)
        right_gap <- if (max(nz) < length(r_sm)) (max(nz) + 1L):length(r_sm) else integer(0)
        if (length(left_gap))  r_sm[left_gap]  <- r_sm[min(nz)]
        if (length(right_gap)) r_sm[right_gap] <- r_sm[max(nz)]
      } else {
        r_sm <- rate
      }
    }
    rate <- as.numeric(r_sm)
  }

  list(breaks = h$breaks, rate = as.numeric(rate))
}

#' Compute mean spike rate in a fixed analysis window
#'
#' Internal helper used by \code{\link{step_stim_resp_metrics}}.
#'
#' @param spk Numeric vector of spike times.
#' @param t0 Window start time.
#' @param win Window length.
#' @param binwidth Bin width for the optional histogram-based computation.
#'   If \code{NA}, the rate is computed directly as count divided by window
#'   length.
#' @param gauss_sd Standard deviation of the Gaussian kernel. Included for
#'   interface consistency with \code{make_rate_trace()}, but smoothing is
#'   disabled internally for this mean-rate calculation.
#'
#' @return Numeric scalar giving mean spike rate in spikes/s.
#'
#' @keywords internal
window_mean_rate <- function(spk, t0, win, binwidth, gauss_sd = 0.02) {
  if (is.na(binwidth)) {
    return(count_in(spk, t0, t0 + win) / win)
  }

  tr <- make_rate_trace(
    spk = spk,
    t0 = t0,
    win = win,
    binwidth = binwidth,
    smooth = FALSE,
    gauss_sd = gauss_sd
  )

  if (is.null(tr$rate) || !length(tr$rate)) {
    return(0)
  }

  sum(tr$rate * diff(tr$breaks), na.rm = TRUE) / win
}
#' @keywords internal
#' @noMd
#' @param x for step_spikerate: Numeric vector of spike times already restricted to the relevant
#'   analysis window.
#' @param t for step_spikerate: Duration of the analysis window in seconds, used to convert spike
#'   count to mean spike rate.#'
#' @return For \code{step_spikerate()}, a numeric scalar giving mean spike rate
#'   in spikes/s, or \code{NA_real_} if the window duration is invalid.
step_spikerate <- function(x, t) {
  t <- to_num(t)

  if (!is.finite(t) || t <= 0)
    return(NA_real_)

  ra <- length(x) / t
  if (!is.finite(ra))
    NA_real_
  else
    as.numeric(ra)
}

#' left Inclusive, right exclusive window counting
#' @keywords internal
#' @noMd
#' @param x for step_spikerate: Numeric vector of spike times already restricted to the relevant
#'   analysis window.
#' @param a,b Left () and right margins of window
count_in <- function(x, a, b) {
  if (!length(x))
    return(0L)
  sum(x >= a & x < b)
}

