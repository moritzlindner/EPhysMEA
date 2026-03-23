#' #' Segment-level variability calculation (internal)
#' #'
#' #' Computes variability statistics of spike trains within predefined
#' #' intervals (`breaks`) for one trial × one channel.
#' #'
#' #' Supports four metric families:
#' #' \itemize{
#' #'   \item \code{"lv"}: Local variation (Shinomoto LV), anchored at interior spikes.
#' #'   \item \code{"cv2"}: Local CV2 (|CV2-1|), anchored at interior spikes.
#' #'   \item \code{"cv"}: Coefficient of variation of ISIs per segment.
#' #'   \item \code{"mad"}: Median absolute deviation of ISIs per segment.
#' #' }
#' #'
#' #' @param times Numeric vector of spike times (sorted, finite).
#' #' @param breaks Numeric vector of interval cut points (see [base::cut()]).
#' #' @param metric Character; one of \code{"lv"}, \code{"cv2"}, \code{"cv"}, \code{"mad"}.
#' #' @param labels Optional labels for the segments.
#' #' @param right Logical; whether intervals are closed on the right.
#' #' @param include.lowest Logical; whether to include the lowest break value.
#' #' @param p_values Logical; if \code{TRUE}, compute KS-test p-values for LV / CV2 metrics.
#' #' @param ... Passed to [base::cut()].
#' #'
#' #' @return A data.frame with columns:
#' #' \describe{
#' #'   \item{seg_id}{Integer index of the segment.}
#' #'   \item{seg_start, seg_end}{Segment boundaries.}
#' #'   \item{seg_label}{Label of the segment.}
#' #'   \item{n_items}{Number of items (spikes or ISIs) used.}
#' #'   \item{value}{Aggregate variability value.}
#' #'   \item{p_value}{Optional KS-test p-value; \code{NA} otherwise.}
#' #' }
#' #'
#' #' @keywords internal
#' segment_variability_by_breaks <- function(times,
#'                                           breaks,
#'                                           metric = c("lv", "cv2", "cv", "mad"),
#'                                           labels = NULL,
#'                                           right = TRUE,
#'                                           include.lowest = FALSE,
#'                                           p_values = FALSE,
#'                                           ...) {
#'   stopifnot(is.numeric(times), is.numeric(breaks))
#'   times  <- sort(times[is.finite(times)])
#'   metric <- match.arg(metric)
#'
#'   n_seg <- length(breaks) - 1L
#'   if (n_seg < 1L) stop("'breaks' must define at least one interval.")
#'
#'   # ---- construct "local items" and their anchors depending on metric family ----
#'   if (metric %in% c("lv", "cv2")) {
#'     # adjacent-pair metrics -> anchor at interior spikes (length n-2)
#'     if (metric == "lv") {
#'       local <- spike_irreg_lv_shinomoto_local(times)     # length n-2
#'     } else {
#'       local <- spike_irreg_cv2_local(times)              # |CV2-1|, length n-2
#'     }
#'     anchor <- if (length(times) >= 3L) times[2:(length(times)-1L)] else numeric(0)
#'
#'     agg_fun <- if (metric == "lv") {
#'       function(x) {
#'         if (!length(x)) return(NA_real_)
#'         m <- mean(x, na.rm = TRUE)
#'         if (is.na(m)) 1 else m
#'       }
#'     } else {
#'       function(x) if (!length(x)) NA_real_ else stats::median(x, na.rm = TRUE)
#'     }
#'
#'     # p-values available only for lv/cv2 (vs Poissonized null)
#'     pval_fun <- if (isTRUE(p_values)) {
#'       if (metric == "lv") {
#'         function(x) {
#'           d <- abs(sqrt((4/3) * x) - 1)          # map LV -> |CV2-1| ∈ [0,1]
#'           d <- d[is.finite(d)]
#'           if (!length(d)) return(NA_real_)
#'           suppressWarnings(stats::ks.test(d, "punif", min = 0, max = 1)$p.value)
#'         }
#'       } else { # cv2 already = |CV2-1|
#'         function(x) {
#'           d <- x[is.finite(x)]
#'           if (!length(d)) return(NA_real_)
#'           suppressWarnings(stats::ks.test(d, "punif", min = 0, max = 1)$p.value)
#'         }
#'       }
#'     } else NULL
#'
#'   } else {
#'     # "cv" / "mad": use ISIs, anchored at ISI midpoints
#'     if (length(times) < 2L) {
#'       anchor <- numeric(0); isi <- numeric(0)
#'     } else {
#'       isi   <- diff(times)
#'       anchor <- (times[-1L] + times[-length(times)]) / 2
#'     }
#'
#'     agg_fun <- if (metric == "cv") {
#'       function(x) {                               # x are ISIs in the segment
#'         x <- x[is.finite(x)]
#'         if (length(x) < 2L) return(NA_real_)
#'         m <- mean(x); s <- stats::sd(x)
#'         if (!is.finite(m) || m == 0) NA_real_ else s / m
#'       }
#'     } else { # "mad": normalized MAD of ISIs (scale-free)
#'       function(x) {
#'         x <- x[is.finite(x)]
#'         if (!length(x)) return(NA_real_)
#'         med <- stats::median(x)
#'         if (!is.finite(med) || med == 0) return(NA_real_)
#'         stats::median(abs(x - med)) / med
#'       }
#'     }
#'
#'     local <- isi
#'     # no formal p-values for cv/mad in this helper → return NA
#'     pval_fun <- if (isTRUE(p_values)) function(x) NA_real_ else NULL
#'   }
#'
#'   f <- cut(anchor, breaks = breaks, labels = labels,
#'            right = right, include.lowest = include.lowest, ...)
#'
#'   split_vals <- split(local, f, drop = FALSE)
#'
#'   agg_vals <- vapply(split_vals, agg_fun, numeric(1))
#'   n_used   <- vapply(split_vals, function(x) sum(is.finite(x)), integer(1))
#'
#'   pval <- if (!is.null(pval_fun)) {
#'     vapply(split_vals, pval_fun, numeric(1))
#'   } else {
#'     rep(NA_real_, length(split_vals))
#'   }
#'
#'   seg_start <- breaks[-length(breaks)]
#'   seg_end   <- breaks[-1L]
#'   seg_lab   <- if (is.null(labels)) levels(f) else labels
#'
#'   data.frame(
#'     seg_id    = seq_len(n_seg),
#'     seg_start = seg_start,
#'     seg_end   = seg_end,
#'     seg_label = seg_lab,
#'     n_items   = n_used,
#'     value     = as.numeric(agg_vals),
#'     p_value   = as.numeric(pval),
#'     stringsAsFactors = FALSE
#'   )
#' }
#'
#' #' SegmentVarInternal (internal S4 method)
#' #'
#' #' Core worker for \code{\link{SegmentVar}}, applied to \link{EPhysEvents}.
#' #' Computes trial × segment × channel arrays of variability values
#' #' (and optional p-values) using \code{segment_variability_by_breaks()}.
#' #'
#' #' @inheritParams SegmentVar
#' #'
#' #' @return A list with elements:
#' #' \itemize{
#' #'   \item \code{values}: 3D array [trial × seg × channel] of metric values.
#' #'   \item \code{p_values}: (optional) 3D array [trial × seg × channel] of p-values.
#' #' }
#' #'
#' #' @importClassesFrom EPhysData EPhysEvents
#' #' @importFrom EPhysData Metadata Channels ChannelMetadata
#' #' @keywords internal
#' setGeneric("SegmentVarInternal",
#'            function(X, ...) standardGeneric("SegmentVarInternal"))
#'
#' setMethod("SegmentVarInternal", "EPhysEvents",
#'           function(X, breaks,
#'                    metric = c("lv", "cv2", "cv", "mad"),
#'                    p_values = FALSE,
#'                    labels = NULL,
#'                    right = TRUE,
#'                    include.lowest = FALSE,
#'                    ...) {
#'
#'             metric <- match.arg(metric)
#'
#'             n_tr  <- nrow(Metadata(X))
#'             ch    <- Channels(X)
#'             n_ch  <- length(ch)
#'             n_seg <- length(breaks) - 1L
#'             if (n_seg < 1L) stop("'breaks' must define at least one interval.")
#'
#'             one_trial_one_channel <- function(times_vec) {
#'               segment_variability_by_breaks(
#'                 times   = times_vec,
#'                 breaks  = breaks,
#'                 metric  = metric,
#'                 labels  = labels,
#'                 right   = right,
#'                 include.lowest = include.lowest,
#'                 p_values = p_values,
#'                 ...
#'               )
#'             }
#'
#'             # arrays [trial × seg × channel]
#'             val_arr <- array(NA_real_, dim = c(n_tr, n_seg, n_ch),
#'                              dimnames = list(
#'                                trial   = as.character(seq_len(n_tr)),
#'                                seg     = paste0("seg", seq_len(n_seg)),
#'                                channel = ch
#'                              ))
#'             p_arr <- if (isTRUE(p_values)) array(NA_real_, dim = c(n_tr, n_seg, n_ch),
#'                                                  dimnames = list(
#'                                                    trial   = as.character(seq_len(n_tr)),
#'                                                    seg     = paste0("seg", seq_len(n_seg)),
#'                                                    channel = ch
#'                                                  )) else NULL
#'
#'             for (tr in seq_len(n_tr)) {
#'               di <- X@Data[[tr]]
#'               for (j in seq_len(n_ch)) {
#'                 chname <- ch[j]
#'                 times_vec <- if (is.null(di) || is.null(di[[chname]])) numeric(0) else di[[chname]]
#'                 res <- one_trial_one_channel(times_vec)
#'                 if (nrow(res)) {
#'                   val_arr[tr, , j] <- res$value
#'                   if (isTRUE(p_values) && "p_value" %in% names(res)) {
#'                     p_arr[tr, , j] <- res$p_value
#'                   }
#'                 }
#'               }
#'             }
#'
#'             val_arr[is.na(val_arr)] <- NA
#'
#'             out <- list(values = val_arr)
#'             if (isTRUE(p_values)) out$p_values <- p_arr
#'             out
#'           })
#'
#' #' Combine p-values across trials (internal)
#' #'
#' #' Applies Fisher’s method to combine segment-wise p-values across trials.
#' #'
#' #' @param p_arr Numeric 3D array [trial × seg × channel] of p-values.
#' #'
#' #' @return A 2D matrix [seg × channel] of combined p-values.
#' #'
#' #' @details
#' #' For each segment × channel, all finite, positive p-values across trials
#' #' are combined using
#' #' \deqn{ X^2 = -2 \sum_i \ln(p_i) }
#' #' which follows a \eqn{\chi^2} distribution with \eqn{2k} degrees of freedom
#' #' (where \eqn{k} is the number of p-values).
#' #'
#' #' @keywords internal
#' combine_pvalues_fisher <- function(p_arr) {
#'   stopifnot(length(dim(p_arr)) == 3)
#'   apply(p_arr, c(2, 3), function(p) {
#'     p <- p[is.finite(p) & p > 0]
#'     if (!length(p)) return(NA_real_)
#'     stat <- -2 * sum(log(p))
#'     1 - stats::pchisq(stat, df = 2 * length(p))
#'   })
#' }
#'
#' #' Segment-wise spike train variability metrics
#' #'
#' #' Computes variability metrics of spike trains segmented by user-specified
#' #' interval boundaries (`breaks`). Returns an \link{EPhysContinuous} object
#' #' with per-trial, per-segment, and per-channel values.
#' #'
#' #' @section Metrics:
#' #' Four different metrics are available:
#' #' \describe{
#' #'   \item{"lv"}{Shinomoto local variation (LV), averaged within each segment.}
#' #'   \item{"cv2"}{Median CV2 (|CV2-1|) across interspike intervals in each segment.}
#' #'   \item{"cv"}{Coefficient of variation of ISIs within each segment
#' #'     (sd / mean of ISIs).}
#' #'   \item{"mad"}{Median absolute deviation (MAD) of ISIs, normalized by the
#' #'     median ISI.}
#' #' }
#' #'
#' #' @param X An \link{EPhysEvents} object containing spike times organized
#' #'   per trial and channel.
#' #' @param breaks Numeric vector of cut points defining the segment
#' #'   boundaries (passed to [base::cut()]).
#' #' @param metric Character; which variability metric to compute. One of
#' #'   \code{"lv"}, \code{"cv2"}, \code{"cv"}, \code{"mad"}.
#' #' @param labels Optional labels for the segments (see [base::cut()]).
#' #' @param right Logical; if \code{TRUE}, intervals are closed on the right
#' #'   (see [base::cut()]).
#' #' @param include.lowest Logical; whether to include the lowest
#' #'   break value (see [base::cut()]).
#' #' @param p_values Logical; if \code{TRUE}, additionally computes
#' #'   segment-level p-values for LV and CV2 metrics using KS tests
#' #'   against the uniform null. For \code{"cv"} and \code{"mad"},
#' #'   p-values are not available and will be returned as \code{NA}.
#' #' @param ... Passed to internal helpers and [base::cut()].
#' #'
#' #' @return An \link{EPhysContinuous} object with:
#' #' \itemize{
#' #'   \item \code{Data}: 3D array \code{[time × trial × channel]} where rows
#' #'     correspond to segment boundaries and columns to trials/channels.
#' #'   \item \code{TimeTrace}: The \code{breaks} vector, marking segment edges.
#' #'   \item \code{Metadata}, \code{Channels}, and other slots copied from \code{X}.
#' #' }
#' #'
#' #' @seealso [spike_irreg_lv_shinomoto_local()], [spike_irreg_cv2_local()] for local metrics.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' breaks <- seq(0, 2, by = 0.5)
#' #' segvar <- SegmentVar(ephys_events, breaks, metric = "cv")
#' #' plot(as.data.frame(segvar))
#' #' }
#' #'
#' setGeneric("SegmentVar",
#'            function(X, ...) standardGeneric("SegmentVar"))
#'
#' #' @export
#' setMethod("SegmentVar", "EPhysEvents",
#'           function(X, breaks,
#'                    metric = c("lv", "cv2", "cv", "mad"),
#'                    labels = NULL,
#'                    right = TRUE,
#'                    include.lowest = FALSE,
#'                    p_values = FALSE,
#'                    ...) {
#'
#'             metric <- match.arg(metric)
#'
#'             res <- SegmentVarInternal(
#'               X, breaks = breaks, metric = metric, p_values = p_values,
#'               labels = labels, right = right, include.lowest = include.lowest, ...
#'             )
#'
#'             val_arr <- res$values  # [trial × seg × channel]
#'             if (is.null(val_arr)) stop("No values computed.")
#'
#'             n_tr  <- dim(val_arr)[1]
#'             n_seg <- dim(val_arr)[2]
#'             n_ch  <- dim(val_arr)[3]
#'
#'             # Build Data array [time × trial × channel] with time = length(breaks)
#'             Data <- array(NA_real_, dim = c(length(breaks), n_tr, n_ch),
#'                           dimnames = list(
#'                             time    = paste0("t", seq_along(breaks)),
#'                             trial   = NULL,
#'                             channel = Channels(X)
#'                           ))
#'             Data[2:(n_seg + 1), , ] <- aperm(val_arr, c(2, 1, 3))  # seg→time rows 2..end
#'
#'             md <- Metadata(X)
#'             dimnames(Data)$trial <- as.character(md$RunUID)
#'
#'             TimeTrace <- as.numeric(breaks)
#'
#'             out <- newEPhysContinuous(
#'               Data             = Data,
#'               TimeTrace        = TimeTrace,
#'               Metadata         = md,
#'               Channels         = Channels(X),
#'               Channel_Metadata = X@Channel_Metadata,
#'               StimulusTrace    = numeric(0),
#'               TimeUnits        = X@TimeUnits %||% "s",
#'               StimulusUnits    = "",
#'               Imported         = X@Imported,
#'               ExamInfo         = X@ExamInfo,
#'               SubjectInfo      = X@SubjectInfo
#'             )
#'             out
#'           })
#'
#' `%||%` <- function(a, b) if (!is.null(a) && length(a) == 1L && nzchar(a)) a else b
