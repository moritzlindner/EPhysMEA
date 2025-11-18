#' Per-channel location metrics from EPhysContinuous
#'
#' Computes per-channel location metrics from a trial-averaged
#' \code{EPhysContinuous} object containing RayleighZ (or similar) values.
#'
#' For each channel:
#' \itemize{
#'   \item \strong{peak_seg}: segment (from \code{TimeTrace}) at maximum value
#'   \item \strong{peak_value}: value at the peak
#'   \item \strong{auc_norm}: peak-normalized AUC over segments, i.e.
#'         the mean of \eqn{y/peak} across segments (with \code{NA} treated as 0).
#' }
#'
#' @inheritParams Normalize
#' @param span Ignored (kept for backward compatibility).
#' @inheritParams EPhysData::`lapply-EPhys`
#'
#' @return A \code{data.frame} with one row per channel and columns:
#'   \code{channel}, \code{peak_seg}, \code{peak_value}, \code{auc_norm}, \code{unit}.
#'
#' @importClassesFrom EPhysData EPhysContinuous
#' @importFrom EPhysData Channels TimeTrace TimeUnits nested2df lapply
#' @importFrom reshape2 melt dcast
#' @export
setGeneric("LocationMetrics", function(X,
                                       span = 0.3,
                                       parallel = !interactive(),
                                       progress = interactive(),
                                       lapply_error = c("stop", "warn")[1]) {
  standardGeneric("LocationMetrics")
})

#' @export
#' @export
setMethod("LocationMetrics", signature(X = "EPhysContinuous"),
          function(X,
                   span = 0.3,
                   parallel = !interactive(),
                   progress = interactive(),
                   lapply_error = c("stop", "warn")[1]) {
            st <- StimulusTrace(X)
            ch <- Channels(X)

            # Worker applied per (run, channel): receives a numeric time-series "y"
            per_channel_fun <- function(y) {
              y[is.na(y)] <- 0
              peak_idx <- which.max(y)
              peak_seg <- st[peak_idx]
              peak_val <- y[peak_idx]

              if (is.finite(peak_val) && peak_val > 0) {
                auc_norm <- mean(y / peak_val)
                if (is.finite(auc_norm)) {
                  auc_norm <- max(0, min(1, auc_norm))
                } else {
                  auc_norm <- NA_real_
                  peak_idx <- NA_real_
                }
              } else {
                auc_norm <- NA_real_
                peak_idx <- NA_real_
              }

              list(
                peak_seg = peak_seg,
                peak_value = as.vector(peak_val),
                auc_norm = auc_norm
              )
            }

            # Compute per (run,channel)
            res_nested <- lapply(
              X        = X,
              FUN      = per_channel_fun,
              parallel = parallel,
              error    = match.arg(lapply_error, c("stop", "warn")),
              progress = progress
            )
            nested2df(X,res_nested)
          })
