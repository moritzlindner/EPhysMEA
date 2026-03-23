#' Low-level: Fisher-z mean of pairwise trial correlations
#'
#' Computes all pairwise Pearson correlations between trials (zero-lag),
#' averages them on the Fisher-z scale, and returns the summary.
#'
#' @param trials_by_time Numeric matrix [trials x time] for one channel
#' @importFrom stats cor
#' @return A list with elements:
#'   \item{mean_trial_r}{Fisher-z averaged correlation (back-transformed)}
#'   \item{r_vals}{Numeric vector of upper-tri pairwise correlations}
#'   \item{z_mean}{Mean Fisher-z}
#'   \item{n_pairs}{Number of trial pairs}
#' @keywords internal
trial_cor_fisher_z <- function(trials_by_time, return_pairs = FALSE) {
  stopifnot(is.matrix(trials_by_time), is.numeric(trials_by_time))
  n_trials <- nrow(trials_by_time)
  if (n_trials < 2L) stop("Need at least 2 trials to compute correlations.")

  # Simple path: full correlation matrix of trials (zero-lag)
  R <- cor(t(trials_by_time), use = "pairwise.complete.obs")
  r_vals <- R[upper.tri(R, diag = FALSE)]

  # Fisher-z averaging with clamping to avoid atanh(±1)
  eps <- 1e-12
  r_vals <- pmax(pmin(r_vals, 1 - eps), -1 + eps)
  z_vals <- atanh(r_vals)
  z_mean <- mean(z_vals)

  if (isTRUE(return_pairs)) {
    list(mean_trial_r = tanh(z_mean), r_vals = r_vals, z_mean = z_mean, n_pairs = length(r_vals))
  } else {
    list(mean_trial_r = tanh(z_mean), n_pairs = length(r_vals))
  }
}

#' Low-level: Permutation p-value by random within-trial shuffling
#'
#' Builds a null distribution by randomly permuting the time-bin order within
#' each trial (destroys time locking but preserves each trial's marginal
#' distribution). Returns a one-sided p-value Pr(null >= observed).
#'
#' @param trials_by_time Numeric matrix [trials x time]
#' @param n_perm Integer number of permutations (default 1000)
#' @param seed Optional integer RNG seed for reproducibility
# @useDynLib HetDataAnalysis, .registration=TRUE
# @importFrom Rcpp sourceCpp
#' @return List with observed mean_trial_r, p_value, and the vector of permuted mean_trial_r
#' @importFrom cli cli_progress_done cli_progress_update cli_progress_bar
#' @keywords internal
# Low-level: Permutation p-value by random within-trial shuffling (Monte Carlo with early stopping)
trial_cor_pvalue_perm <- function(trials_by_time,
                                  n_perm = 1000L,
                                  seed = NULL,
                                  return_perm = FALSE,
                                  min_perm = 75L,
                                  stop_lower = 0.005,
                                  stop_upper = 0.10,
                                  return_pairs = FALSE) {
  stopifnot(is.matrix(trials_by_time), is.numeric(trials_by_time))
  n_trials <- nrow(trials_by_time)
  n_time   <- ncol(trials_by_time)
  if (n_trials < 2L || n_time < 2L) stop("Need at least 2 trials and 2 time bins.")
  if (!is.null(seed)) set.seed(seed)

  # Observed statistic from unshuffled data
  obs <- trial_cor_fisher_z(trials_by_time, return_pairs = return_pairs) #cpp_trial_cor_fisher_z
  obs_r   <- obs$mean_trial_r
  n_pairs <- obs$n_pairs

  exceed <- 0L
  used   <- 0L
  if (isTRUE(return_perm)) perm_stats <- numeric(n_perm)
  # pre-allocate once outside the loop
  shuffled <- matrix(NA_real_, n_trials, n_time)

  cli_progress_bar("Calculating p-value for trial-to-trial correlation.", total = n_perm)
  for (b in seq_len(n_perm)) {
    idx <- vapply(seq_len(n_trials), function(i) {
      repeat {
        p <- sample.int(n_time)
        if (!all(p == seq_len(n_time))) return(p)  # reject identity
      }
    }, integer(n_time))

    rows <- rep(seq_len(n_trials), each = n_time)  # 1,1,...,1, 2,2,...,2, ...
    cols <- as.vector(idx)                          # c(idx[,1], idx[,2], ..., idx[,n_trials])

    sel <- trials_by_time[cbind(rows, cols)]
    shuffled <- matrix(sel, nrow = n_trials, ncol = n_time, byrow = TRUE)

    stat_b <- trial_cor_fisher_z(shuffled, return_pairs = FALSE)$mean_trial_r

    used <- b
    if (isTRUE(return_perm)) perm_stats[b] <- stat_b
    if (is.finite(stat_b) && stat_b >= obs_r) exceed <- exceed + 1L

    if (used >= min_perm) {
      p_hat <- (exceed + 1) / (used + 1)
      if (p_hat < stop_lower || p_hat > stop_upper) break
    }
    cli_progress_update()
  }
  cli_progress_done()

  p_value <- (exceed + 1) / (used + 1)

  if (isTRUE(return_perm)) {
    list(mean_trial_r = obs_r,
         p_value = p_value,
         n_pairs = n_pairs,
         n_perm_used = used,
         permuted = perm_stats[seq_len(used)])
  } else {
    list(mean_trial_r = obs_r,
         p_value = p_value,
         n_pairs = n_pairs,
         n_perm_used = used)
  }
}

#' TrialCorrelationTest generic
#'
#' Compute per-channel similarity of repeated time-locked runs by Pearson
#' correlations summarized on the Fisher-z scale, with Monte Carlo p-values
#' from within-run random shuffling.
#'
#' @param X An \code{EPhysContinuous} object.
#' @param n_perm Integer; maximum number of permutations (default 1000). Early stopping may end sooner.
#' @param seed Optional integer RNG seed for reproducibility.
#' @param min_perm Integer; minimum permutations before early stopping is considered (default 75).
#' @param stop_lower Numeric; lower running p-hat threshold for early stopping (default 0.005).
#' @param stop_upper Numeric; upper running p-hat threshold for early stopping (default 0.10).
#' @param return_perm Logical; currently retained for internal use.
#' @inheritParams EPhysData::group_apply
#'
#' @details
#' For \code{EPhysContinuous} input, repeated trials are defined as repeated
#' runs within the same \code{RecordingID}. For each
#' \code{RecordingID} and \code{Channel}, the method compares all usable runs
#' pairwise. After \code{\link{Bin}}, runs remain separate; binning does not
#' average or merge them. The \code{Repeat} column is not used explicitly by
#' this method.
#'
#' @return
#' A data.frame with one row per (\code{RecordingID}, \code{Channel}) containing:
#' \describe{
#'   \item{RecordingID}{Character; recording identifier copied from \code{Metadata(X)}.}
#'   \item{Channel}{Character; channel name from \code{Channels(X)}.}
#'   \item{mean_trial_r}{Numeric; Fisher-z-averaged zero-lag Pearson correlation
#'     between all usable pairs of repeated runs for this recording and channel.}
#'   \item{p_value}{Numeric; one-sided Monte Carlo permutation p-value
#'     \eqn{\mathrm{Pr}(r_{\mathrm{perm}} \ge r_{\mathrm{obs}})}.}
#'   \item{n_trials}{Integer; number of usable runs/trials entering the correlation.}
#'   \item{n_pairs}{Integer; number of distinct run/trial pairs contributing to \code{mean_trial_r}.}
#'   \item{n_perm_used}{Integer; number of permutations actually evaluated.}
#' }
#' @export
setGeneric("TrialCorrelationTest", function(X,
                                            n_perm      = 1000L,
                                            seed        = NULL,
                                            min_perm    = 75L,
                                            stop_lower  = 0.005,
                                            stop_upper  = 0.10,
                                            return_perm = FALSE,
                                            parallel    = !interactive(),
                                            error       = c("stop","warn")[1],
                                            progress    = interactive())
  standardGeneric("TrialCorrelationTest"))

#' @rdname TrialCorrelationTest
#' @importFrom stats cor var
#' @importFrom matrixStats colVars
setMethod("TrialCorrelationTest",
          signature(X = "EPhysContinuous"),
          function(X,
                   n_perm      = 1000L,
                   seed        = NULL,
                   min_perm    = 75L,
                   stop_lower  = 0.005,
                   stop_upper  = 0.10,
                   return_perm = FALSE,
                   parallel    =  !interactive(),
                   error       = c("stop","warn")[1],
                   progress    = interactive()) {

            if (!is.null(seed)) set.seed(seed)
            err_mode <- match.arg(error, c("stop","warn"))

            # worker over each (RecordingID, Channel)
            # mat: [time × runs_in_recording] -> returns a named numeric vector of metrics
            worker <- function(mat, n_perm, min_perm, stop_lower, stop_upper, return_perm) {
              n_runs <- ncol(mat)
              if (n_runs < 2L)
                return(c(mean_trial_r=0, p_value=1, n_trials=n_runs, n_pairs=0, n_perm_used=0))

              keep <- (colSums(is.na(mat)) == 0L) & (colVars(mat) > 0)
              if (!any(keep))
                return(c(mean_trial_r=0, p_value=1, n_trials=0, n_pairs=0, n_perm_used=0))

              Xtr <- t(mat[, keep, drop=FALSE])  # trials × time, only once
              n_trials <- nrow(Xtr)
              if (n_trials < 2L)
                return(c(mean_trial_r=0, p_value=1, n_trials=n_trials, n_pairs=0, n_perm_used=0))

              res <- trial_cor_pvalue_perm(
                Xtr,
                n_perm      = n_perm,
                seed        = NULL,      # seeding handled outside if desired
                return_perm = return_perm,
                min_perm    = min_perm,
                stop_lower  = stop_lower,
                stop_upper  = stop_upper
              )
              c(mean_trial_r = res$mean_trial_r, p_value = res$p_value,
                n_trials = n_trials, n_pairs = res$n_pairs, n_perm_used = res$n_perm_used)
            }

            # 1) compute per (RecordingID, Channel)
            nested <- group_apply(
              X, worker,
              parallel = parallel, error = err_mode, progress = progress,
              n_perm = n_perm, min_perm = min_perm,
              stop_lower = stop_lower, stop_upper = stop_upper,
              return_perm = return_perm
            )
            # nested[[RecordingID]][[Channel]] -> named numeric vector

            # 2) assemble tidy data.frame (RecordingID, RunUID, Channel, metrics)
            md <- Metadata(X)
            rec_levels <- unique(md$RecordingID)
            keep_rows  <- match(rec_levels, md$RecordingID)

            rec_ids <- names(nested)
            ch_names <- Channels(X)

            rows <- lapply(rec_ids, function(rec) {
              ch_list <- nested[[rec]]
              do.call(rbind, lapply(ch_names, function(ch) {
                v <- ch_list[[ch]]
                if (is.null(v)) v <- c(mean_trial_r = NA_real_, p_value = NA_real_,
                                       n_trials = NA_real_, n_pairs = NA_real_, n_perm_used = NA_real_)
                as.data.frame(
                  c(list(RecordingID = rec, Channel = ch), as.list(v)),
                  check.names = FALSE, stringsAsFactors = FALSE
                )
              }))
            })

            out <- do.call(rbind, rows)
            rownames(out) <- NULL
            out
          }
)
