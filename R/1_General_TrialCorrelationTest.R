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

trial_shuffle_indices <- function(n_trials, n_time) {
  original <- seq_len(n_time)
  vapply(seq_len(n_trials), function(i) {
    repeat {
      indices <- sample.int(n_time)
      if (!identical(indices, original)) return(indices)
    }
  }, integer(n_time))
}

trial_cor_shuffled_r <- function(trials_by_time, rows, indices) {
  selected <- trials_by_time[cbind(rows, as.vector(indices))]
  shuffled <- matrix(selected, nrow = nrow(trials_by_time), byrow = TRUE)
  trial_cor_fisher_z(shuffled)$mean_trial_r
}

trial_cor_cuda_device <- function() {
  # Loading torch can consume R RNG draws before the first permutation.
  seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (!is.null(seed)) {
      assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("GPU = TRUE requires the optional 'torch' package. Install it with install.packages('torch').",
         call. = FALSE)
  }
  if (!torch::torch_is_installed() || !torch::cuda_is_available()) {
    stop("GPU = TRUE requires CUDA to be available in R torch. Check torch::cuda_is_available() or use GPU = FALSE.",
         call. = FALSE)
  }
  "cuda"
}

trial_cor_torch_evaluator <- function(trials_by_time, device) {
  n_trials <- nrow(trials_by_time)
  n_time <- ncol(trials_by_time)
  centred <- trials_by_time - rowMeans(trials_by_time)
  norms <- sqrt(rowSums(centred^2))
  if (any(!is.finite(norms) | norms == 0)) {
    stop("GPU trial correlation requires finite, non-constant trials.", call. = FALSE)
  }

  # The mean and norm of each trial stay the same under permutation.
  data <- torch::torch_tensor(centred / norms, dtype = torch::torch_double(), device = device)
  upper <- torch::torch_tensor(upper.tri(matrix(0, n_trials, n_trials)),
                               dtype = torch::torch_double(), device = device)
  n_pairs <- choose(n_trials, 2L)

  function(indices) {
    n_batch <- dim(indices)[3L]
    positions <- torch::torch_tensor(indices, dtype = torch::torch_long(), device = device)
    positions <- positions$permute(c(3L, 2L, 1L))
    data_batch <- data$unsqueeze(1L)$expand(c(n_batch, n_trials, n_time))
    shuffled <- torch::torch_gather(data_batch, dim = 3L, index = positions)
    r <- torch::torch_bmm(shuffled, shuffled$transpose(2L, 3L))

    # Match the CPU Fisher-z clipping in double precision.
    z <- r$clamp(min = -1 + 1e-12, max = 1 - 1e-12)$atanh() * upper
    scores <- (z$sum(dim = 3L)$sum(dim = 2L) / n_pairs)$tanh()
    max_abs_r <- torch::torch_amax(r$abs() * upper, dim = c(2L, 3L))
    torch::as_array(torch::torch_stack(list(scores, max_abs_r), dim = 2L)$cpu())
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
                                  return_pairs = FALSE,
                                  gpu_device = NULL) {
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
  rows <- rep(seq_len(n_trials), each = n_time)

  evaluate_batch <- if (!is.null(gpu_device)) {
    trial_cor_torch_evaluator(trials_by_time, gpu_device)
  }
  batch_end <- 0L
  # Limit memory for both the shuffled traces and the correlation matrices.
  batch_size <- max(1L, min(64L, floor(2^21 / (as.double(n_trials) * max(n_time, n_trials)))))

  cli_progress_bar("Calculating p-value for trial-to-trial correlation.", total = n_perm)
  on.exit(cli_progress_done(), add = TRUE)
  for (b in seq_len(n_perm)) {
    if (is.null(evaluate_batch)) {
      indices <- trial_shuffle_indices(n_trials, n_time)
      stat_b <- trial_cor_shuffled_r(trials_by_time, rows, indices)
    } else {
      if (b > batch_end) {
        batch_start <- b
        batch_end <- min(n_perm, b + batch_size - 1L)
        n_batch <- batch_end - batch_start + 1L
        batch_indices <- array(0L, dim = c(n_time, n_trials, n_batch))
        batch_seeds <- vector("list", n_batch)
        for (i in seq_len(n_batch)) {
          batch_indices[, , i] <- trial_shuffle_indices(n_trials, n_time)
          batch_seeds[[i]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
        batch_results <- evaluate_batch(batch_indices)
      }

      i <- b - batch_start + 1L
      stat_b <- batch_results[i, 1L]
      # Recheck close comparisons and near-perfect pair correlations with the CPU
      # statistic, where roundoff can otherwise change an exceedance decision.
      if (!is.finite(stat_b) || !is.finite(batch_results[i, 2L]) ||
          abs(stat_b - obs_r) <= 1e-8 || batch_results[i, 2L] >= 1 - 1e-6) {
        stat_b <- trial_cor_shuffled_r(trials_by_time, rows, batch_indices[, , i])
      }
    }

    used <- b
    if (isTRUE(return_perm)) perm_stats[b] <- stat_b
    if (is.finite(stat_b) && stat_b >= obs_r) exceed <- exceed + 1L

    if (used >= min_perm) {
      p_hat <- (exceed + 1) / (used + 1)
      if (p_hat < stop_lower || p_hat > stop_upper) break
    }
    cli_progress_update()
  }
  if (used < batch_end) {
    # Discard RNG draws from the unused tail of the final batch.
    assign(".Random.seed", batch_seeds[[used - batch_start + 1L]], envir = .GlobalEnv)
  }

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
#' @param GPU Logical; if \code{TRUE}, batch permutation correlations on one
#'   CUDA device using the optional \pkg{torch} package. Defaults to \code{FALSE}.
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
#' With \code{GPU = TRUE}, shuffle indices are generated by R and the batched
#' correlations run on CUDA. For comparisons with the CPU backend, use
#' \code{parallel = FALSE} with the same \code{seed}. The GPU backend processes
#' recording/channel groups in one R process, independently of \code{parallel}.
#' An explicit GPU request requires a working CUDA installation of \pkg{torch}.
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
#'   \item{n_perm_used}{Integer; number of permutations used to compute the p-value.}
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
                                            progress    = interactive(),
                                            GPU         = FALSE)
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
                   progress    = interactive(),
                   GPU         = FALSE) {

            if (!is.logical(GPU) || length(GPU) != 1L || is.na(GPU)) {
              stop("TrialCorrelationTest() requires GPU to be TRUE or FALSE.", call. = FALSE)
            }
            gpu_device <- if (GPU) trial_cor_cuda_device()
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
                stop_upper  = stop_upper,
                gpu_device  = gpu_device
              )
              c(mean_trial_r = res$mean_trial_r, p_value = res$p_value,
                n_trials = n_trials, n_pairs = res$n_pairs, n_perm_used = res$n_perm_used)
            }

            # 1) compute per (RecordingID, Channel)
            nested <- group_apply(
              X, worker,
              parallel = if (GPU) FALSE else parallel, error = err_mode, progress = progress,
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
