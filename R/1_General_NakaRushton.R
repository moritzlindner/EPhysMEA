#' Naka–Rushton parameters from EPhysContinuous (with baseline)
#'
#' Fits the 4-parameter Naka–Rushton
#' \deqn{R(I) = \mathrm{Baseline} + R_{\max}\, I^n / (I^n + k^n)}
#' per \strong{channel and Run} , using \code{TimeTrace(X)} as the (Michelson) contrast
#' and \code{X@Data} as the response. Contrasts are clamped to \eqn{[0, 1]} in
#' the internal fitter.
#'
#' @name NakaRushton
#' @aliases
#'   NakaRushton-methods
#'   NakaRushton,EPhysContinuous-method
#'
#' @inheritParams Normalize
#' @param I Numeric vector of contrasts (Michelson); values are clamped to \code{[0,1]} internally. Must contain at least three distinct finite values; repeated contrasts are allowed.
#' @param R Numeric vector of responses corresponding to \code{I} (same length and order); non-finite pairs are removed prior to fitting.
#' @param ... Currently ignored; reserved for future options.
#'
#' @details
#' The Stimulus Trace is treated as contrast. A single-channel 4-parameter Naka–Rushton is fit then fit per channel and run. The returned
#' \code{Rmax} is the amplitude above \code{Baseline} (so \code{Top = Baseline + Rmax}).
#' The semi-saturation contrast \code{C50} equals the fitted \code{k}.
#'
#' \strong{Lower-level function:}
#' \code{nakarushton(I, R)} performs the single-channel fit used internally. It
#' clamps \code{I} to \eqn{[0,1]}, requires at least three distinct contrasts,
#' fits \code{R ~ naka_rushton4(I, Baseline, Rmax, k, n)} with \code{minpack.lm::nlsLM},
#' and uses starting values and bounds from \code{.guess_nr4()}. See
#' \emph{Value} below for its return structure.
#'
#' @return A \code{data.frame} with one row per channel:
#' \describe{
#'   \item{Channel}{Channel name.}
#'   \item{Baseline}{Fitted baseline (offset).}
#'   \item{Rmax}{Amplitude above baseline.}
#'   \item{Top}{Baseline + Rmax.}
#'   \item{C50}{Semi-saturation contrast \code{k}.}
#'   \item{n}{Slope/steepness.}
#'   \item{RMSE}{RMSE on the fit points (after duplicate collapse).}
#'   \item{nRMSE_Rmax}{RMSE normalized by \code{Rmax}.}
#'   \item{R2}{Coefficient of determination.}
#'   \item{n_points}{Number of unique contrast points used in the fit.}
#'   \item{message}{Status message (e.g., \code{"ok"}, \code{"insufficient data"},
#'                  \code{"fit failed: ..."}).}
#' }
#'
#' @seealso \code{\link[=nakarushton]{nakarushton}} for the underlying single-channel fit
#'   and \code{\link{naka_rushton4}} for the model function.
#'
#' @importFrom stats aggregate median
#' @importFrom EPhysData TimeTrace Channels nested2df
#' @examples
#' \dontrun{
#' fit_df <- NakaRushton(X)
#' }
NULL


#' @noRd
#' @export
setGeneric("NakaRushton",
           function(X,
                    ...) standardGeneric("NakaRushton"))

#' @rdname NakaRushton
#' @export
setMethod("NakaRushton", "EPhysContinuous",
          function(X,
                   parallel = !interactive(),
                   progress = interactive(),
                   lapply_error = c("stop","warn")[1],
                   ...) {

            if(!HasStimulus(X)){
              stop("Object needs a stimulus trace containing local Michelson Contrast")
            }

            I   <- as.numeric(StimulusTrace(X))   # contrasts (coerced)

            # --- Fit per (run, channel) directly via lapply(EPhysContinuous, ...) ---
            fits_nested <- lapply(
              X        = X,
              FUN      = function(y) nakarushton(I = I, R = y),
              parallel = parallel,
              error    = match.arg(lapply_error, c("stop","warn")),
              progress = progress
            )
            out<-nested2df(X,fits_nested)
            out
          })



#' @rdname NakaRushton
#' @noRd
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats coef predict quantile
nakarushton <- function(I, R) {
  I <- as.numeric(I)
  R <- as.numeric(R)
  ok <- is.finite(I) & is.finite(R)
  R <- R[ok]

  # Clamp to Michelson domain
  I <- pmin(pmax(I[ok], 0), 1)

  if (length(I) < 3L || length(unique(I)) < 3L) {
    return(list(
      coef = c(Baseline = NA, Rmax = NA, k = NA, n = NA),
      rmse = NA_real_, r2 = NA_real_, converged = FALSE,
      message = "insufficient data"
    ))
  }

  # Soft "hard-limit" at I=0 and I=1 via duplicated pseudo-points
  q_lo <- quantile(I, 0.2, na.rm = TRUE)
  q_hi <- quantile(I, 0.9, na.rm = TRUE)
  R_low  <- mean(R[I <= q_lo], na.rm = TRUE)
  R_high <- mean(R[I >= q_hi], na.rm = TRUE)

  w <- rep(1, length(I))
  w[I <= q_lo] <- 5.0
  w[I >= q_hi] <- 10.0

  dat <- data.frame(I = I, R = R)

  # Model and bounds (k ≤ 0.99 since I ∈ [0,1])
  fmla <- R ~ naka_rushton4(I, Baseline, Rmax, k, n)

  bounds <- .guess_nr4(dat$I, dat$R)  # or .guess_nr4(dat$I, dat$R, R_high = R_high)

  start     <- as.numeric(bounds["start", ])
  lower_all <- as.numeric(bounds["lower_all", ])
  upper_all <- as.numeric(bounds["upper_all", ])
  names(start) <- names(lower_all) <- names(upper_all) <- colnames(bounds)

  fit <- tryCatch(
    nlsLM(
      fmla,
      data   = dat,
      start  = start,
      lower  = lower_all[names(start)],
      upper  = upper_all[names(start)],
      weights = w,
      control = nls.lm.control(maxiter = 1000)
    ),
    error = identity
  )

  if (inherits(fit, "error")) {
    return(list(Baseline = NA, Rmax = NA, k = NA, n = NA,
                rmse = NA_real_, r2 = NA_real_, converged = FALSE,
                message = paste("fit failed:", fit$message)))
  }


  cf <- coef(fit)
  if (cf["k"] == 1){
    idx<-which.min(abs(R - (min(R) + (max(R)-min(R))/2)))
    if(length(idx)>1){
      cf["k"]  <- I[round(length(idx)/2)]
    } else {
      cf["k"]  <- I[idx]
    }
    return(
      list(
        Baseline = unname(cf["Baseline"]),
        Rmax     = unname(cf["Rmax"]),
        k        = unname(cf["k"]),
        n        = unname(cf["n"]),
        rmse = NA,
        r2 = NA,
        converged = TRUE,
        message = "k replaced by x at median for cieling to 1."
      )
    )

  }


  fitted <- predict(fit)

  rmse  <- sqrt(mean((dat$R - fitted)^2))
  denom <- sum((dat$R - mean(dat$R))^2)
  r2    <- if (denom > 0) 1 - sum((dat$R - fitted)^2) / denom else NA_real_


  list(
    Baseline = unname(cf["Baseline"]),
    Rmax     = unname(cf["Rmax"]),
    k        = unname(cf["k"]),
    n        = unname(cf["n"]),
    rmse = rmse,
    r2 = r2,
    converged = TRUE,
    message = "ok"
  )
}


#' Four-parameter Naka–Rushton with baseline
#'
#' Internal helper implementing the baseline + amplitude form:
#' \eqn{R(I) = \mathrm{Baseline} + R_{\max} \, I^n / (I^n + k^n)}.
#'
#' @param I Numeric vector of contrasts.
#' @param Baseline Baseline (offset).
#' @param Rmax Amplitude above baseline.
#' @param k Semi-saturation contrast (C50).
#' @param n Slope/steepness.
#' @return Numeric vector of fitted responses (same length as \code{I}).
#' @rdname NakaRushton
#' @keywords internal
naka_rushton4 <- function(I, Baseline, Rmax, k, n) {
  Baseline + Rmax * I ^ n / (I ^ n + k ^ n)
}

#' Starting-value and bounds helper for 4-parameter Naka–Rushton
#'
#' Internal helper that derives reasonable parameter \emph{bounds} and a strictly
#' interior \emph{start} vector for fitting
#' \eqn{R(I) = \mathrm{Baseline} + R_{\max}\, I^n / (I^n + k^n)} with \code{nlsLM()}.
#' It uses the observed \code{(I, R)} pairs to set:
#' \itemize{
#'   \item \strong{Baseline}: lower bound at \code{0}, upper bound at
#'         \code{max(R)/2} or \code{1.5 * R[1]} (whichever is larger).
#'   \item \strong{Rmax}:     lower bound at \code{min(R)}, upper bound at
#'         \code{max(0, R_high * 1.1)}.
#'   \item \strong{k}:        lower bound at \code{max(.Machine$double.eps, min(I))},
#'         upper bound at \code{0.99}.
#'   \item \strong{n}:        lower/upper bounds at \code{1} and \code{20}.
#' }
#' The initial guess is derived from min/max and half-height heuristics and then
#' projected to lie strictly inside the open interval \code{(lower, upper)} by
#' at least \code{eps}.
#'
#' @param I Numeric vector of contrasts (optionally already rescaled to \code{[0, 1]}).
#' @param R Numeric vector of responses.
#' @param R_high Optional numeric scalar used to set the upper bound for \code{Rmax}.
#'   If \code{NULL} (default), \code{max(R, na.rm = TRUE)} is used.
#' @param eps Numeric scalar (default \code{.Machine$double.eps}) used to ensure
#'   \code{start} is strictly interior to the bounds.
#'
#' @return A \code{data.frame} with three rows and four columns:
#'   rows are \code{"lower_all"}, \code{"start"}, \code{"upper_all"};
#'   columns are the parameter names \code{Baseline}, \code{Rmax}, \code{k}, \code{n}.
#'
#' @details
#' Non-finite values in \code{I} or \code{R} are dropped for the calculations.
#' If the dynamic range of \code{R} is degenerate (no spread), conservative
#' defaults are used (\code{n = 1}, \code{Rmax = max(1e-6, max(R))}, \code{k = median(I)}).
#' Bounds are minimally widened if an upper bound would otherwise be \code{\eqn{\le}} its
#' corresponding lower bound.
#'
#' @section Intended use:
#' Pass the returned components to \code{\link[minpack.lm]{nlsLM}}:
#' \preformatted{
#'   bounds <- .guess_nr4(I, R)
#'   start     <- as.numeric(bounds["start", ])
#'   lower_all <- as.numeric(bounds["lower_all", ])
#'   upper_all <- as.numeric(bounds["upper_all", ])
#'   names(start) <- names(lower_all) <- names(upper_all) <- colnames(bounds)
#'   fit <- nlsLM(fmla, data = dat, start = start, lower = lower_all, upper = upper_all)
#' }
#'
#' @seealso \code{\link[minpack.lm]{nlsLM}}
#' @keywords internal
#' @family Naka–Rushton helpers
#' @noRd
.guess_nr4 <- function(I, R, R_high = NULL, eps = .Machine$double.eps) {
  I <- as.numeric(I); R <- as.numeric(R)
  I <- I[is.finite(I)]; R <- R[is.finite(R)]
  if (length(I) == 0L || length(R) == 0L) {
    stop(".guess_nr4: I and R must contain at least one finite value.")
  }

  if (is.null(R_high) || !is.finite(R_high)) {
    R_high <- max(R, na.rm = TRUE)
  }

  # ----- bounds (match your current choices) -----
  Rmin    <- min(R, na.rm = TRUE)
  RmaxObs <- as.vector(quantile(R, 0.9, na.rm = TRUE))
  dR      <- max(1e-6, RmaxObs - Rmin)

  # # jesses params
  lower_all <- c(
    Baseline = 0,
    Rmax     = Rmin,
    k        = min(I),
    n        = 0.9
  )
  upper_all <- c(
    Baseline = Rmin + (dR*0.1), # max(max(R, na.rm = TRUE) / 2, R_first * 1.5),
    Rmax     = 1.1 * RmaxObs, #max(0, R_high * 1.1),
    k        = 1,
    n        = 10
  )
  # ensure upper > lower; if equal/inverted, widen minimally
  for (nm in names(lower_all)) {
    if (!is.finite(lower_all[[nm]])) lower_all[[nm]] <- 0
    if (!is.finite(upper_all[[nm]])) upper_all[[nm]] <- lower_all[[nm]] + 1
    if (upper_all[[nm]] <= lower_all[[nm]]) {
      upper_all[[nm]] <- lower_all[[nm]] + max(1, abs(lower_all[[nm]])) * 1e-6 + 10 * eps
    }
  }

  Baseline0 <- Rmin
  if (!is.finite(dR) || dR <= 0) {
    Rmax0 <- max(1e-6, dR) #max(1e-6, RmaxObs)     # amplitude / Rmax
    k0 <- stats::median(I, na.rm = TRUE)
    n0 <- 1
  } else {
    Rmax0  <- max(1e-6, dR)
    k0  <-   I[which.min(abs(R - (min(R) + (max(R)-min(R))/2)))]
    if (!is.finite(k0)) k0 <- stats::median(I, na.rm = TRUE)
    n0  <- 1
  }

  start <- c(
    Baseline = Baseline0,
    Rmax     = Rmax0,
    k        = k0,
    n        = n0
  )

  # strictly inside (lower, upper)
  # project first, then nudge off the boundaries by eps
  # if (!(all(start >= lower_all))){
  #   warning("start values below range")
  # }
  # if (!(all(start <= upper_all))){
  #   warning("start values above range")
  # }
  start <- pmin(pmax(start, lower_all), upper_all)
  start <- pmin(pmax(start, lower_all + eps), upper_all - eps)

  # In pathological cases where (upper - lower) <= 2*eps, center between them
  tight <- (upper_all - lower_all) <= (2 * eps)
  if (any(tight)) {
    start[tight] <- (lower_all[tight] + upper_all[tight]) / 2
  }

  # return as a tidy data.frame with consistent names
  out <- rbind(lower_all, start, upper_all)
  rownames(out) <- c("lower_all", "start", "upper_all")
  as.data.frame(out, check.names = FALSE)
}


#' @rdname NakaRushton
#' @keywords internal
nr_predict4 <- function(fit, I_new) {
  B <- fit$coef["Baseline"]
  A <- fit$coef["Rmax"]
  k <- fit$coef["k"]
  n <- fit$coef["n"]
  if (any(!is.finite(c(B, A, k, n))))
    return(rep(NA_real_, length(I_new)))
  B + A * I_new ^ n / (I_new ^ n + k ^ n)
}

#' @rdname NakaRushton
#' @keywords internal
plot_nr_gg4 <-
  function(I,
           R,
           fit,
           n_grid = 400,
           log_x = FALSE,
           title = NULL) {
    I <- as.numeric(I)
    R <- as.numeric(R)
    ok <- is.finite(I) & is.finite(R)
    I <- I[ok]
    R <- R[ok]

    rng <- range(I)
    Igrid <- seq(rng[1], rng[2], length.out = n_grid)
    df_points <- data.frame(I = I, R = R)
    df_curve  <- data.frame(I = Igrid, R = nr_predict4(fit, Igrid))

    p <- ggplot2::ggplot(df_points, ggplot2::aes(I, R)) +
      ggplot2::geom_point() +
      ggplot2::geom_line(data = df_curve, ggplot2::aes(I, R)) +
      ggplot2::labs(
        X = "Contrast (I)",
        y = "Response (R)",
        title = title,
        subtitle = if (isTRUE(fit$converged))
          sprintf(
            "B=%.3g, Rmax=%.3g, k=%.3g, n=%.3g  (R²=%.2f)",
            fit$coef["Baseline"],
            fit$coef["Rmax"],
            fit$coef["k"],
            fit$coef["n"],
            fit$r2
          )
        else
          "fit did not converge"
      )

    if(!is.na(fit$coef["k"])){
      p <- p +
        geom_vline(xintercept = fit$coef["k"], linetype = 2) +
        geom_hline(yintercept = nr_predict4(fit, fit$coef["k"]), linetype = 3)
    }

    if (log_x)
      p <- p + ggplot2::scale_x_log10()
    p
  }
