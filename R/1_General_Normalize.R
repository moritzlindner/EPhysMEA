#' Min–max normalize each (run, channel) trace over time
#'
#' Scales the \code{Data} slot of an \code{EPhysContinuous} object so that for every
#' run × channel pair, the time-series is transformed to the range \eqn{[0, 1]}
#' using \code{(x - min) / (max - min)} computed along the time dimension.
#'
#' The input array is assumed to be \code{[time × run × channel]}. The method
#' preserves dimensions and \code{dimnames}. Flat traces (non-finite range or
#' \code{max == min}) are set to all zeros.
#'
#' @param X An \linkS4class{EPhysContinuous} object.
#' @param na.rm Logical; if \code{TRUE} (default), ignore \code{NA}s when computing
#'   per-trace minima and maxima.
#'
#' @return The same \code{EPhysContinuous} object with \code{Data} normalized per
#'   (run, channel) over time.
#'
#' @examples
#' # X is an EPhysContinuous with Data [time × run × channel]
#' Xn <- Normalize(X)
#'
#' @importClassesFrom EPhysData EPhysContinuous
#' @export
setGeneric("Normalize", function(X, ...) standardGeneric("Normalize"))

#' @rdname Normalize
#' @export
setMethod("Normalize", "EPhysContinuous", function(X, na.rm = TRUE) {
  arr <- X@Data
  stopifnot(is.array(arr), length(dim(arr)) == 3L)

  # apply over (run, channel): each x is a length = time vector
  out <- apply(arr, c(2, 3), function(x) {
    r <- range(x, na.rm = na.rm)
    d <- r[2] - r[1]
    if (!is.finite(d) || d <= 0) {
      rep(0, length(x))            # flat or all-NA trace -> zeros
    } else {
      (x - r[1]) / d
    }
  })

  # Ensure original shape and names [time × run × channel]
  dim(out)      <- dim(arr)
  dimnames(out) <- dimnames(arr)

  X@Data <- out
  X
})
