
#' Helper: coerce numeric or units to plain seconds
#'
#' Safely converts inputs that may carry a \pkg{units} class to plain numeric
#' seconds. If \code{v} has units convertible to seconds, they are converted via
#' \code{units::set_units(v, "s")}. If conversion fails (e.g., length units),
#' the function issues a warning and returns the underlying magnitudes after
#' \code{units::drop_units(v)}. For \code{difftime} input, values are returned
#' in seconds. Non-numeric inputs are coerced with \code{as.numeric()} and may
#' yield \code{NA_real_}.
#'
#' @param v A vector. Typically numeric or a \code{units} vector; may also be
#'   \code{NULL}, \code{difftime}, or any object coercible to numeric.
#'
#' @return A numeric vector. If \code{v} is \code{NULL}, \code{NULL} is returned.
#'
#' @details
#' The fallback of dropping non-time units is intentional to keep downstream
#' code from erroring out, but you may wish to check upstream logic if you see
#' a warning about returning raw magnitudes.
#'
#' @keywords internal
#' @rdname step_stim_resp_metrics
#' @importFrom units set_units drop_units deparse_unit
#'
#' @examples
#' if (requireNamespace("units", quietly = TRUE)) {
#'   x <- units::set_units(1:3, "s")
#'   to_num(x)   # 1, 2, 3
#'
#'   y <- units::set_units(1:3, "ms")
#'   to_num(y)   # 0.001, 0.002, 0.003
#'
#'   z <- units::set_units(1:3, "cm")
#'   to_num(z)   # warns; returns 1, 2, 3 (raw magnitudes)
#' }
#' to_num(c(0.1, NA, 2))
#' to_num(as.difftime(c(1, 2), units = "mins"))  # 60, 120
to_num <- function(v) {
  # NULL → NULL (preserve sentinel use)
  if (is.null(v)) return(v)

  # Zero-length vector → zero-length numeric
  if (length(v) == 0L) return(numeric(0))

  # 'units' inputs: try converting to seconds
  if (inherits(v, "units")) {
    if (deparse_unit(v) == ""){
      return(as.numeric(v))
    }
    out <- tryCatch({
      v_s <- set_units(v, "s")   # may error if not time-convertible
      drop_units(v_s)            # strip units → numeric
    }, error = function(e) {
      # Fallback: try dropping units as-is (raw magnitudes), warn the user
      raw <- try(drop_units(v), silent = TRUE)
      unit_lbl <- tryCatch(deparse_unit(v), error = function(...) "<?>" )
      if (!inherits(raw, "try-error")) {
        warning(
          sprintf("to_num(): cannot convert from '%s' to seconds; returning raw magnitudes (units dropped).", unit_lbl),
          call. = FALSE
        )
        as.numeric(raw)
      } else {
        warning(
          sprintf("to_num(): failed to handle 'units' input; returning NA_real_. Reason: %s",
                  conditionMessage(e)),
          call. = FALSE
        )
        rep(NA_real_, length(v))
      }
    })
    return(as.numeric(out))
  }

  # difftime inputs: return seconds
  if (inherits(v, "difftime")) {
    return(as.numeric(v, units = "secs"))
  }

  # Everything else: best-effort numeric coercion (quietly)
  suppressWarnings(as.numeric(v))
}
