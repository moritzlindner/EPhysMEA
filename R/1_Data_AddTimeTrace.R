#' Add or replace the TimeTrace of an EPhysContainer
#'
#' Fills the \code{TimeTrace} slot using \code{base::seq()} and sets the
#' \code{TimeUnits} slot after validating the unit symbol.
#'
#' @inheritParams base::seq
#' @param x An \code{EPhysContainer}.
#' @param Unit Character. A valid udunits time unit symbol (e.g., \code{"s"}, \code{"ms"}).
#'
#' @return The modified \code{EPhysContainer} object.
#' @details
#' \code{Unit} is checked by attempting \code{units::set_units(1, Unit, mode = "standard")}.
#' The resulting \code{TimeTrace} must be strictly increasing (i.e., \code{diff(TimeTrace) > 0}).
#'
#' @seealso \code{\link[base]{seq}}
#' @examples
#' # Add a 0..1 s trace in 1 ms steps:
#' # obj <- AddTimeTrace(obj, from = 0, to = 1, by = 0.001, Unit = "s")
#' @importClassesFrom EPhysData EPhysContainer
#' @name AddTimeTrace
#' @export
setGeneric(
  name = "AddTimeTrace",
  def  = function(x, from, to, by, Unit) {
    standardGeneric("AddTimeTrace")
  }
)

#' @importFrom units set_units
#' @rdname AddTimeTrace
#' @noRd
setMethod(
  "AddTimeTrace",
  signature(x = "EPhysContainer"),
  function(x, from, to, by, Unit) {

    if (!requireNamespace("units", quietly = TRUE)) {
      stop("Package 'units' is required for AddTimeTrace().")
    }

    ## Validate Unit
    ok <- tryCatch({
      set_units(1, Unit, mode = "standard")
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      stop("`Unit` must be a valid udunits symbol (e.g., 's', 'ms'). Got: ", Unit)
    }

    md<-Metadata(x)
    if ("Diff" %in% colnames(md)){
      if(length(unique(md$Diff))!=1){
        stop("All runs in object must have the same duration (Metadata column Diff)")
      }
    }

    ## Build sequence
    if (missing(from) || missing(to) || missing(by)) {
      stop("`from`, `to`, and `by` are all required.")
    }
    tt <- seq(from = from, to = to, by = by)
    if (!is.numeric(tt)) tt <- as.numeric(tt)
    if (length(tt) < 1L) stop("Sequence is empty; check `from`, `to`, and `by`.")

    ## Enforce strictly increasing (keeps downstream binning safe)
    if (!isTRUE(all(diff(tt) > 0))) {
      stop("TimeTrace must be strictly increasing. Check `from`, `to`, and `by`.")
    }

    ## Write slots
    x@TimeTrace <- tt
    x@TimeUnits <- Unit

    validObject(x)
    x
  }
)
