#' Compute channel order for EPhysContainer
#'
#' Returns the channel name order according to \code{sort_by}.
#'
#' @param X   An \code{EPhysContainer} X.
#' @param sort_by  One of:
#'   \itemize{
#'     \item \code{NULL}: original channel order.
#'     \item integer index vector into \code{Channels(X)}.
#'     \item character vector of channel names (others are appended).
#'     \item single bare column name from \code{ChannelMetadata(X)}.
#'     \item numeric vector of length \code{length(Channels(X))}.
#'   }
#' @return Character vector of channel names in the desired order.
#' @keywords internal
#' @importClassesFrom EPhysData EPhysContainer
#' @importFrom EPhysData Metadata Channels ChannelMetadata
#' @name SortChannels
setGeneric("SortChannels", function(X, sort_by = NULL)
  standardGeneric("SortChannels"))

#' @describeIn SortChannels Method for EPhysContainer
#' @keywords internal
setMethod("SortChannels", signature(X = "EPhysContainer"),
          function(X, sort_by = NULL) {

            ch <- as.character(Channels(X))
            if (length(ch) == 0L) return(character(0))

            # 1) NULL → original order
            if (is.null(sort_by)) return(ch)

            # 3) bare column name from Channel_Metadata → numeric desc
            if (is.character(sort_by) && length(sort_by) == 1L) {
              if (!is.null(ChannelMetadata(X)) && sort_by %in% names(ChannelMetadata(X))) {
                vec <- ChannelMetadata(X,sort_by)
                if (!(is.numeric(vec) || is.logical(vec)) || length(vec) != length(ch) || anyNA(vec)) {
                  stop("Channel_Metadata column '", sort_by, "' must be numeric and alignable to channels.")
                }
                ord_idx <- order(-as.numeric(vec), seq_along(ch)) # stable ties
                return(ch[ord_idx])
              }
            }

            # 4) explicit character vector of channel names
            if (is.character(sort_by) && length(sort_by) >= 1L) {
              unknown <- setdiff(sort_by, ch)
              if (length(unknown)) stop("Unknown channel name(s) in sort_by: ", paste(unknown, collapse = ", "))
              return(c(sort_by, setdiff(ch, sort_by)))
            }

            # 5) integer indices into channels
            if (is.numeric(sort_by) && all(sort_by %% 1 == 0)) {
              idx <- as.integer(sort_by)
              if (any(idx < 1L | idx > length(ch))) stop("sort_by indices out of range.")
              return(c(ch[idx], ch[-idx]))
            }

            # 6) numeric metric per channel, descending
            if (is.numeric(sort_by) && length(sort_by) == length(ch)) {
              ord_idx <- order(-as.numeric(sort_by), seq_along(ch))
              return(ch[ord_idx])
            }

            stop("Unsupported sort_by for SortChannels().")
          }
)


#' Compute metadata row order for EPhysContainer
#'
#' Returns the **row indices** of Metadata(X) in the desired order.
#'
#' @param X  An \code{EPhysContainer} X.
#' @param sort_by One of:
#'   \itemize{
#'     \item \code{NULL}: original row order.
#'     \item integer index vector into rows of \code{Metadata(X)} (others appended).
#'     \item single bare column name from \code{Metadata(X)} (numeric/logical; ordered descending).
#'     \item numeric vector of length \code{nrow(Metadata(X))}: per-row metric (descending).
#'   }
#' @return Integer vector of row indices in the desired order.
#' @keywords internal
#' @name SortRun
setGeneric("SortRun", function(X, sort_by = NULL)
  standardGeneric("SortRun"))

#' @describeIn SortRun Method for EPhysContainer
#' @keywords internal
setMethod("SortRun", signature(X = "EPhysContainer"),
          function(X, sort_by = NULL) {

            md <- tryCatch(Metadata(X), error = function(e) NULL)
            if (is.null(md) || !NROW(md)) return(integer(0))

            n  <- nrow(md)
            ix <- seq_len(n)

            # 1) NULL → original order
            if (is.null(sort_by)) return(ix)

            # 2) bare column name from Metadata → numeric/logical desc
            if (is.character(sort_by) && length(sort_by) == 1L && sort_by %in% names(md)) {
              v <- md[[sort_by]]
              if (!(is.numeric(v) || is.logical(v)))
                stop("Metadata column '", sort_by, "' must be numeric/logical for sorting.")
              if (length(v) != n || anyNA(v))
                stop("Metadata column '", sort_by, "' has wrong length or contains NA.")
              o <- order(-as.numeric(v), ix)  # stable ties by original order
              return(o)
            }

            # 3) numeric metric vector of length n → descending
            if (is.numeric(sort_by) && length(sort_by) == n) {
              o <- order(-as.numeric(sort_by), ix)  # stable ties
              return(o)
            }

            # 4) integer index vector into rows (unspecified appended)
            if (is.numeric(sort_by) && all(sort_by %% 1 == 0)) {
              idx <- as.integer(sort_by)
              if (any(idx < 1L | idx > n)) stop("sort_by row indices out of range.")
              return(c(idx, setdiff(ix, idx)))
            }

            stop("Unsupported sort_by for SortRun().")
          })
