#' Load a Recording's Events and Channel Spikes into an EPhysEvents Object
#'
#' Imports the table of contents (trial/event metadata) and raw channel
#' timestamp data for a single recording, assigns spikes to event windows,
#' and returns a populated \link{EPhysEvents} container.
#'
#' @section Offsets:
#' The \code{start_offset} and \code{end_offset} parameters accept either:
#' \itemize{
#'   \item a single \strong{numeric} value — applied to \emph{all} rows of
#'         \code{Metadata$Start} (or \code{Metadata$Stop}, respectively); or
#'   \item a \strong{named list} with one single numeric value per name — applied
#'         \emph{selectively} to rows where \code{Metadata$Experiment} equals the
#'         corresponding list name. Rows whose \code{Experiment} is not present in the
#'         list receive an offset of \code{0}.
#' }
#' Offsets are interpreted in the same time units as the imported metadata
#' (seconds, by convention in this package).
#'
#' @param file_root Character. Path prefix identifying the recording (used by the
#'   import helpers to locate files).
#' @param type_assignment_file Character. Path to the Experiment type-assignment file used
#'   by \code{Import_ToC()} when constructing the table of contents.
#' @param start_offset Either a single numeric offset applied to every row's
#'   \code{Start} time, or a named list with one single numeric value per name
#'   (event \code{Experiment}), in which case the offset is applied only to matching
#'   rows in \code{Metadata$Experiment}. Default \code{0}.
#' @param end_offset Either a single numeric offset applied to every row's
#'   \code{Stop} time, or a named list with one single numeric value per name
#'   (event \code{Experiment}), in which case the offset is applied only to matching
#'   rows in \code{Metadata$Experiment}. Default \code{0}.
#'
#' @return An \code{EPhysEvents} S4 object with slots:
#' \describe{
#'   \item{\code{Metadata}}{Data frame of event/table-of-contents rows (with
#'         \code{Start}, \code{Stop}, and \code{Diff} updated after applying
#'         offsets).}
#'   \item{\code{Data}}{A list (one element per event row) of named numeric
#'         vectors; each element contains spike timestamps per channel assigned
#'         to that window.}
#'   \item{\code{Channels}}{Character vector of channel names (identical across
#'         all windows).}
#'   \item{\code{TimeUnits}}{Character scalar of the time unit (``s'').}
#'   \item{\code{ExamInfo}, \code{SubjectInfo}}{Lists with exam/subject metadata
#'         (empty by default).}
#'   \item{\code{Imported}}{POSIXct timestamp of the import.}
#' }
#'
#' @details
#' After importing the table of contents via \code{Import_ToC()}, the function
#' applies \code{start_offset} to \code{Metadata$Start} and \code{end_offset} to
#' \code{Metadata$Stop}, then recomputes \code{Metadata$Diff = Stop - Start}.
#' It then imports raw channel timestamps with \code{Import_ChannelData()} and
#' assigns spikes to ToC windows using \code{Assign_Events_To_ToC()}. A
#' consistency check ensures that channel names are identical across all event
#' windows.
#'
#' @examples
#' \dontrun{
#' # 1) Global offsets (apply to all rows)
#' ev1 <- LoadEPhysEvents(
#'   file_root = "path/to/recA",
#'   type_assignment_file = "path/to/types.csv",
#'   SingleTrialTimings = TRUE,
#'   start_offset =  0.050,   # shift all starts by +50 ms
#'   end_offset   = -0.025    # shift all stops  by -25 ms
#' )
#'
#' # 2) Experiment-specific offsets (apply only where Metadata$Experiment matches)
#' ev2 <- LoadEPhysEvents(
#'   file_root = "path/to/recB",
#'   type_assignment_file = "path/to/types.csv",
#'   start_offset = list(Stim = 0.02, Blank = -0.01),
#'   end_offset   = list(Stim = 0.00)  # others get 0 by default
#' )
#' }
#'
#' @seealso \code{\link{EPhysEvents}}, \code{\link{Import_ToC}},
#'   \code{\link{Import_ChannelData}}, \code{\link{Assign_Events_To_ToC}}
#' @keywords IO import electrophysiology spikes events
#' @importClassesFrom EPhysData EPhysEvents
#' @export
LoadEPhysEvents <- function(
    file_root,
    type_assignment_file,
    start_offset       = 0,
    end_offset         = 0
) {
  # --- helper: turn an offset spec into a per-row numeric delta ----------------
  .offset_vector <- function(types, spec, spec_name = "offset") {
    # Case 1: single numeric
    if (is.numeric(spec) && length(spec) == 1L && is.finite(spec)) {
      return(rep(spec, length(types)))
    }
    # Case 2: named list of single numerics
    if (is.list(spec) && length(spec) >= 1L && !is.null(names(spec))) {
      vals_ok <- vapply(spec, function(z) is.numeric(z) && length(z) == 1L && is.finite(z),
                        logical(1))
      if (!all(vals_ok)) {
        stop(sprintf(
          "%s: when given as a list, every element must be a single finite numeric.",
          spec_name
        ))
      }
      # Ensure we can match on Experiment
      if (is.null(types)) stop("Metadata$Experiment is NULL; cannot apply Experiment-specific offsets.")
      if (!is.atomic(types)) stop("Metadata$Experiment must be an atomic vector for matching offsets.")
      out <- numeric(length(types))
      m   <- match(types, names(spec))
      has <- !is.na(m)
      if (any(has)) out[has] <- unlist(spec, use.names = FALSE)[m[has]]
      return(out)
    }
    stop(sprintf(
      "%s must be a single numeric or a named list of single numerics.", spec_name
    ))
  }

  # 1) import table of contents (events)
  Metadata <- Import_ToC(
    file_root,
    type_assignment_file,
    SingleTrialTimings = TRUE
  )

  # Basic sanity checks for required columns
  req_cols <- c("Start", "Stop")
  miss <- setdiff(req_cols, names(Metadata))
  if (length(miss)) {
    stop("Import_ToC() did not return required columns: ", paste(miss, collapse = ", "))
  }
  # Warn if Experiment is missing when list offsets were provided
  if ((is.list(start_offset) || is.list(end_offset)) && !"Experiment" %in% names(Metadata)) {
    stop("Experiment-specific offsets were provided, but Metadata lacks a 'Experiment' column.")
  }

  # 2) apply offsets and compute Diff correctly
  start_delta <- .offset_vector(Metadata$Experiment, start_offset, "start_offset")
  end_delta   <- .offset_vector(Metadata$Experiment, end_offset,   "end_offset")

  Metadata$Start <- Metadata$Start + start_delta
  Metadata$Stop  <- Metadata$Stop  + end_delta
  Metadata$Diff  <- Metadata$Stop  - Metadata$Start

  # Coerce common ID columns (if present)
  to_int <- intersect(c("Repeat", "RecordingID", "Trials", "RunUID"), names(Metadata))
  for (nm in to_int) Metadata[[nm]] <- as.integer(Metadata[[nm]])

  # 3) import raw channel timestamp data
  Data_raw <- Import_ChannelData(file_root)

  # 4) assign events to ToC windows
  Data <- Assign_Events_To_ToC(Metadata, Data_raw)

  # Ensure all event windows carry identical channel name sets
  chan_lists   <- lapply(Data, names)
  first_chan   <- chan_lists[[1L]]
  all_identical<- all(vapply(chan_lists, function(x) identical(x, first_chan), logical(1)))
  if (!all_identical) {
    stop("Channel names are not identical across all event windows.")
  }
  Channels <- first_chan

  # 5) wrap in the new EPhysEvents container
  out <- new(
    "EPhysEvents",
    Metadata    = Metadata,
    Data        = Data,
    Channels    = Channels,
    TimeUnits   = "s",
    ExamInfo    = list(),
    SubjectInfo = list(),
    Imported    = as.POSIXct(Sys.time())
  )
  out
}
