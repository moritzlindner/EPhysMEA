#' Import spike-event data exported by Multi Channel Systems software
#'
#' Imports spike timestamps produced by MC_Rack / Multi Channel Systems (MCS)
#' into an \code{\link[EPhysData:EPhysEvents-class]{EPhysEvents}} object that can
#' be analysed directly with EPhysMEA.
#'
#' Two input formats are supported:
#' \itemize{
#'   \item \strong{MC_Rack raster ASCII}: the tab-delimited files produced by
#'     exporting a Raster Plot. The first three columns are interpreted as time,
#'     sweep number and Spike Sorter unit tag. One file per electrode/channel is
#'     expected. This format preserves the sweep/trial structure and is therefore
#'     the recommended format for repeated-stimulus experiments.
#'   \item \strong{MCS HDF5}: timestamp streams in the MCS RawData HDF5 format,
#'     in particular streams with \code{DataSubType = "NeuralSpike"}. HDF5
#'     timestamp streams do not themselves encode MC_Rack raster sweeps, so one
#'     imported run is created per selected MCS recording.
#' }
#'
#' A named list can be supplied as \code{path} to import several experimental
#' steps (for example several light intensities) into one EPhysEvents object.
#' The list names are used as \code{Metadata$Step} unless \code{step} is supplied.
#'
#' @param path Character path, character vector of MC_Rack raster files, or a
#'   named list of such inputs. A directory is searched for spike raster files.
#'   For multiple stimulus conditions, use a named list, e.g.
#'   \code{list("1" = "low/", "2" = "medium/", "3" = "high/")}.
#' @param format Input format. \code{"auto"} detects MCS HDF5 files by extension
#'   and otherwise uses MC_Rack raster ASCII. Alternatively use
#'   \code{"mcrack_raster"} or \code{"mcs_hdf5"} explicitly.
#' @param experiment Character scalar stored in \code{Metadata$Experiment}.
#' @param step Optional step label(s). If \code{path} is a named list and
#'   \code{step = NULL}, list names are used. Otherwise steps default to
#'   \code{"1"}, \code{"2"}, ... in input order.
#' @param recording_id Integer scalar or one integer per input group. A scalar is
#'   recycled. Use the same value across intensities when they belong to the
#'   same biological recording.
#' @param time_range Optional numeric \code{c(start, stop)} giving the complete
#'   trial/recording time range in seconds. If omitted, \code{TimeTrace} is left
#'   empty because spike timestamps alone do not define the trial boundaries.
#' @param unit_tags Optional integer vector of MC_Rack Spike Sorter unit tags to
#'   retain. \code{NULL} (default) keeps all tags. For unsorted MC_Rack spikes
#'   the tag is normally 0.
#' @param raster_time_unit Time unit used in MC_Rack raster ASCII files:
#'   \code{"auto"} (default), \code{"s"}, \code{"ms"}, or \code{"us"}.
#'   Automatic mode inspects the file header and falls back to seconds with a
#'   warning if no unit can be identified.
#' @param hdf5_recording Optional selector for MCS HDF5 recordings. May be a
#'   numeric index/RecordingID or a character \code{"Recording_x"} name. By
#'   default all recordings containing matching timestamp streams are imported.
#' @param hdf5_stream Optional selector for MCS HDF5 timestamp streams. May be a
#'   numeric index or a character stream name/label. By default all streams with
#'   \code{DataSubType = "NeuralSpike"} are used.
#'
#' @return An \code{EPhysEvents} object. Times are converted to seconds. Each
#'   raster sweep (or selected HDF5 recording) becomes one metadata row/run and
#'   each electrode becomes a channel.
#'
#' @details
#' MC_Rack Raster Plot ASCII export contains a header row followed by a units row;
#' data rows contain time, sweep and (when Spike Sorter information is present)
#' unit tag. \code{ImportMCS()} combines all channel files by sweep and stores the
#' original sweep number in \code{Metadata$Sweep}.
#'
#' MCS RawData HDF5 timestamp values are specified in microseconds and are
#' converted to seconds. The importer reads \code{TimeStampStream/Stream_x},
#' \code{InfoTimeStamp}, and the corresponding \code{TimeStampEntity_x}
#' datasets. Channel names preferentially use \code{SourceChannelLabels}, with
#' fallbacks to entity labels and IDs.
#'
#' @examples
#' \dontrun{
#' # One light intensity: directory containing MC_Rack raster exports
#' x <- ImportMCS("data/intensity_1", time_range = c(-2, 6))
#'
#' # Several intensities in one object; list names become ordered Step levels
#' x <- ImportMCS(
#'   path = list(
#'     "0.01" = "data/I_001",
#'     "0.1"  = "data/I_01",
#'     "1"    = "data/I_1"
#'   ),
#'   experiment = "retinal_step",
#'   recording_id = 1L,
#'   time_range = c(-2, 6)
#' )
#'
#' # Direct MCS DataManager HDF5 timestamp import
#' x_h5 <- ImportMCS("recording.h5", format = "mcs_hdf5")
#' }
#'
#' @seealso \code{\link{ggEventRaster}}, \code{\link{StepStimResponseProfile}}
#' @keywords IO import electrophysiology MCS MC_Rack spikes events
#' @export
ImportMCS <- function(path,
                      format = c("auto", "mcrack_raster", "mcs_hdf5"),
                      experiment = "MCS",
                      step = NULL,
                      recording_id = 1L,
                      time_range = NULL,
                      unit_tags = NULL,
                      raster_time_unit = c("auto", "s", "ms", "us"),
                      hdf5_recording = NULL,
                      hdf5_stream = NULL) {
  format <- match.arg(format)
  raster_time_unit <- match.arg(raster_time_unit)

  if (!is.character(experiment) || length(experiment) != 1L || is.na(experiment) || !nzchar(experiment)) {
    stop("'experiment' must be one non-empty character value.")
  }

  groups <- .mcs_normalize_groups(path, step)
  n_groups <- length(groups$paths)

  if (!(length(recording_id) %in% c(1L, n_groups))) {
    stop("'recording_id' must have length 1 or one value per input group.")
  }
  recording_id <- rep(recording_id, length.out = n_groups)
  if (!is.numeric(recording_id) || any(!is.finite(recording_id))) {
    stop("'recording_id' must contain finite numeric/integer values.")
  }
  recording_id <- as.integer(recording_id)

  if (!is.null(time_range)) {
    if (!is.numeric(time_range) || length(time_range) != 2L ||
        any(!is.finite(time_range)) || time_range[1] >= time_range[2]) {
      stop("'time_range' must be a finite numeric c(start, stop) with start < stop.")
    }
    time_range <- as.numeric(time_range)
  }

  if (!is.null(unit_tags)) {
    if (!is.numeric(unit_tags) || any(!is.finite(unit_tags))) {
      stop("'unit_tags' must be NULL or a finite numeric/integer vector.")
    }
    unit_tags <- unique(as.integer(unit_tags))
  }

  parsed <- vector("list", n_groups)
  for (i in seq_len(n_groups)) {
    fmt_i <- .mcs_resolve_format(groups$paths[[i]], format)
    parsed[[i]] <- switch(
      fmt_i,
      mcrack_raster = .mcs_import_raster_group(
        groups$paths[[i]],
        unit_tags = unit_tags,
        raster_time_unit = raster_time_unit
      ),
      mcs_hdf5 = .mcs_import_hdf5_group(
        groups$paths[[i]],
        recording_selector = hdf5_recording,
        stream_selector = hdf5_stream
      ),
      stop("Unsupported MCS import format: ", fmt_i)
    )
    parsed[[i]]$format <- fmt_i
  }

  .mcs_build_events(
    parsed = parsed,
    steps = groups$steps,
    experiment = experiment,
    recording_id = recording_id,
    time_range = time_range
  )
}

.mcs_normalize_groups <- function(path, step) {
  if (is.list(path)) {
    if (!length(path)) stop("'path' must not be an empty list.")
    ok <- vapply(path, function(z) is.character(z) && length(z) >= 1L && all(!is.na(z) & nzchar(z)), logical(1))
    if (!all(ok)) stop("Every element of list 'path' must be a non-empty character path/vector.")
    paths <- unname(path)
    path_names <- names(path)
  } else {
    if (!is.character(path) || length(path) < 1L || any(is.na(path) | !nzchar(path))) {
      stop("'path' must be a non-empty character path/vector or a list of such inputs.")
    }
    paths <- list(path)
    path_names <- NULL
  }

  n <- length(paths)
  if (is.null(step)) {
    if (!is.null(path_names) && length(path_names) == n && all(nzchar(path_names))) {
      steps <- path_names
    } else {
      steps <- as.character(seq_len(n))
    }
  } else {
    if (!(length(step) %in% c(1L, n))) {
      stop("'step' must have length 1 or one value per input group.")
    }
    steps <- rep(as.character(step), length.out = n)
    if (any(is.na(steps) | !nzchar(steps))) stop("'step' values must be non-empty.")
  }

  list(paths = paths, steps = steps)
}

.mcs_resolve_format <- function(path, format) {
  if (!identical(format, "auto")) return(format)

  files <- path
  if (length(path) == 1L && dir.exists(path)) {
    files <- list.files(path, full.names = TRUE, recursive = FALSE)
  }
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) {
    stop("No input files found for: ", paste(path, collapse = ", "))
  }

  is_h5 <- grepl("\\.(h5|hdf5)$", files, ignore.case = TRUE)
  if (all(is_h5) && length(files) == 1L) return("mcs_hdf5")
  if (any(is_h5) && !all(is_h5)) {
    stop("Automatic format detection found both HDF5 and non-HDF5 files. Specify 'format' explicitly.")
  }
  if (length(files) > 1L && all(is_h5)) {
    stop("One MCS HDF5 file is accepted per input group. Use a list of paths to import multiple files.")
  }
  "mcrack_raster"
}

.mcs_combine_channel_metadata <- function(rows, channels) {
  if (!length(channels)) return(data.frame())
  if (!length(rows)) return(data.frame(row.names = channels))

  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- data.frame(row.names = channels)

  collapse_vals <- function(v) {
    v <- as.character(v)
    v <- unique(v[!is.na(v) & nzchar(v)])
    if (!length(v)) NA_character_ else paste(v, collapse = ";")
  }

  for (col in cols) {
    out[[col]] <- vapply(channels, function(ch) {
      vals <- unlist(lapply(rows, function(z) {
        if (is.data.frame(z) && ch %in% rownames(z) && col %in% names(z)) z[ch, col] else NA
      }), use.names = FALSE)
      collapse_vals(vals)
    }, character(1))
  }
  out
}

.mcs_build_events <- function(parsed, steps, experiment, recording_id, time_range) {
  channels <- unique(unlist(lapply(parsed, function(z) unlist(lapply(z$data, names), use.names = FALSE)), use.names = FALSE))
  if (!length(channels)) stop("No channels were imported.")

  data_all <- list()
  md_all <- list()
  ch_meta_rows <- list()
  row_counter <- 0L

  for (g in seq_along(parsed)) {
    z <- parsed[[g]]
    if (length(z$data) != nrow(z$metadata)) {
      stop("Internal ImportMCS error: data/metadata row mismatch in group ", g, ".")
    }

    for (j in seq_along(z$data)) {
      row_counter <- row_counter + 1L
      src <- z$data[[j]]
      row <- setNames(lapply(channels, function(ch) {
        if (ch %in% names(src)) sort(as.numeric(src[[ch]])) else numeric(0)
      }), channels)

      if (!is.null(time_range)) {
        bad <- vapply(row, function(v) any(v < time_range[1] | v > time_range[2]), logical(1))
        if (any(bad)) {
          stop(
            "Spike timestamp(s) fall outside 'time_range' in group ", g,
            ", run ", j, ", channel(s): ", paste(names(row)[bad], collapse = ", "), "."
          )
        }
      }
      data_all[[row_counter]] <- row
    }

    md <- z$metadata
    md$Experiment <- factor(rep(experiment, nrow(md)), levels = experiment)
    md$Step <- factor(rep(steps[g], nrow(md)), levels = unique(steps), ordered = TRUE)
    md$RecordingID <- rep(as.integer(recording_id[g]), nrow(md))
    md$SourceFormat <- z$format
    if (!is.null(time_range)) md$Duration <- diff(time_range)
    md_all[[g]] <- md

    cm <- z$channel_metadata
    if (is.data.frame(cm) && nrow(cm)) ch_meta_rows <- c(ch_meta_rows, split(cm, seq_len(nrow(cm))))
  }

  metadata <- do.call(rbind, md_all)
  rownames(metadata) <- NULL
  metadata$RunUID <- as.integer(seq_len(nrow(metadata)))

  first <- c("Experiment", "Step", "RecordingID", "RunUID")
  metadata <- metadata[c(first, setdiff(names(metadata), first))]

  metadata$Experiment <- factor(as.character(metadata$Experiment), levels = experiment)
  metadata$Step <- factor(as.character(metadata$Step), levels = unique(steps), ordered = TRUE)
  metadata$RecordingID <- as.integer(metadata$RecordingID)
  metadata$RunUID <- as.integer(metadata$RunUID)

  channel_metadata <- .mcs_combine_channel_metadata(ch_meta_rows, channels)
  if (nrow(channel_metadata) == 0L) channel_metadata <- data.frame(row.names = channels)
  channel_metadata <- channel_metadata[channels, , drop = FALSE]

  time_trace <- if (is.null(time_range)) numeric(0) else as.numeric(time_range)

  out <- EPhysData::newEPhysEvents(
    Data = data_all,
    TimeTrace = time_trace,
    Metadata = metadata,
    Channels = channels,
    Channel_Metadata = channel_metadata,
    StimulusTrace = numeric(0),
    TimeUnits = "s",
    StimulusUnits = "",
    Imported = as.POSIXct(Sys.time()),
    ExamInfo = list(Source = "Multi Channel Systems", Importer = "EPhysMEA::ImportMCS"),
    SubjectInfo = list()
  )
  methods::validObject(out)
  out
}
