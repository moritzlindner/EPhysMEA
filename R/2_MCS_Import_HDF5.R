# Internal helpers for ImportMCS(): MCS RawData HDF5 timestamp streams.

.mcs_h5_attr <- function(x, name, default = NULL) {
  tryCatch({
    val <- hdf5r::h5attr(x, name)
    if (length(val) == 0L) default else val
  }, error = function(e) default)
}

.mcs_h5_names <- function(group) {
  x <- group$ls(recursive = FALSE)
  if (!nrow(x)) character(0) else as.character(x$name)
}

.mcs_h5_has <- function(group, name) {
  name %in% .mcs_h5_names(group)
}

.mcs_h5_as_char <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.raw(x)) x <- rawToChar(x)
  if (is.list(x) && length(x) == 1L) x <- x[[1L]]
  if (is.null(x) || length(x) == 0L) return("")
  x <- as.character(x)
  x <- x[1L]
  if (is.na(x)) return("")
  trimws(gsub("\\x00", "", x, fixed = TRUE))
}

.mcs_h5_info_df <- function(dataset) {
  x <- dataset$read()

  if (is.data.frame(x)) return(x)

  if (is.matrix(x) || is.array(x)) {
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!is.null(colnames(x))) names(out) <- colnames(x)
    return(out)
  }

  if (is.list(x) && !is.null(names(x))) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }

  stop("Unsupported InfoTimeStamp representation in MCS HDF5 file.")
}

.mcs_h5_col <- function(df, name, default = NA) {
  hit <- which(tolower(names(df)) == tolower(name))
  if (!length(hit)) return(rep(default, nrow(df)))
  df[[hit[1L]]]
}

.mcs_h5_select_groups <- function(parent, prefix, selector = NULL, label_attr = "Label") {
  nms <- .mcs_h5_names(parent)
  nms <- nms[grepl(paste0("^", prefix, "_[0-9]+$"), nms)]
  if (!length(nms)) return(character(0))
  nms <- nms[order(as.integer(sub(paste0("^", prefix, "_"), "", nms)))]

  if (is.null(selector)) return(nms)

  selectors <- selector
  keep <- logical(length(nms))

  for (i in seq_along(nms)) {
    g <- parent[[nms[i]]]
    suffix <- suppressWarnings(as.integer(sub(paste0("^", prefix, "_"), "", nms[i])))
    rid <- suppressWarnings(as.integer(.mcs_h5_attr(g, "RecordingID", NA_integer_)))
    lab <- .mcs_h5_as_char(.mcs_h5_attr(g, label_attr, ""))

    for (s in selectors) {
      if (is.numeric(s)) {
        ss <- as.integer(s)
        if ((!is.na(suffix) && suffix == ss) || (!is.na(rid) && rid == ss)) keep[i] <- TRUE
      } else {
        ss <- as.character(s)
        if (identical(nms[i], ss) || identical(lab, ss) ||
            identical(as.character(suffix), ss) || identical(as.character(rid), ss)) {
          keep[i] <- TRUE
        }
      }
    }
  }

  out <- nms[keep]
  if (!length(out)) {
    stop(
      "No ", prefix, " groups matched selector: ",
      paste(selector, collapse = ", "), ". Available: ", paste(nms, collapse = ", ")
    )
  }
  out
}

.mcs_h5_select_streams <- function(ts_parent, selector = NULL) {
  nms <- .mcs_h5_select_groups(ts_parent, "Stream", NULL)
  if (!length(nms)) return(character(0))

  if (is.null(selector)) {
    subtype <- vapply(nms, function(nm) {
      .mcs_h5_as_char(.mcs_h5_attr(ts_parent[[nm]], "DataSubType", ""))
    }, character(1))
    hit <- tolower(subtype) == "neuralspike"
    if (!any(hit)) return(character(0))
    return(nms[hit])
  }

  keep <- logical(length(nms))
  for (i in seq_along(nms)) {
    g <- ts_parent[[nms[i]]]
    suffix <- suppressWarnings(as.integer(sub("^Stream_", "", nms[i])))
    label <- .mcs_h5_as_char(.mcs_h5_attr(g, "Label", ""))
    subtype <- .mcs_h5_as_char(.mcs_h5_attr(g, "DataSubType", ""))

    for (s in selector) {
      if (is.numeric(s)) {
        if (!is.na(suffix) && suffix == as.integer(s)) keep[i] <- TRUE
      } else {
        ss <- as.character(s)
        if (identical(nms[i], ss) || identical(as.character(suffix), ss) ||
            identical(label, ss) || identical(subtype, ss)) keep[i] <- TRUE
      }
    }
  }

  out <- nms[keep]
  if (!length(out)) {
    stop(
      "No TimeStampStream matched 'hdf5_stream'. Available streams: ",
      paste(nms, collapse = ", ")
    )
  }
  out
}

.mcs_h5_entity_channel_name <- function(info_row, entity_id, stream_name, stream_label) {
  src <- .mcs_h5_as_char(.mcs_h5_col(info_row, "SourceChannelLabels", "")[1L])
  lab <- .mcs_h5_as_char(.mcs_h5_col(info_row, "Label", "")[1L])

  candidate <- if (length(src) && nzchar(src)) src else if (length(lab) && nzchar(lab)) lab else paste0("entity", entity_id)
  candidate <- gsub("[,;]+", "_", candidate)
  candidate <- gsub("\\s+", "_", candidate)
  candidate <- gsub("[^[:alnum:]_.-]", "_", candidate)
  candidate <- gsub("_+", "_", candidate)
  candidate <- sub("^_+", "", sub("_+$", "", candidate))
  if (!nzchar(candidate)) candidate <- paste0("entity", entity_id)

  candidate
}

.mcs_import_hdf5_group <- function(path,
                                   recording_selector = NULL,
                                   stream_selector = NULL) {
  if (length(path) != 1L || !is.character(path) || !file.exists(path) || dir.exists(path)) {
    stop("MCS HDF5 import requires exactly one existing file per input group.")
  }

  con <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(con$close_all(), silent = TRUE), add = TRUE)

  protocol <- .mcs_h5_as_char(.mcs_h5_attr(con, "McsHdf5ProtocolType", ""))
  if (!identical(protocol, "RawData")) {
    stop(
      "File is not an MCS RawData HDF5 file (McsHdf5ProtocolType='",
      protocol, "'): ", path
    )
  }

  if (!.mcs_h5_has(con, "Data")) stop("MCS HDF5 file does not contain /Data: ", path)
  data_group <- con[["Data"]]
  rec_names <- .mcs_h5_select_groups(data_group, "Recording", recording_selector)
  if (!length(rec_names)) stop("No Recording_x groups found in MCS HDF5 file: ", path)

  data_out <- list()
  md_out <- list()
  ch_meta_acc <- list()
  out_i <- 0L

  for (rec_name in rec_names) {
    rec <- data_group[[rec_name]]
    if (!.mcs_h5_has(rec, "TimeStampStream")) next

    ts_parent <- rec[["TimeStampStream"]]
    stream_names <- .mcs_h5_select_streams(ts_parent, stream_selector)
    if (!length(stream_names)) next

    row_channels <- list()
    row_ch_meta <- list()

    for (stream_name in stream_names) {
      stream <- ts_parent[[stream_name]]
      if (!.mcs_h5_has(stream, "InfoTimeStamp")) {
        warning("Skipping ", rec_name, "/", stream_name, ": missing InfoTimeStamp.", call. = FALSE)
        next
      }

      info <- .mcs_h5_info_df(stream[["InfoTimeStamp"]])
      if (!nrow(info)) next

      ids <- suppressWarnings(as.integer(.mcs_h5_col(info, "TimeStampEntityID", seq_len(nrow(info)) - 1L)))
      stream_label <- .mcs_h5_as_char(.mcs_h5_attr(stream, "Label", stream_name))
      subtype <- .mcs_h5_as_char(.mcs_h5_attr(stream, "DataSubType", ""))

      for (k in seq_len(nrow(info))) {
        entity_id <- ids[k]
        if (is.na(entity_id)) entity_id <- k - 1L
        ds_name <- paste0("TimeStampEntity_", entity_id)
        if (!.mcs_h5_has(stream, ds_name)) {
          warning("Skipping missing dataset ", rec_name, "/", stream_name, "/", ds_name, ".", call. = FALSE)
          next
        }

        ts <- as.numeric(stream[[ds_name]]$read()) / 1e6
        ts <- sort(ts[is.finite(ts)])

        base_name <- .mcs_h5_entity_channel_name(info[k, , drop = FALSE], entity_id, stream_name, stream_label)
        ch <- base_name
        if (ch %in% names(row_channels)) {
          ch <- paste0(base_name, "_", stream_name, "_", entity_id)
        }
        if (ch %in% names(row_channels)) ch <- make.unique(c(names(row_channels), ch))[length(row_channels) + 1L]

        row_channels[[ch]] <- ts

        source_ids <- .mcs_h5_as_char(.mcs_h5_col(info[k, , drop = FALSE], "SourceChannelIDs", "")[1L])
        source_labels <- .mcs_h5_as_char(.mcs_h5_col(info[k, , drop = FALSE], "SourceChannelLabels", "")[1L])
        entity_label <- .mcs_h5_as_char(.mcs_h5_col(info[k, , drop = FALSE], "Label", "")[1L])

        row_ch_meta[[ch]] <- data.frame(
          SourceFile = basename(path),
          Stream = stream_name,
          StreamLabel = stream_label,
          DataSubType = subtype,
          EntityID = as.integer(entity_id),
          EntityLabel = entity_label,
          SourceChannelIDs = source_ids,
          SourceChannelLabels = source_labels,
          stringsAsFactors = FALSE,
          row.names = ch
        )
      }
    }

    if (!length(row_channels)) next

    out_i <- out_i + 1L
    data_out[[out_i]] <- row_channels

    rec_id <- suppressWarnings(as.integer(.mcs_h5_attr(rec, "RecordingID", sub("^Recording_", "", rec_name))))
    rec_duration_us <- suppressWarnings(as.numeric(.mcs_h5_attr(rec, "Duration", NA_real_)))
    rec_timestamp_us <- suppressWarnings(as.numeric(.mcs_h5_attr(rec, "TimeStamp", NA_real_)))

    md_out[[out_i]] <- data.frame(
      MCSRecording = rec_name,
      MCSRecordingID = rec_id,
      MCSLabel = .mcs_h5_as_char(.mcs_h5_attr(rec, "Label", "")),
      MCSRecordingType = .mcs_h5_as_char(.mcs_h5_attr(rec, "RecordingType", "")),
      MCSRecordingStart = if (is.finite(rec_timestamp_us)) rec_timestamp_us / 1e6 else NA_real_,
      MCSDuration = if (is.finite(rec_duration_us)) rec_duration_us / 1e6 else NA_real_,
      stringsAsFactors = FALSE
    )

    ch_meta_acc <- c(ch_meta_acc, row_ch_meta)
  }

  if (!length(data_out)) {
    stop("No neural spike timestamp entities were imported from MCS HDF5 file: ", path)
  }

  # Collapse repeated channel metadata across MCS recordings while retaining
  # all distinct values encountered for a channel.
  all_channels <- unique(unlist(lapply(data_out, names), use.names = FALSE))
  channel_metadata <- .mcs_combine_channel_metadata(ch_meta_acc, all_channels)

  list(
    data = data_out,
    metadata = do.call(rbind, md_out),
    channel_metadata = channel_metadata
  )
}
