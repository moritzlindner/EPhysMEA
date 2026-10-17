# Internal helpers for ImportMCS(): MC_Rack Raster Plot ASCII export.

.mcs_list_raster_files <- function(path) {
  files <- path
  if (length(path) == 1L && dir.exists(path)) {
    files <- list.files(
      path,
      pattern = "\\.(txt|dat|asc|csv|tsv)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) {
    stop("No MC_Rack raster ASCII files found in: ", paste(path, collapse = ", "))
  }
  sort(normalizePath(files, winslash = "/", mustWork = TRUE))
}

.mcs_raster_time_scale <- function(header_text, raster_time_unit) {
  raster_time_unit <- match.arg(raster_time_unit, c("auto", "s", "ms", "us"))
  if (!identical(raster_time_unit, "auto")) {
    return(switch(raster_time_unit, s = 1, ms = 1e-3, us = 1e-6))
  }

  h <- tolower(header_text)
  if (grepl("microsecond|µs|μs|\\bus\\b", h, perl = TRUE)) return(1e-6)
  if (grepl("millisecond|\\bms\\b", h, perl = TRUE)) return(1e-3)
  if (grepl("second|\\bsec\\b|\\[s\\]|\\(s\\)", h, perl = TRUE)) return(1)

  warning(
    "Could not determine the MC_Rack raster time unit from the header; ",
    "assuming seconds. Set 'raster_time_unit' explicitly if necessary.",
    call. = FALSE
  )
  1
}

.mcs_parse_raster_file <- function(file, unit_tags = NULL,
                                   raster_time_unit = c("auto", "s", "ms", "us")) {
  raster_time_unit <- match.arg(raster_time_unit)
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) stop("MC_Rack raster file is empty: ", file)

  # MC_Rack Raster Plot ASCII exports are tab-delimited. Some exports use
  # semicolons/commas depending on locale; parse all three separators here.
  split_line <- function(x) {
    trimws(strsplit(x, "[\\t;,]+", perl = TRUE)[[1L]])
  }

  parts <- lapply(lines, split_line)
  is_data <- vapply(parts, function(z) {
    if (length(z) < 2L) return(FALSE)
    a <- suppressWarnings(as.numeric(sub(",", ".", z[1L], fixed = TRUE)))
    b <- suppressWarnings(as.numeric(sub(",", ".", z[2L], fixed = TRUE)))
    is.finite(a) && is.finite(b)
  }, logical(1))

  if (!any(is_data)) {
    stop(
      "No numeric spike rows could be detected in MC_Rack raster file: ", file,
      ". Expected at least time and sweep columns."
    )
  }

  first_data <- which(is_data)[1L]
  header_text <- paste(lines[seq_len(max(1L, first_data - 1L))], collapse = " ")
  time_scale <- .mcs_raster_time_scale(header_text, raster_time_unit)

  rows <- parts[is_data]
  ncol_max <- max(vapply(rows, length, integer(1)))
  mat <- matrix(NA_character_, nrow = length(rows), ncol = ncol_max)
  for (i in seq_along(rows)) mat[i, seq_along(rows[[i]])] <- rows[[i]]

  num <- apply(mat, 2L, function(z) suppressWarnings(as.numeric(sub(",", ".", z, fixed = TRUE))))
  if (is.null(dim(num))) num <- matrix(num, ncol = 1L)

  time <- as.numeric(num[, 1L]) * time_scale
  sweep <- as.integer(round(num[, 2L]))
  tag <- if (ncol(num) >= 3L) as.integer(round(num[, 3L])) else rep(0L, length(time))

  ok <- is.finite(time) & !is.na(sweep) & !is.na(tag)
  time <- time[ok]
  sweep <- sweep[ok]
  tag <- tag[ok]

  if (!length(time)) stop("No valid spike rows remained after parsing: ", file)
  if (!is.null(unit_tags)) {
    keep <- tag %in% unit_tags
    time <- time[keep]
    sweep <- sweep[keep]
    tag <- tag[keep]
  }

  data.frame(Time = time, Sweep = sweep, UnitTag = tag, stringsAsFactors = FALSE)
}

.mcs_raster_channel_base <- function(file) {
  x <- tools::file_path_sans_ext(basename(file))
  x <- sub("(?i)([_ .-]*(raster|spikes?|events?))+$", "", x, perl = TRUE)
  x <- trimws(x)
  if (!nzchar(x)) x <- tools::file_path_sans_ext(basename(file))
  x
}

.mcs_import_raster_group <- function(path,
                                     unit_tags = NULL,
                                     raster_time_unit = c("auto", "s", "ms", "us")) {
  raster_time_unit <- match.arg(raster_time_unit)
  files <- .mcs_list_raster_files(path)

  parsed <- lapply(
    files,
    .mcs_parse_raster_file,
    unit_tags = unit_tags,
    raster_time_unit = raster_time_unit
  )

  channel_events <- list()
  channel_meta <- list()

  for (i in seq_along(files)) {
    df <- parsed[[i]]
    base <- .mcs_raster_channel_base(files[[i]])
    tags <- sort(unique(df$UnitTag))

    # Preserve Spike Sorter identities. If only one tag is present, retain the
    # electrode/file name for the common unsorted (tag 0) use case.
    for (tag in tags) {
      ch <- if (length(tags) == 1L) base else paste0(base, "_unit", tag)
      if (ch %in% names(channel_events)) {
        stop(
          "Duplicate channel name '", ch, "' inferred from MC_Rack files. ",
          "Rename the files so that each electrode/channel has a unique name."
        )
      }
      sub <- df[df$UnitTag == tag, , drop = FALSE]
      channel_events[[ch]] <- split(sub$Time, sub$Sweep)
      channel_meta[[ch]] <- data.frame(
        SourceFile = basename(files[[i]]),
        Electrode = base,
        UnitTag = as.integer(tag),
        stringsAsFactors = FALSE,
        row.names = ch
      )
    }
  }

  # If unit filtering removed every row from every file, return a useful error.
  if (!length(channel_events)) {
    stop("No spike events remain after applying 'unit_tags'.")
  }

  sweeps <- sort(unique(as.integer(unlist(lapply(channel_events, names), use.names = FALSE))))
  if (!length(sweeps)) stop("No sweep identifiers were found in MC_Rack raster data.")

  channels <- names(channel_events)
  data <- lapply(sweeps, function(sw) {
    setNames(lapply(channels, function(ch) {
      x <- channel_events[[ch]][[as.character(sw)]]
      if (is.null(x)) numeric(0) else sort(as.numeric(x))
    }), channels)
  })

  metadata <- data.frame(
    Sweep = as.integer(sweeps),
    stringsAsFactors = FALSE
  )

  list(
    data = data,
    metadata = metadata,
    channel_metadata = do.call(rbind, channel_meta)
  )
}
