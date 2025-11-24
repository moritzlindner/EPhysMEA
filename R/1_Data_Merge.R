## -------------------------------------------------------------------------
## Low-level helper on EPhysContainer: common checks & metadata merging
## -------------------------------------------------------------------------

.merge_EPhysContainer_common <- function(objects,
                                         match_by = c("Step", "Experiment", "Run"),
                                         prefixes = NULL) {
  if (!is.list(objects) || length(objects) == 0L) {
    stop("'objects' must be a non-empty list.")
  }
  warning(".merge_EPhysContainer_common is untested")
  ## Drop NULL entries defensively
  objects <- Filter(Negate(is.null), objects)
  if (length(objects) == 0L) {
    stop("'objects' contains only NULL entries.")
  }

  ## All objects must inherit from EPhysContainer and share the same concrete class
  if (!all(vapply(objects, function(x) methods::is(x, "EPhysContainer"), logical(1)))) {
    stop("All elements of 'objects' must inherit from 'EPhysContainer'.")
  }
  classes    <- vapply(objects, function(x) class(x)[1L], character(1))
  base_class <- classes[1L]
  if (!all(classes == base_class)) {
    stop("All elements of 'objects' must have the same concrete class.")
  }

  ## Helper: build matching key from Metadata
  match_by   <- unique(as.character(match_by))
  key_from_md <- function(md) {
    if (!is.data.frame(md)) stop("Metadata() must return a data.frame.")
    missing <- setdiff(match_by, names(md))
    if (length(missing)) {
      stop("Required 'match_by' columns missing from Metadata: ",
           paste(missing, collapse = ", "))
    }
    do.call(paste, c(lapply(md[match_by], as.character), sep = "\r"))
  }

  md_ref <- Metadata(objects[[1L]])
  key_ref <- key_from_md(md_ref)
  if (anyDuplicated(key_ref)) {
    stop("Rows in Metadata(objects[[1]]) are not unique with respect to 'match_by'.")
  }
  n_runs <- nrow(md_ref)

  ## Align all other objects to the first by match_by columns
  idx_list <- vector("list", length(objects))
  idx_list[[1L]] <- seq_len(n_runs)

  for (j in seq_along(objects)[-1L]) {
    md_j  <- Metadata(objects[[j]])
    key_j <- key_from_md(md_j)
    if (anyDuplicated(key_j)) {
      stop("Rows in Metadata of object ", j, " are not unique with respect to 'match_by'.")
    }
    idx <- match(key_ref, key_j)
    if (anyNA(idx)) {
      missing_keys <- unique(key_ref[is.na(idx)])
      stop("Some runs in object 1 are not present in object ", j,
           " when matching by ", paste(match_by, collapse = ", "), ":\n",
           paste(missing_keys, collapse = ", "))
    }
    idx_list[[j]] <- idx
  }

  ## Common invariants: TimeTrace, TimeUnits, StimulusTrace, StimulusUnits
  tt_ref   <- TimeTrace(objects[[1L]])
  tu_ref   <- TimeUnits(objects[[1L]])
  stim_ref <- StimulusTrace(objects[[1L]])
  su_ref   <- StimulusUnits(objects[[1L]])

  for (j in seq_along(objects)[-1L]) {
    if (!identical(TimeUnits(objects[[j]]), tu_ref)) {
      stop("TimeUnits differ between objects 1 and ", j, ".")
    }
    if (!isTRUE(all.equal(TimeTrace(objects[[j]]), tt_ref))) {
      stop("TimeTrace differs between objects 1 and ", j, ".")
    }
    if (!identical(StimulusUnits(objects[[j]]), su_ref)) {
      stop("StimulusUnits differ between objects 1 and ", j, ".")
    }
    if (!isTRUE(all.equal(StimulusTrace(objects[[j]]), stim_ref))) {
      stop("StimulusTrace differs between objects 1 and ", j, ".")
    }
  }

  ## Channel names and prefixes ----------------------------------------------
  ch_list <- lapply(objects, Channels)
  n_ch    <- vapply(ch_list, length, integer(1))

  ## Prefixes: use list names when present, otherwise "set1", "set2", ...
  if (is.null(prefixes)) {
    prefixes <- names(objects)
    if (is.null(prefixes)) prefixes <- rep("", length(objects))
    empty <- which(!nzchar(prefixes))
    if (length(empty)) {
      prefixes[empty] <- paste0("set", seq_along(empty))
    }
  } else {
    if (length(prefixes) != length(objects)) {
      stop("'prefixes' must have the same length as 'objects'.")
    }
  }

  new_ch_list <- mapply(function(pref, ch) {
    if (!length(ch)) return(character(0L))
    if (nzchar(pref)) paste(pref, ch, sep = "_") else ch
  }, pref = prefixes, ch = ch_list, SIMPLIFY = FALSE)

  new_channels <- unlist(new_ch_list, use.names = FALSE)

  ## For each object: map old channel -> new channel
  rename_map_list <- Map(function(old, new) {
    if (!length(old)) return(setNames(character(0), character(0)))
    setNames(new, old)
  }, old = ch_list, new = new_ch_list)

  ## Global positions of each object's channels in the merged vector
  offsets  <- c(0L, cumsum(n_ch))
  pos_list <- lapply(seq_along(n_ch), function(j) {
    if (n_ch[j] == 0L) integer(0L) else seq.int(offsets[j] + 1L, offsets[j + 1L])
  })

  ## Align Channel_Metadata rows to channel vector ----------------------------
  align_chmeta_to_channels <- function(chm, channels) {
    if (!is.data.frame(chm) || nrow(chm) == 0L) {
      return(data.frame(row.names = channels))
    }
    rn <- rownames(chm)
    if (!is.null(rn) && all(channels %in% rn)) {
      chm <- chm[channels, , drop = FALSE]
      rownames(chm) <- channels
      return(chm)
    }
    if ("Channel" %in% names(chm)) {
      if (!all(channels %in% chm$Channel)) {
        stop("Channel_Metadata is missing rows for some channels.")
      }
      chm <- chm[match(channels, chm$Channel), , drop = FALSE]
      rownames(chm) <- channels
      return(chm)
    }
    if (nrow(chm) == length(channels)) {
      rownames(chm) <- channels
      return(chm)
    }
    stop("Unable to align Channel_Metadata to channels (incomplete or mismatched rows).")
  }

  ## Build merged Channel_Metadata (existing columns only) --------------------
  chm_aligned_list <- mapply(function(obj, old_ch, new_ch) {
    chm <- obj@Channel_Metadata
    chm <- align_chmeta_to_channels(chm, old_ch)
    if (nrow(chm)) {
      rownames(chm) <- new_ch
      if ("Channel" %in% names(chm)) chm$Channel <- new_ch
    } else {
      chm <- data.frame(row.names = new_ch)
    }
    chm
  }, obj = objects, old_ch = ch_list, new_ch = new_ch_list, SIMPLIFY = FALSE)

  all_cols <- unique(unlist(lapply(chm_aligned_list, names)))
  chm_completed <- lapply(chm_aligned_list, function(df) {
    if (!length(all_cols)) return(df)
    missing <- setdiff(all_cols, names(df))
    for (m in missing) df[[m]] <- NA
    df[all_cols]
  })
  chm_merged <- do.call(rbind, chm_completed)
  if (nrow(chm_merged)) {
    chm_merged <- chm_merged[new_channels, , drop = FALSE]
  } else {
    chm_merged <- data.frame(row.names = new_channels)
  }

  ## SubjectInfo / ExamInfo: diverging fields -> Channel_Metadata -------------
  `%||%` <- function(a, b) if (is.null(a)) b else a

  add_info_diff_cols_multi <- function(chm, info_list, prefix, pos_list) {
    info_list <- lapply(info_list, function(x) x %||% list())
    keys <- unique(unlist(lapply(info_list, names)))
    if (!length(keys)) return(chm)

    n_total <- nrow(chm)
    n_obj   <- length(info_list)

    for (nm in keys) {
      vals <- lapply(info_list, function(info) info[[nm]])

      ## Skip if all values are identical (including all NULL)
      first <- vals[[1L]]
      all_identical <- TRUE
      if (length(vals) > 1L) {
        for (j in 2:length(vals)) {
          if (!isTRUE(identical(vals[[j]], first))) {
            all_identical <- FALSE
            break
          }
        }
      }
      if (all_identical) next

      colname <- paste0(prefix, ".", nm)
      if (colname %in% names(chm)) {
        colname <- make.unique(c(names(chm), colname))[ncol(chm) + 1L]
      }

      vec <- rep(NA, n_total)
      for (j in seq_len(n_obj)) {
        vj  <- vals[[j]]
        pos <- pos_list[[j]]
        if (!length(pos)) next
        if (is.null(vj) || length(vj) == 0L) next
        if (length(vj) == 1L) {
          vec[pos] <- vj
        } else if (length(vj) == length(pos)) {
          vec[pos] <- vj
        } else {
          ## Fallback: recycle first element
          vec[pos] <- vj[1L]
        }
      }

      any_char <- any(vapply(vals, function(v) is.character(v) || is.factor(v), logical(1)))
      if (any_char) {
        chm[[colname]] <- factor(vec)
      } else {
        chm[[colname]] <- vec
      }
    }
    chm
  }

  subj_list <- lapply(objects, function(o) o@SubjectInfo)
  exam_list <- lapply(objects, function(o) o@ExamInfo)

  chm_merged <- add_info_diff_cols_multi(chm_merged, subj_list, "SubjectInfo", pos_list)
  chm_merged <- add_info_diff_cols_multi(chm_merged, exam_list, "ExamInfo",  pos_list)

  ## Build merged SubjectInfo / ExamInfo with only non-diverging fields -------
  merge_info_non_diverging <- function(info_list) {
    info_list <- lapply(info_list, function(x) x %||% list())
    keys <- unique(unlist(lapply(info_list, names)))
    out  <- vector("list", length(keys))
    names(out) <- keys
    keep <- logical(length(keys))
    for (i in seq_along(keys)) {
      nm   <- keys[[i]]
      vals <- lapply(info_list, function(info) info[[nm]])
      first <- vals[[1L]]
      all_identical <- TRUE
      if (length(vals) > 1L) {
        for (j in 2:length(vals)) {
          if (!isTRUE(identical(vals[[j]], first))) {
            all_identical <- FALSE
            break
          }
        }
      }
      if (all_identical) {
        out[[i]] <- first
        keep[i]  <- TRUE
      }
    }
    out[keep]
  }

  subj_merged <- merge_info_non_diverging(subj_list)
  exam_merged <- merge_info_non_diverging(exam_list)

  imported <- objects[[1L]]@Imported

  list(
    base_class         = base_class,
    match_by           = match_by,
    metadata           = md_ref,
    idx_list           = idx_list,
    channels           = new_channels,
    old_channels       = ch_list,
    new_channels_by_obj = new_ch_list,
    rename_map_list    = rename_map_list,
    positions_list     = pos_list,
    channel_metadata   = chm_merged,
    subject_info       = subj_merged,
    exam_info          = exam_merged,
    time_trace         = tt_ref,
    time_units         = tu_ref,
    stimulus_trace     = stim_ref,
    stimulus_units     = su_ref,
    imported           = imported
  )
}

## -------------------------------------------------------------------------
## EPhysContinuous: merge along channel dimension
## -------------------------------------------------------------------------

#' Merge multiple EPhysContinuous objects along channels
#'
#' @param objects  List of EPhysContinuous objects.
#' @param match_by Columns in Metadata used to align runs.
#' @param prefixes Optional character vector of prefixes for channel names
#'   (one per object). If NULL, list names / "set1", "set2", … are used.
#' @export
MergeContinuous <- function(objects,
                            match_by = c("Step", "Experiment", "Run"),
                            prefixes = NULL) {
  if (!is.list(objects) || length(objects) == 0L) {
    stop("'objects' must be a non-empty list.")
  }
  if (!all(vapply(objects, function(x) methods::is(x, "EPhysContinuous"), logical(1)))) {
    stop("All elements of 'objects' must be 'EPhysContinuous'.")
  }

  common <- .merge_EPhysContainer_common(objects,
                                         match_by = match_by,
                                         prefixes = prefixes)

  if (!identical(common$base_class, "EPhysContinuous")) {
    warning("Merged container base_class is ", common$base_class,
            "; proceeding but this was expected to be 'EPhysContinuous'.")
  }

  md_ref   <- common$metadata
  idx_list <- common$idx_list
  new_ch   <- common$channels
  old_ch   <- common$old_channels
  pos_list <- common$positions_list

  arr_ref <- objects[[1L]]@Data
  if (!is.array(arr_ref) || length(dim(arr_ref)) != 3L) {
    stop("EPhysContinuous@Data must be a 3D numeric array [time × run × channel].")
  }

  n_time <- dim(arr_ref)[1L]
  n_runs <- nrow(md_ref)
  n_tot  <- length(new_ch)

  out <- array(NA_real_, dim = c(n_time, n_runs, n_tot))

  dn         <- dimnames(arr_ref)
  time_names <- if (length(dn) >= 1L) dn[[1L]] else NULL
  run_names  <- if (length(dn) >= 2L) dn[[2L]] else NULL

  for (j in seq_along(objects)) {
    obj   <- objects[[j]]
    arr_j <- obj@Data
    if (!is.array(arr_j) || length(dim(arr_j)) != 3L) {
      stop("All objects must have 3D numeric Data arrays.")
    }
    if (!isTRUE(all.equal(dim(arr_j)[1L], n_time))) {
      stop("Time dimension of Data differs for object ", j, ".")
    }

    idx_runs <- idx_list[[j]]
    if (dim(arr_j)[2L] < length(idx_runs)) {
      stop("Object ", j, " has fewer runs than required by 'match_by' alignment.")
    }
    arr_j_aligned <- arr_j[, idx_runs, , drop = FALSE]
    if (!identical(dim(arr_j_aligned)[2L], n_runs)) {
      stop("After alignment, run dimension mismatch for object ", j, ".")
    }

    pos    <- pos_list[[j]]
    n_ch_j <- length(old_ch[[j]])
    if (!identical(length(pos), n_ch_j) ||
        !identical(dim(arr_j_aligned)[3L], n_ch_j)) {
      stop("Channel dimension of object ", j,
           " does not match length(Channels(object)).")
    }

    out[, , pos] <- arr_j_aligned
  }

  dimnames(out) <- list(
    time    = time_names,
    run     = if (!is.null(run_names)) run_names else rownames(md_ref),
    channel = new_ch
  )

  newEPhysContinuous(
    Data             = out,
    TimeTrace        = common$time_trace,
    Metadata         = md_ref,
    Channels         = new_ch,
    Channel_Metadata = common$channel_metadata,
    StimulusTrace    = common$stimulus_trace,
    TimeUnits        = common$time_units,
    StimulusUnits    = common$stimulus_units,
    ExamInfo         = common$exam_info,
    SubjectInfo      = common$subject_info,
    Imported         = common$imported
  )
}

## -------------------------------------------------------------------------
## EPhysEvents: merge along channel dimension
## -------------------------------------------------------------------------

#' Merge multiple EPhysEvents objects along channels
#'
#' @param objects  List of EPhysEvents objects.
#' @param match_by Columns in Metadata used to align runs.
#' @param prefixes Optional character vector of prefixes for channel names
#'   (one per object). If NULL, list names / "set1", "set2", … are used.
#' @export
MergeEvents <- function(objects,
                        match_by = c("Step", "Experiment", "Run"),
                        prefixes = NULL) {
  if (!is.list(objects) || length(objects) == 0L) {
    stop("'objects' must be a non-empty list.")
  }
  if (!all(vapply(objects, function(x) methods::is(x, "EPhysEvents"), logical(1)))) {
    stop("All elements of 'objects' must be 'EPhysEvents'.")
  }

  common <- .merge_EPhysContainer_common(objects,
                                         match_by = match_by,
                                         prefixes = prefixes)

  if (!identical(common$base_class, "EPhysEvents")) {
    warning("Merged container base_class is ", common$base_class,
            "; proceeding but this was expected to be 'EPhysEvents'.")
  }

  md_ref          <- common$metadata
  idx_list        <- common$idx_list
  new_ch          <- common$channels
  old_ch          <- common$old_channels
  new_byobj       <- common$new_channels_by_obj
  rename_map_list <- common$rename_map_list

  n_runs      <- nrow(md_ref)
  merged_data <- vector("list", length = n_runs)

  for (r in seq_len(n_runs)) {
    per_obj_lists <- vector("list", length(objects))

    for (j in seq_along(objects)) {
      obj    <- objects[[j]]
      data_j <- obj@Data

      idx_run <- idx_list[[j]][r]
      if (idx_run < 1L || idx_run > length(data_j)) {
        stop("Run index ", idx_run, " out of range for object ", j, ".")
      }

      run_list <- data_j[[idx_run]]
      if (!is.list(run_list)) {
        stop("Data[[", idx_run, "]] of object ", j,
             " is not a list of per-channel event vectors.")
      }

      nm <- names(run_list)
      if (is.null(nm) || !length(nm)) {
        ## Fall back: assume order equals Channels(object)
        ch_expected <- old_ch[[j]]
        if (length(run_list) != length(ch_expected)) {
          stop("Cannot infer channel names for object ", j,
               " in run ", r, " (unnamed channels with unexpected length).")
        }
        names(run_list) <- new_byobj[[j]]
      } else {
        missing <- setdiff(nm, names(rename_map_list[[j]]))
        if (length(missing)) {
          stop("Some channels in Data[[", idx_run, "]] of object ", j,
               " are not listed in Channels(object): ",
               paste(missing, collapse = ", "))
        }
        names(run_list) <- unname(rename_map_list[[j]][nm])
      }

      per_obj_lists[[j]] <- run_list
    }

    merged_data[[r]] <- do.call(c, per_obj_lists)
  }

  if ("RunUID" %in% names(md_ref)) {
    names(merged_data) <- as.character(md_ref$RunUID)
  } else {
    names(merged_data) <- sprintf("Run%03d", seq_len(n_runs))
  }

  newEPhysEvents(
    Data             = merged_data,
    TimeTrace        = common$time_trace,
    Metadata         = md_ref,
    Channels         = new_ch,
    Channel_Metadata = common$channel_metadata,
    StimulusTrace    = common$stimulus_trace,
    TimeUnits        = common$time_units,
    StimulusUnits    = common$stimulus_units,
    ExamInfo         = common$exam_info,
    SubjectInfo      = common$subject_info,
    Imported         = common$imported
  )
}
