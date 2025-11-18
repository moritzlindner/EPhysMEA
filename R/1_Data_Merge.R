#' Merge two EPhysContainer objects along channels
#'
#' Concatenate channels from two \code{EPhysContainer} objects after aligning
#' **Runs** (rows in \code{Metadata}) by key columns given in \code{match_by}
#' (default \code{c("Step","Experiment")}). Supports two data layouts:
#' (1) list-of-arrays where the 3rd dim is channels, and
#' (2) list-of-lists where the 2nd-level list enumerates channels.
#' Both inputs must use the same layout.
#'
#' @param x,y \code{EPhysContainer} objects.
#' @param match_by Character vector of column names in the Runs (Metadata) used
#'   to align Runs. Defaults to \code{c("Step","Experiment")}.
#' @return A merged \code{EPhysContainer}.
#' @importClassesFrom EPhysData EPhysContainer
setGeneric("Merge", function(x, y, match_by = c("Step","Experiment"))
  standardGeneric("Merge"))

setMethod("Merge", signature(x="EPhysContainer", y="EPhysContainer"),
          function(x, y, match_by = c("Step","Experiment")) {
            .key <- function(md, by) {
              if (!all(by %in% names(md)))
                stop("match_by columns not found in Runs (Metadata): ",
                     paste(setdiff(by, names(md)), collapse=", "))
              do.call(paste, c(lapply(md[by], as.character), sep = "\r"))
            }

            align_by_key <- function(runs_x, runs_y, by) {
              kx <- .key(runs_x, by); ky <- .key(runs_y, by)
              if (anyDuplicated(kx)) stop("Keys in 'x' Runs are not unique.")
              if (anyDuplicated(ky)) stop("Keys in 'y' Runs are not unique.")
              idx <- match(kx, ky)
              if (anyNA(idx)) stop("Some Runs in 'x' were not found in 'y' by match_by.")
              idx
            }

            # Determine data layout
            layout_of <- function(d) {
              if (!is.list(d) || length(d) == 0L) return("list")  # fall back; will catch later
              d1 <- d[[ min(1L, length(d)) ]]
              if (is.array(d1)) return("list_of_arrays")
              if (is.list(d1))  return("list_of_lists")
              "unknown"
            }

            # Align Channel_Metadata rows to channels (by rownames or 'Channel' column)
            align_chmeta_to_channels <- function(chm, channels, for_y = FALSE, rename_map = NULL) {
              if (!is.data.frame(chm) || nrow(chm) == 0L) {
                out <- data.frame(row.names = channels)
                return(out)
              }
              # If y-side and a rename_map is given, rename its rownames/Channel before aligning
              if (for_y && !is.null(rename_map)) {
                if (!is.null(rownames(chm))) {
                  rn <- rownames(chm)
                  rn[match(names(rename_map), rn)] <- unname(rename_map[names(rename_map)])
                  rownames(chm) <- rn
                } else if ("Channel" %in% names(chm)) {
                  chm$Channel <- ifelse(chm$Channel %in% names(rename_map),
                                        rename_map[chm$Channel], chm$Channel)
                }
              }

              if (!is.null(rownames(chm)) && setequal(rownames(chm), channels)) {
                chm <- chm[match(channels, rownames(chm)), , drop = FALSE]
              } else if ("Channel" %in% names(chm) && setequal(chm$Channel, channels)) {
                chm <- chm[match(channels, chm$Channel), , drop = FALSE]
                rownames(chm) <- channels
              } else if (nrow(chm) == length(channels)) {
                rownames(chm) <- channels
              } else {
                stop("Unable to align Channel_Metadata to channels (incomplete or mismatched rows).")
              }
              chm
            }

            # helper: merge two info lists; if values differ, mark as "EPhysContainer"
            merge_info_lists_with_sentinel <- function(ax, ay, sentinel = "EPhysContainer") {
              `%||%` <- function(a,b) if (is.null(a)) b else a
              keys <- union(names(ax %||% list()), names(ay %||% list()))
              out <- vector("list", length(keys)); names(out) <- keys
              for (nm in keys) {
                vx <- ax[[nm]]; vy <- ay[[nm]]
                if (isTRUE(identical(vx, vy))) out[[nm]] <- vx else out[[nm]] <- sentinel
              }
              out
            }

            # Union two data.frames by rows (align columns, fill missing with NA)
            rbind_union <- function(a, b) {
              cols <- union(names(a), names(b))
              add_missing <- function(df) {
                missing <- setdiff(cols, names(df))
                for (m in missing) df[[m]] <- NA
                df[cols]
              }
              a <- add_missing(a)
              b <- add_missing(b)
              out <- rbind(a, b)
              rownames(out) <- c(rownames(a), rownames(b))
              out
            }

            # Add difference columns from info lists; character -> factor
            add_info_diff_cols <- function(chm, chx, chy_renamed, info_x, info_y, prefix) {
              nx <- length(chx); ny <- length(chy_renamed)
              pos_x <- seq_len(nx)
              pos_y <- nx + seq_len(ny)
              keys <- union(names(info_x %||% list()), names(info_y %||% list()))
              for (nm in keys) {
                vx <- info_x[[nm]]
                vy <- info_y[[nm]]
                # treat NULL as NA; skip if identical (including both NULL)
                if (isTRUE(identical(vx, vy))) next
                colname <- make.unique(c(names(chm), paste0(prefix, ".", nm)))[ncol(chm) + 1L]
                vec <- rep(NA, nx + ny)
                vec[pos_x] <- if (length(vx)) vx else NA
                vec[pos_y] <- if (length(vy)) vy else NA
                # coerce character to factor
                if (is.character(vec) || is.character(vx) || is.character(vy)) {
                  chm[[colname]] <- factor(vec)
                } else {
                  chm[[colname]] <- vec
                }
              }
              chm
            }
            `%||%` <- function(a,b) if (is.null(a)) b else a

            ## --- checks & alignment -------------------------------------------------
            # Same layout?
            lay_x <- layout_of(x@Data)
            lay_y <- layout_of(y@Data)
            if (!lay_x %in% c("list_of_arrays", "list_of_lists"))
              stop("Unsupported data layout in 'x'. Must be list-of-arrays or list-of-lists.")
            if (lay_y != lay_x)
              stop("Both inputs must use the same data layout. Got x=", lay_x, ", y=", lay_y, ".")

            # Align Y rows to X by match_by
            idx_y <- align_by_key(x@Metadata, y@Metadata, match_by)

            # Basic invariants that must match
            if (!isTRUE(all.equal(x@TimeTrace, y@TimeTrace)))
              stop("TimeTrace differs between 'x' and 'y'.")
            if (!identical(x@TimeUnits, y@TimeUnits))
              stop("TimeUnits differ between 'x' and 'y'.")
            # Stimulus: allow both empty; if present, must match
            if (length(x@StimulusTrace) != length(y@StimulusTrace) ||
                (length(x@StimulusTrace) > 0 && !isTRUE(all.equal(x@StimulusTrace, y@StimulusTrace))))
              stop("StimulusTrace differs between 'x' and 'y'.")
            if (!identical(x@StimulusUnits, y@StimulusUnits))
              stop("StimulusUnits differ between 'x' and 'y'.")

            # Reorder y to match x
            md_x <- x@Metadata
            md_y <- y@Metadata[idx_y, , drop = FALSE]
            data_x <- x@Data
            data_y <- y@Data[idx_y]
stop("this method was found not to be required currently, maybe helpful in future. needs to make unique channel names w reasonable prefix and rename in cannels column as well cas matrix/list names")
            # Prepare channel names and potential renaming for y
            chx <- x@Channels
            chy <- y@Channels
            new_all <- make.unique(c(chx, chy), sep = ".y")
            chy_renamed <- new_all[(length(chx) + 1L):length(new_all)]
            new_channels <- c(chx, chy_renamed)
            # map original y channel -> renamed
            y_ch_map <- setNames(chy_renamed, chy)

            ## --- merge Data ---------------------------------------------------------
            if (lay_x == "list_of_lists") {
              # Each trial is a list of channels
              merged_data <- Map(function(cx, cy) {
                # rename y channel names, then concatenate and reorder to new_channels
                names(cy) <- y_ch_map[names(cy)]
                out <- c(cx, cy)
                out[new_channels]  # ensure order
              }, data_x, data_y)

            } else if (lay_x == "list_of_arrays") {
              # Each trial is an array [time × something × channel]
              merged_data <- mapply(function(ax, ay) {
                if (!is.array(ax) || !is.array(ay) || length(dim(ax)) < 3L || length(dim(ay)) < 3L)
                  stop("Expected per-trial arrays with at least 3 dimensions.")
                dx <- dim(ax); dy <- dim(ay)
                if (!identical(dx[1:2], dy[1:2]))
                  stop("Non-channel dimensions of per-trial arrays do not match between 'x' and 'y'.")

                kx <- dx[3L]; ky <- dy[3L]
                out <- array(NA_real_, dim = c(dx[1L], dx[2L], kx + ky))
                out[, , seq_len(kx)] <- ax
                out[, , kx + seq_len(ky)] <- ay

                # set dimnames, with channel names = new_channels
                dnx <- dimnames(ax)
                dny <- dimnames(ay)
                # derive time/trial names from x if present
                tnames <- if (!is.null(dnx)) dnx[[1]] else NULL
                rnames <- if (!is.null(dnx) && length(dnx) >= 2L) dnx[[2]] else NULL
                dimnames(out) <- list(tnames, rnames, new_channels)
                out
              }, data_x, data_y, SIMPLIFY = FALSE)

            } else {
              stop("Unexpected layout guard.")
            }

            ## --- merge Channel_Metadata ---------------------------------------------
            chm_x <- align_chmeta_to_channels(x@Channel_Metadata, chx)
            chm_y <- align_chmeta_to_channels(y@Channel_Metadata, chy, for_y = TRUE, rename_map = y_ch_map)

            # tag provenance
            chm_x$Source <- factor(rep("x", length(chx)), levels = c("x","y"))
            chm_y$Source <- factor(rep("y", length(chy)), levels = c("x","y"))

            # union & bind
            chm <- rbind_union(chm_x, chm_y)
            chm <- chm[ new_channels, , drop = FALSE ]  # ensure row order

            # Differences from SubjectInfo and ExamInfo -> columns
            s_x <- x@SubjectInfo; s_y <- y@SubjectInfo
            r_x <- x@ExamInfo
            r_y <- y@ExamInfo

            chm <- add_info_diff_cols(chm, chx, chy_renamed, s_x, s_y, prefix = "SubjectInfo")
            chm <- add_info_diff_cols(chm, chx, chy_renamed, r_x, r_y, prefix = "ExamInfo")

            # Ensure rownames are channels
            rownames(chm) <- new_channels

            # Build merged SubjectInfo / ExamInfo using sentinel for differing fields
            subj_merged <- merge_info_lists_with_sentinel(x@SubjectInfo, y@SubjectInfo)
            exam_merged <- merge_info_lists_with_sentinel(x@ExamInfo, y@ExamInfo)

            ## --- construct result ----------------------------------------------------
            out <- new("EPhysContainer",
                       Metadata        = md_x,
                       Data            = merged_data,
                       ExamInfo        = exam_merged,
                       SubjectInfo     = subj_merged,
                       Imported        = x@Imported,
                       TimeTrace       = x@TimeTrace,
                       Channels        = new_channels,
                       Channel_Metadata= chm,
                       StimulusTrace   = x@StimulusTrace,
                       TimeUnits       = x@TimeUnits,
                       StimulusUnits   = x@StimulusUnits
            )
            validObject(out)
            out
          }
)
