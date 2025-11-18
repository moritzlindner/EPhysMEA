#' Trial aggregation helpers for EPhys* objects
#'
#' @title TrialAggregation: Average or flatten repeated trials by RecordingID
#'
#' @name TrialAggregation
#' @description
#' Two complementary helpers to aggregate repeated trials per recording:
#' \itemize{
#'   \item \code{AverageTrials()} — For \strong{EPhysContinuous}: averages binned/continuous data
#'         across trials within each \code{RecordingID}.
#'   \item \code{FlattenTrials()} — For \strong{EPhysEvents}: concatenates (optionally sorts)
#'         spike/event timestamps across trials within each \code{RecordingID}.
#' }
#'
#' @section Grouping:
#' Both functions group by \code{Metadata$RecordingID}. If this column is absent,
#' all trials are treated as a single group and a warning is issued.
#'
#' @return
#' \describe{
#'   \item{\code{AverageTrials(X)}}{A new \code{EPhysContinuous} with one averaged
#'         trial per \code{RecordingID}.}
#'   \item{\code{FlattenTrials(X)}}{A new \code{EPhysEvents} with one row per
#'         \code{RecordingID}; each channel holds the concatenated (optionally
#'         sorted) spike times from all trials in the group.}
#' }
#' @importClassesFrom EPhysData EPhysContinuous EPhysEvents
#' @importFrom EPhysData Metadata Channels
#' @importFrom methods selectMethod new
#' @seealso \linkS4class{EPhysContinuous}, \linkS4class{EPhysEvents}
#' @examples
#' \dontrun{
#'   x_avg  <- AverageTrials(x_cont)    # continuous data
#'   x_flat <- FlattenTrials(x_events)  # event data
#' }
NULL


#' Average Trials in EPhysContinuous (grouped by RecordingID)
#'
#' Computes the mean spike count across trials *within each RecordingID* for
#' each time bin and channel, and returns a new \code{EPhysContinuous} object
#' containing one averaged trial per unique \code{RecordingID}.
#'
#' @rdname TrialAggregation
#' @inheritParams AddStimulus
#' @inheritParams EPhysData::group_apply
#' @inheritDotParams base::rowMeans
#' @importFrom EPhysData group_apply nested2array
setGeneric("AverageTrials",
           function(X,
                    parallel = FALSE,
                    progress = interactive(),
                    ...)
             standardGeneric("AverageTrials"))
#' @noRd
#' @export
setMethod("AverageTrials",
          signature(X = "EPhysContinuous"),
          function(X,
                   parallel = FALSE,
                   progress = interactive(),
                   ...) {
            nested <-
              group_apply(X,
                          parallel = parallel,
                          progress = progress,
                          function(mat, ...) {
                            rowMeans(mat, ...)
                          })

            # 2) Pack back to array matching EPhysContinuous@Data order: [time × trial × channel]
            out <- nested2array(nested)

            # 3) Build condensed Metadata: one row per RecordingID (first occurrence)
            md         <- Metadata(X)
            rec_levels <- unique(md$RecordingID)
            keep_rows  <- match(rec_levels, md$RecordingID)
            new_md     <- md[keep_rows, , drop = FALSE]

            dimnames(out)   <- NULL

            # 5) Construct the averaged object
            newEPhysContinuous(
              Data             = out,
              TimeTrace        = X@TimeTrace,
              Metadata         = new_md,
              Channels         = Channels(X),
              Channel_Metadata = X@Channel_Metadata,
              StimulusTrace    = X@StimulusTrace,
              TimeUnits        = X@TimeUnits,
              StimulusUnits    = X@StimulusUnits,
              ExamInfo         = X@ExamInfo,
              SubjectInfo      = X@SubjectInfo,
              Imported         = X@Imported
            )
          })

#' Flatten Trials in EPhysEvents (grouped by RecordingID)
#'
#' Concatenates spike/event timestamps per channel across trials that share the
#' same \code{RecordingID}, returning one row per unique \code{RecordingID}.
#'
#' @rdname TrialAggregation
#' @inheritParams Bin
#' @param sort Logical; if \code{TRUE} (default), sort timestamps within each channel after concatenation.
#' @export
setGeneric("FlattenTrials", function(X, ...) standardGeneric("FlattenTrials"))

#' @rdname TrialAggregation
#' @noRd
#' @export
setMethod("FlattenTrials",
          signature(X = "EPhysEvents"),
          function(X, sort = TRUE) {
            md  <- Metadata(X)
            rec <- md$RecordingID
            rec_levels <-
              unique(rec)                       # preserve order of first appearance

            # 1) Concatenate runs within each (RecordingID, Channel)
            #    mat: named list of per-run timestamp vectors -> single numeric vector
            nested <- group_apply(X,
                                  function(mat, sort) {
                                    vals <- unlist(mat, use.names = FALSE)
                                    if (length(vals)) {
                                      if (!is.numeric(vals))
                                        stop("Non-numeric timestamps encountered.")
                                      if (isTRUE(sort))
                                        vals <- sort(vals, na.last = NA)
                                      as.numeric(vals)
                                    } else {
                                      numeric(0L)
                                    }
                                  },
                                  sort = sort)

            # 2) Reorder to match first-occurrence RecordingID order
            nested <- nested[as.character(rec_levels)]

            # 3) One metadata row per RecordingID (first occurrence)
            keep_rows <- match(rec_levels, rec)
            new_md    <- md[keep_rows, , drop = FALSE]

            # 4) Top-level names: use RunUID of the kept rows
            names(nested) <- as.character(new_md$RunUID)

            # 5) Construct new EPhysEvents
            newEPhysEvents(
              Data             = nested,
              # list[[RunUID]][[Channel]] -> numeric vector
              TimeTrace        = X@TimeTrace,
              Metadata         = new_md,
              Channels         = Channels(X),
              Channel_Metadata = X@Channel_Metadata,
              StimulusTrace    = X@StimulusTrace,
              TimeUnits        = X@TimeUnits,
              StimulusUnits    = X@StimulusUnits,
              ExamInfo         = X@ExamInfo,
              SubjectInfo      = X@SubjectInfo,
              Imported         = X@Imported
            )
          })



#' Container-level dispatcher for trial aggregation
#'
#' Calls the appropriate aggregation depending on the subclass:
#' \itemize{
#'   \item \strong{EPhysEvents}: dispatches to \code{FlattenTrials(X, ...)} (concatenates timestamps).
#'   \item \strong{EPhysContinuous}: dispatches to \code{AverageTrials(X, ...)} (averages bins).
#' }
#'
#' @rdname TrialAggregation
#' @param X An \code{EPhysContainer} (or subclass) object.
#' @param ... Passed to the subclass-specific method:
#'   \itemize{
#'     \item For \code{EPhysEvents}: forwarded to \code{FlattenTrials()} (e.g., \code{sort=}).
#'     \item For \code{EPhysContinuous}: forwarded to \code{AverageTrials()} and then \code{mean()}.
#'   }
#' @noRd
#' @export
setMethod("FlattenTrials",
          signature(X = "EPhysContainer"),
          function(X, ...) {

            # EPhysEvents → FlattenTrials(X, ...)
            if (is(X, "EPhysEvents")) {
              # Call the more specific method explicitly to avoid accidental recursion
              ev_mth <- selectMethod("FlattenTrials", signature = c(X = "EPhysEvents"))
              return(ev_mth(X, ...))
            }

            # EPhysContinuous → AverageTrials(X, ...)
            if (is(X, "EPhysContinuous")) {
              dots <- list(...)
              if ("sort" %in% names(dots)) {
                # 'sort' is not meaningful for averaging; drop with a heads-up
                warning("Argument 'sort' is ignored for EPhysContinuous; using AverageTrials() which forwards ... to mean().")
                dots$sort <- NULL
              }
              return(do.call(AverageTrials, c(list(X = X), dots)))
            }

            stop("Unsupported subclass for FlattenTrials(): ",
                 paste(class(X), collapse = "/"),
                 ". Expected an EPhysEvents or EPhysContinuous object.")
          })
