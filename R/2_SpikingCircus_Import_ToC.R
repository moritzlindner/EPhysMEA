#' Import Event Table with Optional Single-Trial Splitting
#'
#' Loads event data from a MATLAB `_events.mat` file, processes it via a type map,
#' and—if requested—splits each event into sub‐events based on `evt01` timestamps.
#'
#' @param file_root A character string specifying the file root
#'                  (path + base name without extension).
#' @param type_assignment_file Path to a tab-delimited text file with two columns:
#'        Experiment names (col 1) and corresponding numeric IDs (col 2).
#' @param SingleTrialTimings Logical; if `TRUE`, each event will be subdivided
#'        by the `evt01` timestamps into smaller “trials.” Default `FALSE`.
#'
#' @return A single data frame with columns:
#'   \describe{
#'     \item{Intensity}{Event intensity.}
#'     \item{Trials}{Original trial index from the file.}
#'     \item{Experiment}{Factor representing individual experiments (labels from the map).}
#'     \item{Repeat}{Repeat index inferred from fractional part of raw type.}
#'     \item{Start}{Start time of the interval.}
#'     \item{Stop}{Stop time of the interval.}
#'     \item{Diff}{Duration (`Stop - Start`).}
#'     \item{Trial}{If `SingleTrialTimings=TRUE`, the sub‐trial index (1,2,…) within each event;
#'                   otherwise absent.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # only main events:
#'   df_main <- Import_ToC("data/session1", "types.tsv", FALSE)
#'
#'   # only sub‐trials:
#'   df_sub  <- Import_ToC("data/session1", "types.tsv", TRUE)
#' }
#'
#' @importFrom R.matlab readMat
#' @export
Import_ToC <- function(file_root,
                       type_assignment_file,
                       SingleTrialTimings = FALSE) {
  if(!is.character(file_root) || length(file_root)!=1){
    stop("'file_root' must be a character string scalar.")
  }
  events_file <- paste0(file_root, "_events.mat")
  if (!file.exists(events_file)) {
    warning("Events file not found: ", events_file)
    return(data.frame())
  }

  mat_data <- readMat(events_file)
  if (!"Events" %in% names(mat_data)) {
    warning("No 'Events' object found in file: ", events_file)
    return(data.frame())
  }

  events <- mat_data$Events[, 1, ]
  req <- c("intensity", "trials", "type", "start", "stop")
  if (!all(req %in% names(events))) {
    missing <- req[!req %in% names(events)]
    warning("Missing fields in Events: ", paste(missing, collapse = ", "))
    return(data.frame())
  }

  # pull vectors
  intensity <- as.vector(events$intensity)
  trials    <- as.vector(events$trials)
  type_raw  <- as.vector(events$type)
  start     <- as.vector(events$start)
  stop      <- as.vector(events$stop)

  step      <- as.ordered(intensity)

  # compute repeat & clean type
  intp  <- floor(type_raw)
  rep_i <- round((type_raw - intp) * 10)
  rep_i[rep_i == 0] <- 1
  type_clean <- intp

  # build the main events DF
  df_ev <- data.frame(
    Intensity = intensity,
    Step = step,
    Runs    = trials,
    Experiment      = type_clean,
    Repeat    = rep_i,
    Start     = start,
    Stop      = stop,
    Diff      = stop - start,
    stringsAsFactors = FALSE
  )
  df_ev <- subset(df_ev, Stop != 0)
  df_ev$RecordingID<-as.integer(1:nrow(df_ev))

  # apply type mapping
  map <- read.table(type_assignment_file,
                    sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE)
  names(map) <- c("TypeName", "TypeID")

  notfound<-(!(df_ev$Experiment %in% unique(map$TypeID)))
  if (sum(notfound)>0){
    warnstr<-unique(df_ev$Experiment[notfound])
    warning("The following Experiment identifieres could not be matched: " ,paste(warnstr,","))
  }

  df_ev$Experiment <- factor(df_ev$Experiment,
                       levels = map$TypeID,
                       labels = map$TypeName)

  df_ev <- df_ev[!notfound,]

  if (!SingleTrialTimings) {
    # no Trial column at all
    return(df_ev)
  }

  # extract evt01 timestamps
  evt_all <- numeric(0)
  if ("evt01" %in% names(events)) {
    tmp <- events$evt01
    evt_all <- if (is.list(tmp) && length(tmp) == 1) as.vector(tmp[[1]]) else as.vector(tmp)
  }

  # build sub-trial rows only
  df_sub <- do.call(rbind, lapply(seq_len(nrow(df_ev)), function(i) {
    er <- df_ev[i, ]
    stamps <- evt_all[evt_all >= er$Start & evt_all <= er$Stop]
    if (length(stamps) < 2) {
      er$Run <- 1
      er
    } else {
      ends <- c(stamps[-1], er$Stop)
      do.call(rbind, lapply(seq_along(stamps), function(k) {
        row <- er
        row$Start <- stamps[k]
        row$Stop  <- ends[k]
        row$Diff  <- ends[k] - stamps[k]
        if(all(k%%1==0)){
          row$Run <- as.integer(k)
        } else {
          row$Run <- k
          warning("Run IDs are not convertible to integer!")
        }
        row
      }))
    }
  }))

  df_sub$RunUID<-as.integer(1:nrow(df_sub))

  rownames(df_sub) <- NULL
  df_sub
}
