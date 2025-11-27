#' Compute receptive fields and r50 metrics for EPhysEvents
#'
#' @description
#' For each run (row of \code{Metadata}) and each channel, this method:
#' \enumerate{
#'   \item builds a receptive-field map with \code{\link{recfield_response_map}}
#'         (signed \eqn{\Delta}rate per pixel, averaged over exposures), and
#'   \item extracts peak/radius anisotropy metrics with \code{\link{recfield_r50_radial}}.
#' }
#' It returns two tidy data.frames: a long-form stack of RF maps and a per-(Run,Channel)
#' table of r50 metrics.
#'
#' @param X An \code{EPhysEvents} X.
#' @inheritParams recfield_response_map
#' @param pixel_width Numeric scalar (e.g., µm per pixel) for \code{recfield_r50_radial()}.
#' @param ReturnMap Logical; if \code{TRUE}, also returns the individual heatmaps.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @param ... Additional arguments forwarded to \code{\link{recfield_response_map}}
#'   (and ultimately to \code{\link{step_stim_resp_metrics}}), e.g. \code{trans_win},
#'   \code{binwidth}, \code{smooth}, etc.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{maps}}{data.frame with columns \code{Run}, \code{Channel}, \code{x}, \code{y}, \code{value}.}
#'   \item{\code{metrics}}{data.frame with columns \code{Run}, \code{Channel}, \code{peak_amp}, \code{r50},
#'         \code{symmetry_ratio}, \code{peak_at_margin}.}
#' }
#' @importFrom cli cli_progress_done cli_progress_update cli_progress_bar
#' @examples
#' \donttest{
#' # Assuming you already have: E (EPhysEvents), StimSeq, change_times
#' res <- Recfield(E, StimSeq = StimSeq, change_times = change_times,
#'                 pixel_width = 10, baserate_win = 1.0, trans_win = 0.4)
#' head(res$maps)
#' head(res$metrics)
#' }
#'
#' @export
setGeneric("Recfield", function(X, StimSeq, change_times, pixel_width,
                                baserate_win = NULL, ReturnMap = FALSE,
                                parallel = parallel,
                                error = c("stop", "warn")[1],
                                progress = interactive(), ...) {
  standardGeneric("Recfield")
})

#' @rdname Recfield
#' @export
setMethod("Recfield", signature(X = "EPhysEvents"),
          function(X, StimSeq, change_times, pixel_width,
                   baserate_win = NULL, ReturnMap = FALSE,
                   parallel = parallel,
                   error = c("stop", "warn")[1],
                   progress = interactive(), ...) {

            # --- basic checks ---
            if (length(dim(StimSeq)) != 3L) stop("StimSeq must be a 3D array [rows × cols × frames].")
            if (!is.numeric(change_times) || length(change_times) == 0L) {
              stop("'change_times' must be a numeric vector of onsets (seconds).")
            }
            if (!is.numeric(pixel_width) || length(pixel_width) != 1L || !is.finite(pixel_width) || pixel_width <= 0) {
              stop("'pixel_width' must be a single positive number.")
            }

            n_runs <- nrow(Metadata(X))
            channels <- Channels(X)

            single_worker <- function(spikes,
                                      StimSeq,
                                      change_times,
                                      baserate_win,
                                      ReturnMap,
                                      ...) {
              rf_map <- recfield_response_map(
                spike_times  = spikes,
                StimSeq      = StimSeq,
                change_times = change_times,
                baserate_win = baserate_win,
                ...
              )

              m <- recfield_r50_radial(rf_map, pixel_width = pixel_width)

              if (ReturnMap) {
                nr <- nrow(rf_map)
                nc <- ncol(rf_map)
                map_df <- data.frame(
                  x = rep(seq_len(nc), each = nr),
                  y = rep(seq_len(nr), times = nc),
                  value = as.vector(rf_map),
                  check.names = FALSE
                )
                return(list(r50=m,map=map_df))
              } else {
                return(m)
              }
            }

            res_list<-lapply(
              X,
              single_worker,
              StimSeq,
              change_times,
              baserate_win,
              ReturnMap,
              parallel = parallel,
              error = error,
              progress = progress
            )
            r50<-lapply(res_list, function(x){lapply(x, function(y){as.data.frame(y$r50)})})
            r50<-nested2df(X,r50)
            if(ReturnMap){
              rf_map<-lapply(res_list, function(x){lapply(x, function(y){y$map})}) #as.data.frame(y$map)
              rf_map<-nested_maps_to_df(X,rf_map)
              return(list(maps    = rf_map,
                          metrics = r50))
            } else {
              r50
            }
          })

#' @keywords internal
nested_maps_to_df <- function(X, nestedlist) {
  md <- Metadata(X)
  ch <- Channels(X)

  run_dfs <- lapply(seq_along(nestedlist), function(i) {
    run <- nestedlist[[i]]
    if (!is.null(names(run))) run <- run[ch]

    rows <- lapply(names(run), function(chan_name) {
      v <- run[[chan_name]]  # long df x,y,value
      v$Channel     <- chan_name
      v$RunUID      <- md$RunUID[i]
      v$RecordingID <- md$RecordingID[i]
      v
    })
    do.call(rbind, rows)
  })

  do.call(rbind, run_dfs)
}
