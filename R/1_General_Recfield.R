#' Compute receptive fields and r50 metrics for EPhysEvents
#'
#' @description
#' For each run (row of \code{Metadata}) and each channel, this method:
#' \enumerate{
#'   \item builds a receptive-field map with \code{\link{recfield_response_map}}, and
#'   \item extracts peak/radius anisotropy metrics with \code{\link{recfield_r50_radial}}.
#' }
#' It returns either a per-(Run,Channel) table of r50 metrics or, if
#' \code{ReturnMap = TRUE}, additionally a long-form stack of RF maps.
#'
#' A fixed response profile can be supplied in one of two ways. This profile is
#' used to determine which response component contributes to the frame-wise
#' \code{Response_Amplitude} values from which the receptive-field map is built:
#' positive onset responses, positive offset responses, or suppressive onset
#' responses.
#' \enumerate{
#'   \item by \code{profile_win}, a single reference step interval used to infer
#'   \code{profile_polarity} and \code{profile_bias} separately for each channel;
#'   \item by explicit vectors \code{profile_polarity} and \code{profile_bias},
#'   each of length \code{length(Channels(X))}, giving one fixed profile
#'   decision per channel.
#' }
#' The resolved per-channel values are passed to
#' \code{\link{recfield_response_map}}.
#'
#' @param X An \code{EPhysEvents} instance.
#' @inheritParams recfield_response_map
#' @param pixel_width Numeric scalar (e.g., µm per pixel) for
#'   \code{recfield_r50_radial()}.
#' @param profile_win Optional numeric or \code{units} length-2 vector defining a
#'   reference stimulus step from which \code{profile_polarity} and
#'   \code{profile_bias} are inferred for each channel. It is converted with
#'   \code{to_num()} and must lie fully outside the interval covered by
#'   \code{change_times}.
#' @param profile_polarity Optional character vector of length
#'   \code{length(Channels(X))}, giving one fixed polarity decision per channel.
#'   Entries must be one of \code{"positive"}, \code{"negative"}, or
#'   \code{"ambiguous"}. Used only if \code{profile_win} is \code{NULL}.
#' @param profile_bias Optional character vector of length
#'   \code{length(Channels(X))}, giving one fixed bias decision per channel.
#'   Entries must be one of \code{"onset"}, \code{"offset"}, or
#'   \code{"ambiguous"}. Used only if \code{profile_win} is \code{NULL}.
#' @param ReturnMap Logical; if \code{TRUE}, also returns the individual heatmaps.
#' @inheritParams EPhysData::`lapply-EPhys`
#' @param ... Additional arguments forwarded to \code{\link{recfield_response_map}}
#'   and hence to the relevant StepStimResponseProfile low-level helpers.
#'
#' @return
#' If \code{ReturnMap = FALSE}, a wide \code{data.frame} with one row per
#' run × channel containing the receptive-field metrics.
#'
#' If \code{ReturnMap = TRUE}, a named list with:
#' \describe{
#'   \item{\code{maps}}{data.frame with columns \code{RunUID},
#'   \code{RecordingID}, \code{Channel}, \code{x}, \code{y}, \code{value}.}
#'   \item{\code{metrics}}{data.frame with one row per run × channel and columns
#'   including \code{peak_amp}, \code{r50}, \code{symmetry_ratio}, and
#'   \code{peak_at_margin}.}
#' }
#'
#' @importFrom cli cli_progress_done cli_progress_update cli_progress_bar
#' @examples
#' \donttest{
#' # Assuming you already have: E (EPhysEvents), StimSeq, change_times
#' res <- Recfield(
#'   E,
#'   StimSeq = StimSeq,
#'   change_times = change_times,
#'   pixel_width = 10,
#'   profile_polarity = rep("positive", length(Channels(E))),
#'   profile_bias = rep("onset", length(Channels(E)))
#' )
#' head(res)
#' }
#'
#' @export
setGeneric("Recfield", function(X, StimSeq, change_times, pixel_width,
                                baserate_win = NULL,
                                profile_win = NULL,
                                profile_polarity = NULL,
                                profile_bias = NULL,
                                ReturnMap = FALSE,
                                parallel = parallel,
                                error = c("stop", "warn")[1],
                                progress = interactive(), ...) {
  standardGeneric("Recfield")
})

#' @rdname Recfield
#' @export
setMethod("Recfield", signature(X = "EPhysEvents"),
          function(X, StimSeq, change_times, pixel_width,
                   baserate_win = NULL,
                   profile_win = NULL,
                   profile_polarity = NULL,
                   profile_bias = NULL,
                   ReturnMap = FALSE,
                   parallel = parallel,
                   error = c("stop", "warn")[1],
                   progress = interactive(), ...) {

            # --- basic checks ---
            if (length(dim(StimSeq)) != 3L) {
              stop("StimSeq must be a 3D array [rows × cols × frames].")
            }
            if (!is.numeric(change_times) || length(change_times) == 0L) {
              stop("'change_times' must be a numeric vector of frame times (seconds).")
            }
            if (!is.numeric(pixel_width) || length(pixel_width) != 1L ||
                !is.finite(pixel_width) || pixel_width <= 0) {
              stop("'pixel_width' must be a single positive number.")
            }

            change_times <- to_num(change_times)
            profile_win <- to_num(profile_win)

            channels <- Channels(X)

            if (!is.null(profile_win) && (!is.null(profile_polarity) || !is.null(profile_bias))) {
              stop("Provide either 'profile_win' or explicit 'profile_polarity'/'profile_bias', not both.")
            }

            if (!is.null(profile_win)) {
              if (!is.numeric(profile_win) || length(profile_win) != 2L || any(!is.finite(profile_win))) {
                stop("'profile_win' must be a finite numeric (or units) length-2 vector.")
              }
              profile_win <- sort(profile_win)
              if (!(profile_win[2] > profile_win[1])) {
                stop("'profile_win' must satisfy start < end.")
              }
              if (!(profile_win[2] <= min(change_times) || profile_win[1] >= max(change_times))) {
                stop("'profile_win' must lie fully outside the time interval covered by 'change_times'.")
              }
            } else {
              if (!is.null(profile_polarity)) {
                if (!is.character(profile_polarity) || length(profile_polarity) != length(channels)) {
                  stop("'profile_polarity' must be a character vector of length length(Channels(X)).")
                }
                if (!all(profile_polarity %in% c("positive", "negative", "ambiguous"))) {
                  stop("Entries of 'profile_polarity' must be one of: 'positive', 'negative', 'ambiguous'.")
                }
                names(profile_polarity) <- channels
              }

              if (!is.null(profile_bias)) {
                if (!is.character(profile_bias) || length(profile_bias) != length(channels)) {
                  stop("'profile_bias' must be a character vector of length length(Channels(X)).")
                }
                if (!all(profile_bias %in% c("onset", "offset", "ambiguous"))) {
                  stop("Entries of 'profile_bias' must be one of: 'onset', 'offset', 'ambiguous'.")
                }
                names(profile_bias) <- channels
              }

              if (xor(is.null(profile_polarity), is.null(profile_bias))) {
                stop("Provide both 'profile_polarity' and 'profile_bias', or neither.")
              }
            }

            single_worker <- function(spikes,
                                      StimSeq,
                                      change_times,
                                      baserate_win,
                                      profile_win,
                                      profile_polarity,
                                      profile_bias,
                                      ReturnMap,
                                      ...) {
              rf_map <- recfield_response_map(
                spike_times = spikes,
                StimSeq = StimSeq,
                change_times = change_times,
                baserate_win = baserate_win,
                profile_win = profile_win,
                profile_polarity = profile_polarity,
                profile_bias = profile_bias,
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
                return(list(r50 = m, map = map_df))
              } else {
                return(m)
              }
            }

            if (!is.null(profile_win)) {
              res_list <- EPhysData::lapply(
                X,
                single_worker,
                StimSeq,
                change_times,
                baserate_win,
                profile_win,
                NULL,
                NULL,
                ReturnMap,
                parallel = parallel,
                error = error,
                progress = progress
              )
            } else if (!is.null(profile_polarity) && !is.null(profile_bias)) {

              chan_pol <- profile_polarity[channels]
              chan_bias <- profile_bias[channels]

              res_list <- vector("list", length(X@Data))
              names(res_list) <- names(X@Data)

              if (progress) {
                pb_id <- cli::cli_progress_bar(
                  name = "Recfield",
                  total = sum(vapply(X@Data, length, integer(1L)))
                )
                on.exit(try(cli::cli_progress_done(id = pb_id), silent = TRUE), add = TRUE)
              }

              for (i in seq_along(X@Data)) {
                run_obj <- X@Data[[i]]

                run_res <- vector("list", length(run_obj))
                names(run_res) <- names(run_obj)

                for (j in seq_along(run_obj)) {
                  nm <- names(run_obj)[j]

                  run_res[[j]] <- tryCatch(
                    single_worker(
                      spikes = run_obj[[j]],
                      StimSeq = StimSeq,
                      change_times = change_times,
                      baserate_win = baserate_win,
                      profile_win = NULL,
                      profile_polarity = chan_pol[[nm]],
                      profile_bias = chan_bias[[nm]],
                      ReturnMap = ReturnMap,
                      ...
                    ),
                    error = function(e) {
                      msg <- sprintf(
                        "Recfield() error at run %s (index %d), channel '%s': %s",
                        if (is.null(names(X@Data))) as.character(i) else names(X@Data)[i],
                        i,
                        nm,
                        conditionMessage(e)
                      )

                      if (identical(error, "stop")) {
                        stop(msg, call. = FALSE)
                      } else {
                        warning(msg, call. = FALSE)
                        NULL
                      }
                    }
                  )

                  if (progress) {
                    cli::cli_progress_update(id = pb_id)
                  }
                }

                res_list[[i]] <- run_res
              }
            } else {
              res_list <- EPhysData::lapply(
                X,
                single_worker,
                StimSeq,
                change_times,
                baserate_win,
                NULL,
                NULL,
                NULL,
                ReturnMap,
                parallel = parallel,
                error = error,
                progress = progress
              )
            }

            if (ReturnMap) {
              r50 <- lapply(res_list, function(x) lapply(x, function(y) as.data.frame(y$r50)))
              r50 <- nested2df(X, r50)

              rf_map <- lapply(res_list, function(x) lapply(x, function(y) y$map))
              rf_map <- nested_maps_to_df(X, rf_map)

              return(list(
                maps = rf_map,
                metrics = r50
              ))
            } else {
              nested2df(X, res_list)
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
