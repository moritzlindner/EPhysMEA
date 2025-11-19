#' Import Spike Channel Variables from a Sorted Spikes MATLAB File
#'
#' Loads all variables from a `_sortedspikes.mat` file whose names start with `"sig."`,
#' and returns them as a named list. This is commonly used to import spike waveforms
#' or timestamps saved under such variable names.
#'
#' @param file_root A character string indicating the root of the file path (excluding extension).
#'                  The function will look for `[file_root]_sortedspikes.mat`.
#'
#' @return A named list containing all variables whose names start with `"sig."`.
#'         If no matching variables are found, a warning is issued and an empty list is returned.
#'
#' @details
#' The function uses `R.matlab::readMat()` to load the entire MATLAB file into memory
#' and filters for variable names that start with `"sig."`. Note that `readMat()` does not support
#' header-only loading, so the full file is read even if only a subset of variables is returned.
#'
#' @examples
#' \dontrun{
#'   spike_data <- Import_ChannelData("data/session1")
#'   names(spike_data)  # Shows all channels like "sig.1", "sig.2", etc.
#' }
#'
#' @importFrom R.matlab readMat
#' @export
Import_ChannelData <- function(file_root) {
  spike_file <- paste0(file_root, "_sortedspikes.mat")

  if (!file.exists(spike_file)) {
    warning("Spike data file not found: ", spike_file)
    return(data.frame())
  }

  # Use R.matlab to read header only (to inspect variable names)
  # Note: R.matlab does not support `whos` equivalent, so we use readMat but stop loading large data
  # This implementation loads full data, then drops it; you can optimize with .mat version >= v7.3 and rhdf5
  mat_data <- readMat(spike_file)
  all_vars <- names(mat_data)

  # Identify variables that start with "sig_"
  sig_vars <- all_vars[startsWith(all_vars, "sig.")]

  if (length(sig_vars) == 0) {
    warning('No spike channel variables (starting with "sig_") found in file: ', spike_file)
    return(data.frame())
  }

  mat_data <- lapply(mat_data[sig_vars], function(x){as.vector(x)})

  return(mat_data)
}
