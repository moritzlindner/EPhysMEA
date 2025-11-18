#' Build a stimulus plot for an EPhysContainer
#'
#' Build a 1D stimulus trace panel (a small line plot) from an
#' \code{EPhysContainer}. If no valid stimulus is present, returns \code{NULL}.
#'
#' The function first tries an internal \code{HasStimulus(X)} helper (if defined).
#' If that is not available, it falls back to checking that \code{StimulusTrace(X)}
#' is non-\code{NULL}, has the same length as \code{TimeTrace(X)}, and contains at
#' least one non-\code{NA} value. The x-axis label uses \code{TimeUnits(X)}.
#' When \code{NoRowLabels = TRUE}, y-axis text and ticks are suppressed so the
#' panel stacks cleanly beneath other plots.
#'
#' @param X An \code{EPhysContainer} object.
#' @param NoRowLabels Logical; if \code{TRUE}, hide y-axis labels/ticks in the stimulus panel.
#'
#' @return A \code{ggplot} object with the stimulus trace, or \code{NULL} if no valid stimulus exists.
#'
#' @examples
#' \dontrun{
#' sp <- ggStimulusPlot(X)         # returns ggplot or NULL
#' if (!is.null(sp)) print(sp)
#'
#' # Example stacking with patchwork (if you have another plot 'p'):
#' # library(patchwork)
#' # if (!is.null(sp)) p / sp + plot_layout(heights = c(3, 1))
#' }
#'
#' @seealso \code{\link{TimeTrace}}, \code{\link{StimulusTrace}}, \code{\link{TimeUnits}}
#'
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs theme element_blank margin scale_x_continuous
#' @importClassesFrom EPhysData EPhysContainer
#' @importFrom EPhysData Metadata Channels ChannelMetadata TimeTrace StimulusTrace HasStimulus
#' @name ggStimulusPlot
NULL
setGeneric("ggStimulusPlot", function(X, NoRowLabels = FALSE)
  standardGeneric("ggStimulusPlot"))

setMethod(
  "ggStimulusPlot", signature(X = "EPhysContainer"),
  function(X, NoRowLabels = FALSE) {
    # Prefer internal helper if present; otherwise, derive from accessors
    has_stim <- HasStimulus(X)

    if (isTRUE(has_stim)) {
      stim_df <- data.frame(
        Time     = to_num(TimeTrace(X)),
        Stimulus = drop_units(StimulusTrace(X))
      )

      stim_plot <-
        ggplot(stim_df, aes(x = Time, y = Stimulus)) +
        geom_line(linewidth = 0.4, linetype = 1, na.rm = TRUE) +
        theme_minimal() +
        labs(x = paste0("Time (", TimeUnits(X), ")"), y = "Stim") +
        theme(
          panel.grid  = element_blank(),
          plot.margin = margin(0, 0, 0, 0)
        ) +
        scale_x_continuous(expand = c(0, 0))

      if (isTRUE(NoRowLabels)) {
        stim_plot <- stim_plot + theme(
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
    }

    stim_plot
  }
)
