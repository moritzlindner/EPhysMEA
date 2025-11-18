<!-- README.md is generated from README.Rmd. Please edit that file -->

  # EPhysMEA R Package

  ## Introduction

  The `EPhysMEA` package provides analysis tools for multichannel
electrophysiology recordings stored in the `EPhysData` classes
`EPhysEvents` and `EPhysContinuous`. It adds methods for binning and
trial aggregation, tuning and receptive-field analysis, Naka–Rushton
fitting, spike-train irregularity metrics, and ggplot2-based
visualisation.

`EPhysMEA` is intended as a companion to the `EPhysData` package: raw
data import and basic container handling are done in `EPhysData`, while
`EPhysMEA` focuses on downstream analysis.

## Installation

You can install the development version of `EPhysMEA` from GitHub:

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("EPhysData", quietly = TRUE)) {
  remotes::install_github("moritzlindner/EPhysData")
}
remotes::install_github("moritzlindner/EPhysMEA")

library(EPhysData)
library(EPhysMEA)

## Example: basic rate and variability metrics ---------------------------

# events: EPhysEvents object (runs × channels of spike/event times)
# cont:   EPhysContinuous object (time × runs × channels of traces)

# Mean firing rate per channel

rate_df <- ChannelMeanRate(events)
head(rate_df)
