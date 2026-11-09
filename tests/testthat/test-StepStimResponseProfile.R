test_that("Liang preference preserves the sign of both response transitions", {
  score <- function(on, off, polarity) {
    step_stim_response_liang_index(
      x = c(seq_len(on) / (on + 1), 1 + seq_len(off) / (off + 1)),
      step_onset = 0,
      step_offset = 1,
      baseline_rate = 2,
      profile_polarity = polarity,
      binwidth = 1
    )$Liang_OnOff_PreferenceIndex
  }

  expect_equal(score(4, 2, "positive"), 1)
  expect_equal(score(0, 2, "negative"), -1)
  expect_equal(score(2, 4, "positive"), -1)
  expect_equal(score(2, 0, "ambiguous"), 1)
  expect_equal(score(1, 0, "negative"), 1 / 3)
  expect_identical(score(2, 2, "ambiguous"), NA_real_)
})

test_that("Liang preference omits positive rebounds in suppressive responses", {
  score <- function(polarity) {
    step_stim_response_liang_index(
      x = c(0.25, 0.6, 0.75, 0.9, 1.125, 1.375, 1.625, 1.875),
      step_onset = 0,
      step_offset = 1,
      baseline_rate = 4,
      profile_polarity = polarity,
      binwidth = 0.5
    )$Liang_OnOff_PreferenceIndex
  }

  expect_equal(score("negative"), -1)
  expect_equal(score("ambiguous"), 0)
})

test_that("the recording-level Liang preference uses the average firing traces", {
  spike_train <- function(on, off) {
    c(0.25, 0.75, 1 + seq_len(on) / (on + 1), 2 + seq_len(off) / (off + 1))
  }
  md <- data.frame(
    Experiment = factor(c("test", "test")),
    Step = ordered(c("1", "1"), levels = "1"),
    RecordingID = c(1L, 1L),
    RunUID = c(1L, 2L)
  )
  x <- EPhysData::newEPhysEvents(
    Data = list(
      list(A = spike_train(6, 2), B = spike_train(2, 0)),
      list(B = spike_train(2, 0), A = spike_train(2, 4))
    ),
    TimeTrace = numeric(0),
    Metadata = md,
    Channels = c("A", "B"),
    Channel_Metadata = data.frame(row.names = c("A", "B")),
    TimeUnits = "s"
  )
  profile <- function(average) {
    StepStimResponseProfile(
      x,
      step_range = c(1, 2),
      average_by_recordingID = average,
      baseline_win = 1,
      on_win = 1,
      off_win = 1,
      binwidth = 1,
      parallel = FALSE
    )
  }

  per_run <- profile(FALSE)
  expect_equal(per_run$Liang_OnOff_PreferenceIndex[per_run$Channel == "A"], c(1, -1))
  expect_equal(per_run$Liang_OnOff_PreferenceIndex[per_run$Channel == "B"], c(1, 1))

  averaged <- profile(TRUE)
  expect_equal(averaged$Liang_OnOff_PreferenceIndex[averaged$Channel == "A"], 1 / 3)
  expect_equal(averaged$Liang_OnOff_PreferenceIndex[averaged$Channel == "B"], 1)
})
