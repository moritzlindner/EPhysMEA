test_that("MC_Rack raster ASCII parsing preserves sweeps and converts time", {
  td <- tempfile("mcs_raster_")
  dir.create(td)

  writeLines(
    c(
      "MC_Rack Raster Plot",
      "Time [ms]\tSweep\tUnit",
      "100\t1\t0",
      "250\t1\t0",
      "50\t2\t0"
    ),
    file.path(td, "electrode_11.txt")
  )

  writeLines(
    c(
      "MC_Rack Raster Plot",
      "Time [ms]\tSweep\tUnit",
      "125\t1\t0",
      "75\t2\t0"
    ),
    file.path(td, "electrode_12.txt")
  )

  x <- ImportMCS(td, experiment = "retinal_step", time_range = c(0, 1))

  expect_s4_class(x, "EPhysEvents")
  expect_equal(EPhysData::Channels(x), c("electrode_11", "electrode_12"))

  md <- EPhysData::Metadata(x)
  expect_equal(md$Sweep, c(1L, 2L))
  expect_true(is.ordered(md$Step))
  expect_equal(as.character(md$Experiment), rep("retinal_step", 2))

  expect_equal(x@Data[[1]]$electrode_11, c(0.100, 0.250))
  expect_equal(x@Data[[2]]$electrode_11, 0.050)
  expect_equal(x@Data[[1]]$electrode_12, 0.125)
  expect_equal(x@Data[[2]]$electrode_12, 0.075)
})

test_that("named input groups become ordered experimental steps", {
  root <- tempfile("mcs_steps_")
  dir.create(root)
  low <- file.path(root, "low")
  high <- file.path(root, "high")
  dir.create(low)
  dir.create(high)

  writeLines(
    c("Time [s]\tSweep\tUnit", "0.10\t1\t0", "0.20\t2\t0"),
    file.path(low, "electrode_11.txt")
  )
  writeLines(
    c("Time [s]\tSweep\tUnit", "0.15\t1\t0", "0.25\t2\t0"),
    file.path(high, "electrode_11.txt")
  )

  x <- ImportMCS(
    list("0.1" = low, "1" = high),
    experiment = "retinal_step",
    recording_id = 7L,
    time_range = c(0, 1)
  )

  md <- EPhysData::Metadata(x)
  expect_equal(levels(md$Step), c("0.1", "1"))
  expect_equal(as.character(md$Step), c("0.1", "0.1", "1", "1"))
  expect_equal(md$RecordingID, rep(7L, 4))
  expect_equal(md$RunUID, 1:4)
})

test_that("Spike Sorter unit tags can be retained or filtered", {
  td <- tempfile("mcs_units_")
  dir.create(td)

  f <- file.path(td, "electrode_21.txt")
  writeLines(
    c(
      "Time [s]\tSweep\tUnit",
      "0.10\t1\t0",
      "0.20\t1\t1",
      "0.30\t2\t0",
      "0.40\t2\t1"
    ),
    f
  )

  all_units <- ImportMCS(td, time_range = c(0, 1))
  expect_equal(EPhysData::Channels(all_units), c("electrode_21_unit0", "electrode_21_unit1"))

  unit1 <- ImportMCS(td, unit_tags = 1L, time_range = c(0, 1))
  expect_equal(EPhysData::Channels(unit1), "electrode_21")
  expect_equal(unit1@Data[[1]]$electrode_21, 0.20)
  expect_equal(unit1@Data[[2]]$electrode_21, 0.40)
})

test_that("semicolon-delimited decimal-comma raster export is accepted", {
  td <- tempfile("mcs_decimal_comma_")
  dir.create(td)

  writeLines(
    c(
      "Time [ms];Sweep;Unit",
      "100,5;1;0",
      "250,5;1;0"
    ),
    file.path(td, "electrode_31.txt")
  )

  x <- ImportMCS(td, time_range = c(0, 1))
  expect_equal(x@Data[[1]]$electrode_31, c(0.1005, 0.2505))
})
