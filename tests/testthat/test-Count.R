test_that("Count returns event counts per metadata row and channel", {
  md <- data.frame(
    Experiment = factor(c("test", "test")),
    Step = ordered(c("1", "1"), levels = "1"),
    RecordingID = as.integer(c(1, 1)),
    RunUID = as.integer(c(1, 2))
  )

  x <- EPhysData::newEPhysEvents(
    Data = list(
      list(A = c(0.10, 0.20), B = numeric(0)),
      list(A = 0.30, B = c(0.40, 0.50, 0.60))
    ),
    TimeTrace = numeric(0),
    Metadata = md,
    Channels = c("A", "B"),
    Channel_Metadata = data.frame(row.names = c("A", "B")),
    TimeUnits = "s"
  )

  counts <- Count(x)

  expect_true(is.matrix(counts))
  expect_type(counts, "integer")
  expect_equal(dim(counts), c(2L, 2L))
  expect_equal(rownames(counts), c("1", "2"))
  expect_equal(colnames(counts), c("A", "B"))
  expect_equal(
    unname(counts),
    matrix(c(2L, 0L, 1L, 3L), nrow = 2L, byrow = TRUE)
  )
})

test_that("Count works after time and channel subsetting", {
  md <- data.frame(
    Experiment = factor(c("test", "test")),
    Step = ordered(c("1", "1"), levels = "1"),
    RecordingID = as.integer(c(1, 1)),
    RunUID = as.integer(c(1, 2))
  )

  x <- EPhysData::newEPhysEvents(
    Data = list(
      list(A = c(0.10, 0.20), B = 0.12),
      list(A = c(0.05, 0.30), B = c(0.10, 0.25))
    ),
    TimeTrace = numeric(0),
    Metadata = md,
    Channels = c("A", "B"),
    Channel_Metadata = data.frame(row.names = c("A", "B")),
    TimeUnits = "s"
  )

  x_window <- EPhysData::Subset(
    x,
    Channels = "A",
    TimeRange = c(0, 0.15)
  )

  counts <- Count(x_window)

  expect_equal(dim(counts), c(2L, 1L))
  expect_equal(colnames(counts), "A")
  expect_equal(as.integer(counts[, "A"]), c(1L, 1L))
})
