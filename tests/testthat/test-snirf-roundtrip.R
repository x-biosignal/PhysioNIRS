test_that("full and compact SNIRF round trips are lossless", {
  x <- make_snirf_experiment()
  expected <- SummarizedExperiment::assay(x)
  for (compact in c(FALSE, TRUE)) {
    path <- write_fixture_file(x, compact)
    on.exit(unlink(path), add = TRUE)
    y <- readSNIRF(path)
    expect_identical(dim(y), c(17L, 6L))
    expect_lte(max(abs(SummarizedExperiment::assay(y) - expected)), 1e-12)
    expect_equal(
      SummarizedExperiment::rowData(y)$time_seconds,
      SummarizedExperiment::rowData(x)$time_seconds,
      tolerance = 1e-12
    )
    expect_equal(PhysioCore::samplingRate(y), 10, tolerance = 1e-12)
    expect_identical(
      as.list(measurementList(y)),
      as.list(measurementList(x))
    )
    expect_identical(
      S4Vectors::metadata(y)$snirf$metadata_tags$CustomTag,
      "preserved"
    )
  }
})

test_that("irregular time is preserved without an approximate rate", {
  time <- c(0, 0.1, 0.205, 0.31, 0.42)
  x <- make_snirf_experiment(n_time = length(time), time = time)
  path <- write_fixture_file(x)
  on.exit(unlink(path), add = TRUE)
  y <- readSNIRF(path)
  expect_identical(S4Vectors::metadata(y)$snirf$time_sampling, "irregular")
  expect_true(is.na(PhysioCore::samplingRate(y)))
  expect_equal(SummarizedExperiment::rowData(y)$time_seconds, time)
  expect_error(
    writeSNIRF(x, tempfile(fileext = ".snirf"), compact_time = TRUE),
    "uniformly"
  )
})

test_that("time x measurement orientation never depends on relative sizes", {
  x <- make_snirf_experiment(n_time = 3L)
  path <- write_fixture_file(x)
  on.exit(unlink(path), add = TRUE)
  y <- readSNIRF(path)
  expect_identical(dim(y), c(3L, 6L))
  expect_identical(
    unname(SummarizedExperiment::assay(y)),
    unname(SummarizedExperiment::assay(x))
  )
})

test_that("stimulus tables and event semantics round trip", {
  x <- make_snirf_experiment()
  path <- write_fixture_file(x)
  on.exit(unlink(path), add = TRUE)
  y <- readSNIRF(path)
  stim <- S4Vectors::metadata(y)$snirf$stim
  expect_identical(stim[[1L]]$data_labels,
                   c("starttime", "duration", "value", "difficulty"))
  expect_identical(stim[[1L]]$data[[1L, 4L]], 11)
  event_data <- methods::slot(PhysioCore::getEvents(y), "events")
  events <- data.frame(
    onset = event_data$onset,
    duration = event_data$duration,
    type = event_data$type,
    value = event_data$value,
    stringsAsFactors = FALSE
  )
  expect_equal(events$onset, c(0.2, 0.8), tolerance = 1e-12)
  expect_equal(events$duration, c(0.1, 0.2), tolerance = 1e-12)
  expect_identical(events$type, c("left", "right"))
  expect_identical(events$value, c("1", "2"))
})

test_that("edited events replace stale stimulus tables with a warning", {
  x <- make_snirf_experiment()
  x <- PhysioCore::setEvents(
    x, PhysioCore::PhysioEvents(0.3, 0.4, "edited", "3")
  )
  path <- tempfile(fileext = ".snirf")
  on.exit(unlink(path), add = TRUE)
  expect_warning(writeSNIRF(x, path), "extra columns were dropped")
  y <- readSNIRF(path)
  stim <- S4Vectors::metadata(y)$snirf$stim
  expect_length(stim, 1L)
  expect_identical(stim[[1L]]$name, "edited")
  expect_identical(ncol(stim[[1L]]$data), 3L)
})

test_that("writer is collision-safe and restores overwrite backups", {
  x <- make_snirf_experiment()
  path <- write_fixture_file(x)
  on.exit(unlink(path), add = TRUE)
  before <- readBin(path, "raw", n = file.info(path)$size)
  expect_error(writeSNIRF(x, path), "exists")
  old <- options(PhysioNIRS.write_fail_after_backup = TRUE)
  on.exit(options(old), add = TRUE)
  expect_error(writeSNIRF(x, path, overwrite = TRUE), "Injected")
  after <- readBin(path, "raw", n = file.info(path)$size)
  expect_identical(after, before)
  leftovers <- list.files(
    dirname(path), pattern = "^\\.physionirs-(backup-)?"
  )
  expect_length(leftovers, 0L)
})
