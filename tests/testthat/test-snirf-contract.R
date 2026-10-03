test_that("path, index, flag, and enum validation is exact", {
  expect_error(readSNIRF("missing.SNIRF"), "\\.snirf")
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  expect_error(readSNIRF(path, nirs_index = 0), "positive whole")
  expect_error(readSNIRF(path, data_index = 1.2), "positive whole")
  expect_error(writeSNIRF(make_snirf_experiment(), tempfile(fileext = ".h5")),
               "exactly")
  expect_error(writeSNIRF(make_snirf_experiment(),
                          tempfile(fileext = ".snirf"), overwrite = NA),
               "non-missing")
  expect_error(sourceDetectorDistances(readSNIRF(path), dimension = "3"),
               "exactly")
  expect_error(
    sourceDetectorDistances(readSNIRF(path), dimension = c("auto", "2d")),
    "exactly"
  )
  expect_error(sourceDetectorDistances(readSNIRF(path), unit = "metres"),
               "exactly")
})

test_that("unindexed and indexed roots cannot be ambiguous", {
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  rhdf5::h5createGroup(path, "/nirs1")
  expect_error(readSNIRF(path), "Conflicting SNIRF roots.*nirs.*nirs1")
})

test_that("an indexed root is selected only at its exact contiguous index", {
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  handle <- rhdf5::H5Fopen(path)
  expect_true(rhdf5::H5Ocopy(handle, "nirs", handle, "nirs1"))
  rhdf5::H5Fclose(handle)
  rhdf5::h5delete(path, "/nirs")
  expect_identical(dim(readSNIRF(path, nirs_index = 1L)), c(17L, 6L))
  expect_error(readSNIRF(path, nirs_index = 2L), "does not exist")

  skipped <- write_fixture_file()
  on.exit(unlink(skipped), add = TRUE)
  handle <- rhdf5::H5Fopen(skipped)
  expect_true(rhdf5::H5Ocopy(handle, "nirs", handle, "nirs2"))
  rhdf5::H5Fclose(handle)
  rhdf5::h5delete(skipped, "/nirs")
  expect_error(readSNIRF(skipped), "contiguous")
})

test_that("required strings use scalar dataspaces", {
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  tree <- rhdf5::h5ls(path, all = TRUE)
  required <- c(
    "formatVersion", "SubjectID", "MeasurementDate", "MeasurementTime",
    "LengthUnit", "TimeUnit", "FrequencyUnit"
  )
  expect_true(all(tree$rank[tree$name %in% required] == 0L))
  rhdf5::h5delete(path, "/formatVersion")
  PhysioNIRS:::.snirf_write_vector(path, "/formatVersion", "1.0")
  expect_error(readSNIRF(path), "scalar HDF5 dataspace")
})

test_that("compact time and unknown units are rejected when malformed", {
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  rhdf5::h5delete(path, "/nirs/data1/time")
  PhysioNIRS:::.snirf_write_vector(
    path, "/nirs/data1/time", c(0, -0.1)
  )
  expect_error(readSNIRF(path), "spacing must be positive")

  path2 <- write_fixture_file()
  on.exit(unlink(path2), add = TRUE)
  rhdf5::h5write("minutes", path2, "/nirs/metaDataTags/TimeUnit")
  expect_error(readSNIRF(path2), "Unsupported SNIRF TimeUnit")
})

test_that("milliseconds convert exactly for time and stimulus fields", {
  path <- write_fixture_file()
  on.exit(unlink(path), add = TRUE)
  time <- rhdf5::h5read(path, "/nirs/data1/time")
  rhdf5::h5write(time * 1000, path, "/nirs/data1/time")
  rhdf5::h5write("ms", path, "/nirs/metaDataTags/TimeUnit")
  for (i in 1:2) {
    dataset <- paste0("/nirs/stim", i, "/data")
    data <- rhdf5::h5read(path, dataset)
    data[1:2, ] <- data[1:2, , drop = FALSE] * 1000
    rhdf5::h5write(data, path, dataset)
  }
  x <- readSNIRF(path)
  expect_equal(
    SummarizedExperiment::rowData(x)$time_seconds,
    seq.int(0, 16L) / 10,
    tolerance = 1e-12
  )
  events <- methods::slot(PhysioCore::getEvents(x), "events")
  expect_equal(events$onset, c(0.2, 0.8), tolerance = 1e-12)
  expect_equal(events$duration, c(0.1, 0.2), tolerance = 1e-12)
})

test_that("indexed and vectorized measurement lists normalize identically", {
  indexed <- write_fixture_file()
  vectorized <- write_fixture_file()
  on.exit(unlink(c(indexed, vectorized)), add = TRUE)
  expected <- measurementList(readSNIRF(indexed))
  convert_to_vectorized(vectorized)
  actual <- measurementList(readSNIRF(vectorized))
  expect_identical(as.list(actual), as.list(expected))
  expect_identical(
    S4Vectors::metadata(readSNIRF(vectorized))$snirf$measurement_encoding,
    "vectorized"
  )
})

test_that("mixed, incomplete, and out-of-range measurement lists fail", {
  mixed <- write_fixture_file()
  on.exit(unlink(mixed), add = TRUE)
  rhdf5::h5createGroup(mixed, "/nirs/data1/measurementLists")
  expect_error(readSNIRF(mixed), "Mixed indexed and vectorized")

  missing <- write_fixture_file()
  on.exit(unlink(missing), add = TRUE)
  rhdf5::h5delete(missing, "/nirs/data1/measurementList6")
  expect_error(readSNIRF(missing), "length does not match")

  outside <- write_fixture_file()
  on.exit(unlink(outside), add = TRUE)
  rhdf5::h5write(3L, outside,
                 "/nirs/data1/measurementList1/sourceIndex")
  expect_error(readSNIRF(outside), "outside probe geometry")
})

test_that("unknown data types and duplicate triplets remain ordered and unique", {
  x <- make_snirf_experiment(duplicate = TRUE)
  path <- write_fixture_file(x)
  on.exit(unlink(path), add = TRUE)
  table <- measurementList(readSNIRF(path))
  expect_identical(table$measurement_index, seq_len(6L))
  expect_identical(table$data_type, c(1L, 1L, 1L, 1L, 99999L, 101L))
  expect_identical(anyDuplicated(table$channel_label), 0L)
})

test_that("measurementList is defensive and checks colData identity", {
  x <- readSNIRF(write_fixture_file())
  copy <- measurementList(x)
  copy$source_index[[1L]] <- 999L
  expect_identical(measurementList(x)$source_index[[1L]], 1L)
  SummarizedExperiment::colData(x)$label[[1L]] <- "changed"
  expect_error(measurementList(x), "disagree")
})

test_that("measurementList rejects incomplete or duplicate labels", {
  x <- make_snirf_experiment()
  table <- S4Vectors::metadata(x)$snirf$measurement_list
  table$channel_label[[2L]] <- table$channel_label[[1L]]
  S4Vectors::metadata(x)$snirf$measurement_list <- table
  SummarizedExperiment::colData(x)$label[[2L]] <- table$channel_label[[2L]]
  expect_error(measurementList(x), "complete and unique")
})

test_that("writer rejects unsupported metadata before creating output", {
  x <- make_snirf_experiment()
  S4Vectors::metadata(x)$snirf$metadata_tags$Nested <- list(value = 1)
  path <- tempfile(fileext = ".snirf")
  expect_error(writeSNIRF(x, path), "finite atomic")
  expect_false(file.exists(path))
})

test_that("writer revalidates measurement identity against the probe", {
  x <- make_snirf_experiment()
  table <- S4Vectors::metadata(x)$snirf$measurement_list
  table$source_index[[1L]] <- -1L
  S4Vectors::metadata(x)$snirf$measurement_list <- table
  path <- tempfile(fileext = ".snirf")
  expect_error(writeSNIRF(x, path), "positive whole")
  expect_false(file.exists(path))

  x <- make_snirf_experiment()
  table <- S4Vectors::metadata(x)$snirf$measurement_list
  table$wavelength_nm[[1L]] <- 999
  S4Vectors::metadata(x)$snirf$measurement_list <- table
  path <- tempfile(fileext = ".snirf")
  expect_error(writeSNIRF(x, path), "probe lookup")
  expect_false(file.exists(path))
})

test_that("non-finite and rank-wrong data are rejected", {
  nonfinite <- write_fixture_file()
  on.exit(unlink(nonfinite), add = TRUE)
  data <- rhdf5::h5read(nonfinite, "/nirs/data1/dataTimeSeries")
  data[[1L, 1L]] <- Inf
  rhdf5::h5write(data, nonfinite, "/nirs/data1/dataTimeSeries")
  expect_error(readSNIRF(nonfinite), "finite")

  wrong_rank <- write_fixture_file()
  on.exit(unlink(wrong_rank), add = TRUE)
  rhdf5::h5delete(wrong_rank, "/nirs/data1/dataTimeSeries")
  PhysioNIRS:::.snirf_write_vector(
    wrong_rank, "/nirs/data1/dataTimeSeries", seq_len(17L * 6L)
  )
  expect_error(readSNIRF(wrong_rank), "rank-2")
})
