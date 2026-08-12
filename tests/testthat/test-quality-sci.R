test_that("SCI returns the fixed source-bound schema and separates coupling", {
  x <- make_quality_nirs()
  result <- scalpCouplingIndex(
    x, window_seconds = 10, step_seconds = 5
  )

  expect_s3_class(result, "nirs_quality")
  expect_identical(
    names(result),
    c(
      "schema_version", "metric", "assay_name", "identity_kind",
      "channel_id", "source_index", "detector_index", "wavelength_nm",
      "window_start_sample", "window_end_sample", "window_start_time",
      "window_end_time", "score", "pass", "threshold", "parameters",
      "source_fingerprint", "identity_fingerprint", "implementation"
    )
  )
  expect_identical(result$metric, "scalp_coupling_index")
  expect_identical(result$identity_kind, "measurement")
  expect_equal(dim(result$score), c(4L, 11L))
  expect_identical(result$score[1L, ], result$score[2L, ])
  expect_identical(result$score[3L, ], result$score[4L, ])
  expect_true(all(result$pass[1:2, ]))
  expect_false(any(result$pass[3:4, ]))
  expect_match(capture.output(print(result)), "4 channels x 11 windows")
})

test_that("SCI handles multiple wavelengths and exact complete windows", {
  x <- make_quality_nirs(
    n_time = 400L, wavelengths = c(730, 760, 810, 850),
    coupled = c(TRUE, FALSE)
  )
  result <- scalpCouplingIndex(x, window_seconds = NULL)

  expect_equal(dim(result$score), c(8L, 1L))
  expect_identical(result$window_start_sample, 1L)
  expect_identical(result$window_end_sample, 400L)
  expect_equal(length(unique(result$score[1:4, 1L])), 1L)
  expect_equal(length(unique(result$score[5:8, 1L])), 1L)
})

test_that("SCI is invariant to offsets and positive scale", {
  base <- make_quality_nirs(coupled = c(TRUE))
  changed <- make_quality_nirs(
    coupled = c(TRUE), offset = 1234, scale = 7.5
  )
  score1 <- scalpCouplingIndex(base, window_seconds = 10)$score
  score2 <- scalpCouplingIndex(changed, window_seconds = 10)$score

  expect_equal(score1, score2, tolerance = 1e-9)
})

test_that("SCI assigns zero to zero-variance pair members", {
  x <- make_quality_nirs(coupled = c(TRUE))
  value <- quality_assay(x)
  value[, 2L] <- 3
  SummarizedExperiment::assay(x, "OD", withDimnames = FALSE) <- value
  result <- scalpCouplingIndex(x, window_seconds = 10)

  expect_identical(unname(result$score), matrix(0, 2L, 6L))
  expect_false(any(result$pass))
})

test_that("SCI threshold equality passes", {
  x <- make_quality_nirs(coupled = c(TRUE))
  initial <- scalpCouplingIndex(x, window_seconds = NULL)
  equality <- scalpCouplingIndex(
    x, window_seconds = NULL, threshold = initial$score[[1L]]
  )

  expect_true(all(equality$pass))
})

test_that("SCI rejects malformed windows, bands, time, and values", {
  x <- make_quality_nirs(n_time = 100L, coupled = c(TRUE))
  expect_error(
    scalpCouplingIndex(x, window_seconds = 1.05),
    "exact positive sample"
  )
  expect_error(
    scalpCouplingIndex(x, window_seconds = NULL, step_seconds = 1),
    "must be NULL"
  )
  expect_error(scalpCouplingIndex(x, l_freq = 1.5, h_freq = 1.5))
  expect_error(scalpCouplingIndex(x, h_freq = 5))
  expect_error(scalpCouplingIndex(x, threshold = 1.1))
  expect_error(scalpCouplingIndex(x, order = 0))

  bad_time <- x
  SummarizedExperiment::rowData(bad_time)$time_seconds[[5L]] <- 0.405
  expect_error(scalpCouplingIndex(bad_time), "uniform")

  bad_value <- x
  value <- quality_assay(bad_value)
  value[1L, 1L] <- Inf
  SummarizedExperiment::assay(
    bad_value, "OD", withDimnames = FALSE
  ) <- value
  expect_error(scalpCouplingIndex(bad_value), "finite")
})

test_that("SCI rejects incomplete and duplicate wavelength identity", {
  x <- make_quality_nirs(coupled = c(TRUE))
  single <- x[, 1L, drop = FALSE]
  metadata <- S4Vectors::metadata(single)
  metadata$snirf$measurement_list <-
    metadata$snirf$measurement_list[1L, , drop = FALSE]
  metadata$snirf$measurement_list$measurement_index <- 1L
  S4Vectors::metadata(single) <- metadata
  expect_error(scalpCouplingIndex(single), "at least two")

  duplicate <- x
  metadata <- S4Vectors::metadata(duplicate)
  metadata$snirf$measurement_list$wavelength_nm[[2L]] <-
    metadata$snirf$measurement_list$wavelength_nm[[1L]]
  S4Vectors::metadata(duplicate) <- metadata
  expect_error(scalpCouplingIndex(duplicate), "duplicate|distinct")
})

test_that("quality fingerprints reject score and source mutations", {
  x <- make_quality_nirs(coupled = c(TRUE))
  result <- scalpCouplingIndex(x, window_seconds = 10)

  changed <- result
  changed$score[[1L]] <- changed$score[[1L]] - 0.1
  expect_error(pruneChannels(x, changed), "stale, malformed")

  changed <- result
  changed$pass[[1L]] <- !changed$pass[[1L]]
  attr(changed, "quality_fingerprint") <-
    PhysioNIRS:::.nirs_quality_fingerprint(changed)
  expect_error(pruneChannels(x, changed), "stale, malformed")

  changed <- result
  changed$window_start_time[[1L]] <-
    changed$window_start_time[[1L]] + 0.01
  attr(changed, "quality_fingerprint") <-
    PhysioNIRS:::.nirs_quality_fingerprint(changed)
  expect_error(pruneChannels(x, changed), "stale, malformed")

  edited <- x
  value <- quality_assay(edited)
  value[1L, 1L] <- value[1L, 1L] + 1e-6
  SummarizedExperiment::assay(
    edited, "OD", withDimnames = FALSE
  ) <- value
  expect_error(pruneChannels(edited, result), "another source")
})
