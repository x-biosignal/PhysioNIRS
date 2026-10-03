test_that("TDDR is identity on constants and preserves means", {
  constant <- make_motion_nirs(values = matrix(2.5, 128, 4))
  out <- tddr(constant)
  expect_identical(motion_assay(out, "OD_tddr"), matrix(
    2.5, 128, 4,
    dimnames = dimnames(motion_assay(constant, "OD"))
  ))

  x <- make_motion_nirs()
  value <- motion_assay(x, "OD")
  value[55:60, ] <- value[55:60, ] + 2
  SummarizedExperiment::assay(x, "OD", withDimnames = FALSE) <- value
  corrected <- motion_assay(tddr(x), "OD_tddr")
  expect_equal(colMeans(corrected), colMeans(value), tolerance = 1e-13)
  expect_true(all(is.finite(corrected)))
})

test_that("TDDR matches the pinned reference up to its DC convention", {
  time <- seq(0, 12.7, by = 0.1)
  source <- sin(2 * pi * 0.1 * time)
  source[40:45] <- source[40:45] + 2
  result <- PhysioNIRS:::.nirs_tddr_channel(
    source, fs = 10, cutoff_hz = 0.5, tune = 4.685, max_iter = 50L
  )$value
  pinned_reference <- c(
    0.22029013111501081, 0.26069931347728714,
    0.30071915642678609, 0.33999729647761467,
    0.37820661879393913, 0.41507693245278232,
    0.45040683909522811, 0.48405500199893170
  )
  expect_equal(
    result[seq_along(pinned_reference)] - mean(result),
    pinned_reference - 0.23908607435190424,
    tolerance = 1e-9
  )
  expect_equal(mean(result), mean(source), tolerance = 1e-14)
})

test_that("TDDR suppresses steps while retaining clean haemodynamics", {
  time <- seq(0, 99.9, by = 0.1)
  clean <- sin(2 * pi * 0.05 * time)
  motion <- clean
  motion[120:length(motion)] <- motion[120:length(motion)] + 3
  values <- cbind(motion, motion, motion, motion)
  x <- make_motion_nirs(values = values)
  corrected <- motion_assay(tddr(x), "OD_tddr")[, 1]
  before_jump <- abs(mean(motion[120:150]) - mean(motion[90:119]))
  after_jump <- abs(mean(corrected[120:150]) - mean(corrected[90:119]))
  expect_lt(after_jump, 0.5 * before_jump)

  clean_x <- make_motion_nirs(values = cbind(clean, clean, clean, clean))
  clean_corrected <- motion_assay(tddr(clean_x), "OD_tddr")[, 1]
  amplitude_ratio <- stats::sd(clean_corrected) / stats::sd(clean)
  expect_lt(abs(amplitude_ratio - 1), 0.1)
  expect_gt(stats::cor(clean_corrected, clean), 0.99)
})

test_that("TDDR validates cutoff and deterministic output", {
  x <- make_motion_nirs()
  one <- tddr(x)
  two <- tddr(x)
  expect_identical(
    motion_assay(one, "OD_tddr"),
    motion_assay(two, "OD_tddr")
  )
  expect_error(tddr(x, cutoff_hz = 5.1), "Nyquist")
  expect_silent(tddr(x, cutoff_hz = 5))
  expect_error(tddr(x, tune = 0), "`tune`")
})
