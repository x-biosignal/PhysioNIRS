test_that("explicit physiology bands retain pass and attenuate stop signals", {
  fs <- 10
  time <- seq.int(0, 400 * fs - 1L) / fs
  pass <- sin(2 * pi * 0.1 * time)
  stop_signal <- sin(2 * pi * 1.0 * time)
  values <- matrix(
    rep(pass + stop_signal, 6L),
    nrow = length(time), ncol = 6L
  )
  x <- make_short_nirs(
    n_time = length(time), fs = fs, values = values
  )
  before <- serialize(x, NULL, version = 3)
  out <- expect_no_warning(physiologyBandpass(
    x,
    band = "mayer",
    range_hz = c(0.07, 0.13)
  ))
  filtered <- short_assay(out, "OD_mayer")[, 1L]
  middle <- 501:(length(time) - 500L)
  pass_gain <- abs(sum(filtered[middle] * pass[middle]) /
                     sum(pass[middle]^2))
  stop_gain <- abs(sum(filtered[middle] * stop_signal[middle]) /
                     sum(stop_signal[middle]^2))
  expect_equal(pass_gain, 1, tolerance = 0.03)
  expect_lt(stop_gain, 0.01)
  expect_identical(serialize(x, NULL, version = 3), before)
  expect_identical(
    tail(PhysioCore::provenance(out)$activity, 1L),
    "physiologyBandpass"
  )
})

test_that("constant sources yield finite numerical zero", {
  x <- make_short_nirs(values = matrix(4, nrow = 100L, ncol = 6L))
  out <- physiologyBandpass(x, range_hz = c(0.07, 0.13))
  expect_identical(short_assay(out, "OD_mayer"), matrix(
    0, nrow = 100L, ncol = 6L,
    dimnames = dimnames(short_assay(x, "OD"))
  ))
})

test_that("small variation on a large offset is not treated as constant", {
  fs <- 10
  time <- seq.int(0, 100 * fs - 1L) / fs
  signal <- 1e6 + 1e-7 * sin(2 * pi * 0.1 * time)
  x <- make_short_nirs(
    n_time = length(time), fs = fs,
    values = matrix(rep(signal, 6L), ncol = 6L)
  )
  out <- physiologyBandpass(x, range_hz = c(0.07, 0.13))
  filtered <- short_assay(out, "OD_mayer")
  expect_true(all(is.finite(filtered)))
  expect_gt(max(abs(filtered)), 0)
})

test_that("default ranges warn once and explicit ranges suppress warning", {
  x <- make_short_nirs()
  state <- get(".nirs_physiology_state", asNamespace("PhysioNIRS"))
  state$default_warning_emitted <- FALSE
  expect_warning(
    physiologyBandpass(x, band = "respiration"),
    "Default adult"
  )
  expect_no_warning(
    physiologyBandpass(
      x, band = "cardiac", output_assay = "OD_cardiac"
    )
  )
  expect_no_warning(
    physiologyBandpass(
      x, band = "mayer", range_hz = c(0.08, 0.12),
      output_assay = "OD_explicit"
    )
  )
})

test_that("band controls and source contracts fail loudly", {
  x <- make_short_nirs()
  expect_error(physiologyBandpass(x, band = "may"), "exactly")
  expect_error(
    physiologyBandpass(x, range_hz = c(0.1, 0.1)),
    "strictly increasing"
  )
  expect_error(
    physiologyBandpass(x, range_hz = matrix(c(0.1, 0.2))),
    "finite numeric"
  )
  expect_error(
    physiologyBandpass(x, range_hz = c(0.1, 5)),
    "Nyquist"
  )
  expect_error(physiologyBandpass(x, order = 0), "order")
  expect_error(
    physiologyBandpass(x, output_assay = "OD"),
    "already exists"
  )

  short <- make_short_nirs(n_time = 10L)
  expect_error(
    physiologyBandpass(short, range_hz = c(0.07, 0.13)),
    "at least"
  )
  irregular <- x
  SummarizedExperiment::rowData(irregular)$time_seconds[[10L]] <- 0.95
  expect_error(physiologyBandpass(irregular), "uniformly")
})
