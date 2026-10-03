test_that("PHOEBE-style peak power copies conservative pair scores", {
  x <- make_quality_nirs()
  result <- signalQualityIndex(
    x, method = "peak_power", window_seconds = 10
  )

  expect_identical(result$metric, "phoebe_peak_power")
  expect_identical(result$score[1L, ], result$score[2L, ])
  expect_identical(result$score[3L, ], result$score[4L, ])
  expect_true(mean(result$score[1L, ]) > mean(result$score[3L, ]))
  expect_identical(result$threshold, 0.1)
  expect_match(
    result$implementation$periodogram,
    "Hamming"
  )
})

test_that("peak power supports three and four wavelength groups", {
  for (wavelengths in list(c(730, 760, 850), c(730, 760, 810, 850))) {
    x <- make_quality_nirs(
      wavelengths = wavelengths, coupled = c(TRUE, FALSE)
    )
    result <- signalQualityIndex(
      x, method = "peak_power", window_seconds = 10
    )
    n <- length(wavelengths)
    expect_equal(length(unique(result$score[seq_len(n), 1L])), 1L)
    expect_equal(
      length(unique(result$score[n + seq_len(n), 1L])), 1L
    )
  }
})

test_that("SNR follows an independent Hamming periodogram ratio", {
  fs <- 20
  n <- 400L
  time <- seq.int(0, n - 1L) / fs
  signal <- 2 * sin(2 * pi * 1 * time) +
    0.2 * sin(2 * pi * 3 * time)
  values <- cbind(signal, signal)
  x <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03,
    n_time = n, assay_name = "OD", values = values
  )
  SummarizedExperiment::rowData(x)$time_seconds <- time
  methods::slot(x, "samplingRate") <- fs
  result <- signalQualityIndex(
    x,
    method = "snr",
    window_seconds = 20,
    cardiac_range_hz = c(0.8, 1.2),
    noise_range_hz = matrix(c(2.8, 3.2), nrow = 1L)
  )

  expect_equal(
    unname(result$score[, 1L]), rep(20, 2L), tolerance = 0.15
  )
  expect_identical(result$metric, "snr_db")
  expect_identical(result$threshold, 0)
})

test_that("SNR accepts pair-level haemoglobin identity", {
  x <- make_short_hb(n_time = 400L, fs = 10)
  result <- signalQualityIndex(
    x,
    assay_name = "HbO",
    method = "snr",
    window_seconds = 20,
    cardiac_range_hz = c(0.05, 0.15),
    noise_range_hz = matrix(c(0.3, 1), nrow = 1L)
  )

  expect_identical(result$identity_kind, "pair")
  expect_true(all(is.na(result$wavelength_nm)))
  expect_true(all(is.finite(result$score)))
})

test_that("signal quality rejects invalid methods and spectral bands", {
  x <- make_quality_nirs(n_time = 400L, coupled = c(TRUE))
  expect_error(
    signalQualityIndex(x, method = "peak", window_seconds = 10)
  )
  expect_error(
    signalQualityIndex(
      x, method = "snr", window_seconds = 10,
      noise_range_hz = c(2, 3)
    ),
    "two-column matrix"
  )
  expect_error(
    signalQualityIndex(
      x, method = "snr", window_seconds = 10,
      noise_range_hz = matrix(c(1, 2), nrow = 1L)
    ),
    "must not overlap"
  )
  expect_error(
    signalQualityIndex(
      x, method = "snr", window_seconds = 10,
      noise_range_hz = matrix(c(2, 3, 2.5, 4), ncol = 2L,
                              byrow = TRUE)
    ),
    "disjoint"
  )
  expect_error(
    signalQualityIndex(
      x, method = "peak_power", window_seconds = 10,
      cardiac_range_hz = c(0.7, 5)
    ),
    "below Nyquist"
  )
})

test_that("zero power and zero variance fail conservatively", {
  x <- make_quality_nirs(n_time = 400L, coupled = c(TRUE))
  value <- quality_assay(x)
  value[] <- 1
  SummarizedExperiment::assay(x, "OD", withDimnames = FALSE) <- value

  peak <- signalQualityIndex(
    x, method = "peak_power", window_seconds = 10
  )
  expect_identical(unname(peak$score), matrix(0, 2L, 4L))
  expect_error(
    signalQualityIndex(
      x, method = "snr", window_seconds = 10,
      noise_range_hz = matrix(c(2, 4), nrow = 1L)
    ),
    "strictly positive"
  )
})
