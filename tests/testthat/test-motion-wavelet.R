test_that("periodized db2 coefficients and inverse match the oracle", {
  x <- seq_len(8)
  result <- PhysioNIRS:::.nirs_db2_dwt_periodic(x)
  expect_equal(
    result$approximation,
    c(4.760278777324327, 3.725002596914244, 6.553429721660434,
      10.417133026816707),
    tolerance = 1e-13
  )
  expect_equal(
    result$detail,
    c(-1.035276180410083, 0, 0, 3.863703305156273),
    tolerance = 1e-13
  )
  expect_equal(
    PhysioNIRS:::.nirs_db2_idwt_periodic(
      result$approximation, result$detail
    ),
    x,
    tolerance = 1e-13
  )

  coefficients <- PhysioNIRS:::.nirs_shift_invariant_db2(
    sin(seq_len(128) / 9), 4L
  )
  reconstructed <- PhysioNIRS:::.nirs_inverse_shift_invariant_db2(
    coefficients
  )
  expect_equal(reconstructed, sin(seq_len(128) / 9), tolerance = 1e-12)
})

test_that("wavelet correction is identity when no coefficients are rejected", {
  x <- make_motion_nirs(n_time = 129L)
  out <- waveletMotionCorrect(x, iqr = 1e9)
  expect_equal(
    motion_assay(out, "OD_wavelet"),
    motion_assay(x, "OD"),
    tolerance = 1e-12
  )
  provenance <- S4Vectors::metadata(out)$provenance
  expect_true(length(provenance) >= 1L)
})

test_that("wavelet correction reduces isolated spikes and retains length", {
  time <- seq(0, 25.6, by = 0.1)
  clean <- sin(2 * pi * 0.08 * time)
  spiked <- clean
  spiked[120] <- spiked[120] + 8
  x <- make_motion_nirs(values = cbind(spiked, spiked, clean, clean))
  out <- waveletMotionCorrect(x, iqr = 1.5)
  corrected <- motion_assay(out, "OD_wavelet")
  expect_identical(dim(corrected), dim(motion_assay(x, "OD")))
  expect_lt(abs(corrected[120, 1] - clean[120]), 4)
  expect_true(all(is.finite(corrected)))
})

test_that("wavelet correction validates scale, wavelet, and output", {
  short <- make_motion_nirs(n_time = 32L)
  expect_error(
    waveletMotionCorrect(short, min_level = 4L),
    "too short"
  )
  x <- make_motion_nirs()
  expect_error(waveletMotionCorrect(x, iqr = 0), "`iqr`")
  expect_error(
    waveletMotionCorrect(x, output_assay = "OD"),
    "already exists"
  )
})
