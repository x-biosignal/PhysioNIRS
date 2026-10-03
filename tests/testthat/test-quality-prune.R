test_that("mark is default, auditable, and leaves assays unchanged", {
  x <- make_quality_nirs()
  quality <- scalpCouplingIndex(x, window_seconds = 10)
  before <- quality_assay(x)
  marked <- pruneChannels(x, quality)

  expect_identical(quality_assay(marked), before)
  expect_identical(quality_assay(x), before)
  columns <- SummarizedExperiment::colData(marked)
  expect_true(all(c(
    "nirs_quality_pass", "nirs_quality_metrics",
    "nirs_quality_failed_windows", "nirs_quality_fingerprint"
  ) %in% names(columns)))
  expect_identical(
    as.logical(columns$nirs_quality_pass),
    c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_length(S4Vectors::metadata(marked)$nirs$quality, 1L)
  expect_error(pruneChannels(marked, quality), "already exist")
})

test_that("drop propagates decisions to complete wavelength pairs", {
  x <- make_quality_nirs()
  quality <- scalpCouplingIndex(x, window_seconds = 10)
  dropped <- pruneChannels(x, quality, action = "drop")

  expect_equal(ncol(dropped), 2L)
  expect_identical(
    as.character(measurementList(dropped)$channel_label),
    as.character(measurementList(x)$channel_label[1:2])
  )
  expect_identical(
    as.integer(measurementList(dropped)$measurement_index), 1:2
  )
  expect_equal(ncol(x), 4L)
})

test_that("mark preserves channel decisions while drop propagates pairs", {
  fs <- 20
  time <- 0:399 / fs
  good <- 2 * sin(2 * pi * time) + 0.1 * sin(2 * pi * 3 * time)
  bad <- 0.1 * sin(2 * pi * time) + 2 * sin(2 * pi * 3 * time)
  x <- make_cw_nirs(
    wavelengths = c(760, 850),
    distances = c(0.03, 0.04),
    n_time = length(time),
    assay_name = "OD",
    values = cbind(good, bad, good, good)
  )
  SummarizedExperiment::rowData(x)$time_seconds <- time
  methods::slot(x, "samplingRate") <- fs
  quality <- signalQualityIndex(
    x,
    method = "snr",
    window_seconds = 20,
    cardiac_range_hz = c(0.8, 1.2),
    noise_range_hz = matrix(c(2.8, 3.2), nrow = 1L),
    threshold = 0
  )
  marked <- pruneChannels(x, quality, action = "mark")
  dropped <- pruneChannels(x, quality, action = "drop")

  expect_identical(
    as.logical(SummarizedExperiment::colData(marked)$nirs_quality_pass),
    c(TRUE, FALSE, TRUE, TRUE)
  )
  expect_identical(
    as.integer(measurementList(dropped)$source_index),
    c(2L, 2L)
  )
})

test_that("drop rejects all-channel removal and corrected provenance", {
  x <- make_quality_nirs(coupled = c(FALSE))
  quality <- scalpCouplingIndex(x, window_seconds = 10)
  expect_error(
    pruneChannels(x, quality, action = "drop"),
    "remove every channel"
  )

  corrected <- make_quality_nirs()
  metadata <- S4Vectors::metadata(corrected)
  metadata$nirs$assays$OD$motion_corrected <- TRUE
  S4Vectors::metadata(corrected) <- metadata
  quality <- scalpCouplingIndex(corrected, window_seconds = 10)
  expect_error(
    pruneChannels(corrected, quality, action = "drop"),
    "stored correction provenance"
  )
})

test_that("multiple metrics obey require_all without positional recycling", {
  x <- make_quality_nirs()
  sci <- scalpCouplingIndex(x, window_seconds = 10)
  peak <- signalQualityIndex(
    x, method = "peak_power", window_seconds = 10,
    threshold = 0
  )
  all_required <- pruneChannels(
    x, list(sci = sci, peak = peak), require_all = TRUE
  )
  any_required <- pruneChannels(
    x, list(sci = sci, peak = peak), require_all = FALSE
  )

  expect_identical(
    as.logical(SummarizedExperiment::colData(all_required)$nirs_quality_pass),
    c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_true(all(
    SummarizedExperiment::colData(any_required)$nirs_quality_pass
  ))
  expect_error(
    pruneChannels(x, list(a = sci, b = sci)),
    "distinct metric"
  )
  expect_error(pruneChannels(x, list(sci)))
  expect_error(pruneChannels(x, sci, require_all = NA))
  expect_error(pruneChannels(x, sci, action = "dro"))
})

test_that("quality result order and assay contracts remain exact", {
  x <- make_quality_nirs()
  quality <- scalpCouplingIndex(x, window_seconds = 10)
  reordered <- x[, c(3L, 4L, 1L, 2L), drop = FALSE]
  metadata <- S4Vectors::metadata(reordered)
  metadata$snirf$measurement_list <-
    metadata$snirf$measurement_list[c(3L, 4L, 1L, 2L), ,
                                    drop = FALSE]
  metadata$snirf$measurement_list$measurement_index <- 1:4
  S4Vectors::metadata(reordered) <- metadata
  expect_error(pruneChannels(reordered, quality), "another source")

  no_contract <- x
  metadata <- S4Vectors::metadata(no_contract)
  metadata$nirs$assays$OD <- NULL
  S4Vectors::metadata(no_contract) <- metadata
  expect_error(scalpCouplingIndex(no_contract), "complete governed contract")
})
