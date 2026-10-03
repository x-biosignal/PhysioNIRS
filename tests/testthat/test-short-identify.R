test_that("short channels use strict metre geometry and stable identity", {
  x <- make_short_nirs(distances = c(0.008, 0.01, 0.03))
  map <- identifyShortChannels(x)

  expect_s3_class(map, "nirs_short_channels")
  expect_identical(
    names(map),
    c(
      "channel_index", "channel_id", "source_index", "detector_index",
      "distance_m", "midpoint_x_m", "midpoint_y_m", "midpoint_z_m",
      "is_short", "identity_kind"
    )
  )
  expect_identical(map$channel_index, 1:6)
  expect_identical(
    map$channel_id,
    c(
      "S1_D1_wl760", "S1_D1_wl850", "S2_D2_wl760",
      "S2_D2_wl850", "S3_D3_wl760", "S3_D3_wl850"
    )
  )
  expect_equal(map$distance_m, rep(c(0.008, 0.01, 0.03), each = 2))
  expect_identical(map$is_short, rep(c(TRUE, FALSE, FALSE), each = 2))
  expect_true(all(map$midpoint_z_m == 0))
  expect_true(nzchar(attr(map, "source_fingerprint")))
  expect_true(nzchar(attr(map, "probe_fingerprint")))
  expect_identical(attr(map, "threshold_m"), 0.01)
})

test_that("native units and 3-D geometry are governed", {
  x <- make_short_nirs()
  meta <- S4Vectors::metadata(x)
  meta$snirf$probe$sourcePos2D <- meta$snirf$probe$sourcePos2D * 1000
  meta$snirf$probe$detectorPos2D <- meta$snirf$probe$detectorPos2D * 1000
  meta$snirf$probe$LengthUnit <- "mm"
  S4Vectors::metadata(x) <- meta
  expect_equal(
    identifyShortChannels(x)$distance_m,
    rep(c(0.008, 0.009, 0.03), each = 2)
  )

  meta <- S4Vectors::metadata(x)
  meta$snirf$probe$sourcePos3D <- cbind(
    meta$snirf$probe$sourcePos2D, c(0, 0, 0)
  )
  meta$snirf$probe$detectorPos3D <- cbind(
    meta$snirf$probe$detectorPos2D, c(6, 8, 40)
  )
  S4Vectors::metadata(x) <- meta
  map <- identifyShortChannels(x)
  expect_equal(map$distance_m, rep(c(0.01, sqrt(145) / 1000, 0.05), each = 2))
  expect_equal(map$midpoint_z_m, rep(c(0.003, 0.004, 0.02), each = 2))
})

test_that("pair-collapsed maps reconcile stored distance", {
  hb <- make_short_hb()
  map <- identifyShortChannels(hb)
  expect_identical(map$identity_kind, rep("pair", 2L))
  expect_identical(map$channel_id, c("S1_D1", "S2_D2"))
  expect_identical(map$is_short, c(TRUE, FALSE))

  bad <- hb
  SummarizedExperiment::colData(bad)$source_detector_distance_m[[1L]] <- 0.02
  expect_error(
    identifyShortChannels(bad),
    "disagrees with probe coordinates"
  )
})

test_that("maps are source-bound and reject malformed geometry", {
  x <- make_short_nirs()
  map <- identifyShortChannels(x)

  changed <- x
  value <- short_assay(changed, "OD")
  value[1L, 1L] <- 100
  SummarizedExperiment::assay(
    changed, "OD", withDimnames = FALSE
  ) <- value
  expect_error(
    shortSeparationRegress(changed, short = map),
    "stale"
  )
  changed <- x
  meta <- S4Vectors::metadata(changed)
  meta$snirf$probe$detectorPos2D[1L, 2L] <- 0.007
  S4Vectors::metadata(changed) <- meta
  expect_error(
    shortSeparationRegress(changed, short = map),
    "stale"
  )
  changed <- x
  meta <- S4Vectors::metadata(changed)
  meta$nirs$assays$OD$log_convention <- "base10"
  S4Vectors::metadata(changed) <- meta
  expect_error(
    shortSeparationRegress(changed, short = map),
    "stale"
  )
  corrupt <- map
  corrupt$is_short[[1L]] <- FALSE
  expect_error(
    shortSeparationRegress(x, short = corrupt),
    "stale"
  )
  corrupt <- map
  attr(corrupt, "unexpected") <- "not governed"
  expect_error(
    shortSeparationRegress(x, short = corrupt),
    "stale"
  )
  expect_error(identifyShortChannels(x, threshold_m = 0), "threshold_m")
  expect_error(identifyShortChannels(x, threshold_m = matrix(0.01)), "threshold")
})

test_that("malformed measurement identity cannot fall through to pair mode", {
  x <- make_short_nirs()
  columns <- SummarizedExperiment::colData(x)
  columns$channel_id <- paste0("pair_", seq_len(ncol(x)))
  columns$source_detector_distance_m <- rep(0.03, ncol(x))
  SummarizedExperiment::colData(x) <- columns
  meta <- S4Vectors::metadata(x)
  meta$snirf$measurement_list$measurement_index[[1L]] <- 2L
  S4Vectors::metadata(x) <- meta
  expect_error(
    identifyShortChannels(x),
    "measurement_index"
  )
})
