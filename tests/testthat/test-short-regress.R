test_that("nearest regression matches independent centred OLS", {
  x <- make_short_nirs(distances = c(0.008, 0.03, 0.04))
  map <- identifyShortChannels(x)
  source <- short_assay(x, "OD")
  before <- serialize(x, NULL, version = 3)
  out <- shortSeparationRegress(x, short = map, method = "nearest")
  corrected <- short_assay(out, "OD_ssr")

  for (target in 3:6) {
    compatible <- if (target %% 2L) 1L else 2L
    predictor <- source[, compatible]
    beta <- stats::cov(source[, target], predictor) /
      stats::var(predictor)
    expected <- source[, target] - (predictor - mean(predictor)) * beta
    expect_equal(corrected[, target], expected, tolerance = 1e-11)
    expect_equal(mean(corrected[, target]), mean(source[, target]),
                 tolerance = 1e-13)
  }
  expect_identical(corrected[, 1:2], source[, 1:2])
  expect_identical(serialize(x, NULL, version = 3), before)
  expect_true("OD_ssr" %in% SummarizedExperiment::assayNames(out))
  expect_identical(
    S4Vectors::metadata(out)$nirs$assays$OD_ssr$method,
    "short_separation_regression"
  )
  expect_identical(
    tail(PhysioCore::provenance(out)$activity, 1L),
    "shortSeparationRegress"
  )
})

test_that("all-short regression uses deterministic full-rank SVD", {
  x <- make_short_nirs()
  source <- short_assay(x, "OD")
  out <- shortSeparationRegress(x, method = "all")
  corrected <- short_assay(out, "OD_ssr")
  design_760 <- cbind(1, source[, c(1L, 3L)])
  design_850 <- cbind(1, source[, c(2L, 4L)])
  for (target in c(5L, 6L)) {
    design <- if (target == 5L) design_760 else design_850
    beta <- solve(crossprod(design), crossprod(design, source[, target]))
    expected <- source[, target] -
      sweep(design[, -1L, drop = FALSE], 2L,
            colMeans(design[, -1L, drop = FALSE]), "-") %*% beta[-1L]
    expect_equal(corrected[, target], as.numeric(expected), tolerance = 1e-10)
  }

  collinear <- x
  value <- short_assay(collinear, "OD")
  value[, 3L] <- value[, 1L]
  value[, 4L] <- value[, 2L]
  SummarizedExperiment::assay(
    collinear, "OD", withDimnames = FALSE
  ) <- value
  expect_error(
    shortSeparationRegress(collinear, method = "all"),
    "rank-deficient"
  )
})

test_that("lag sign and unchanged edge rows are exact", {
  n <- 200L
  fs <- 10
  short_760 <- sin(seq_len(n) / 7) + seq_len(n) / 1000
  short_850 <- cos(seq_len(n) / 9) + seq_len(n) / 1100
  target_760 <- c(99, short_760[-n])
  target_850 <- c(-99, short_850[-n])
  values <- cbind(short_760, short_850, target_760, target_850)
  x <- make_short_nirs(
    n_time = n, fs = fs, distances = c(0.008, 0.03),
    values = values
  )
  out <- shortSeparationRegress(x, lag_seconds = 0.1)
  corrected <- short_assay(out, "OD_ssr")

  expect_identical(unname(corrected[1L, 3:4]), unname(values[1L, 3:4]))
  expect_lt(stats::sd(corrected[-1L, 3L]), 1e-10)
  expect_lt(stats::sd(corrected[-1L, 4L]), 1e-10)
  expect_error(
    shortSeparationRegress(x, lag_seconds = 0.15),
    "whole sample"
  )
  expect_error(
    shortSeparationRegress(x, lag_seconds = 19.8),
    "n_time - 2"
  )
})

test_that("wavelength compatibility and nearest ties are explicit", {
  x <- make_short_nirs()
  meta <- S4Vectors::metadata(x)
  meta$snirf$measurement_list$wavelength_nm[c(1L, 3L)] <- 780
  S4Vectors::metadata(x) <- meta
  expect_error(
    shortSeparationRegress(x),
    "No compatible short channel"
  )

  x <- make_short_nirs()
  meta <- S4Vectors::metadata(x)
  meta$snirf$probe$sourcePos2D[3L, ] <- c(0.5, 0)
  meta$snirf$probe$detectorPos2D[3L, ] <- c(0.5, 0.03)
  meta$snirf$probe$sourcePos2D[1L, ] <- c(0, 0.0005)
  meta$snirf$probe$detectorPos2D[1L, ] <- c(0, 0.0085)
  meta$snirf$probe$sourcePos2D[2L, ] <- c(1, 0)
  meta$snirf$probe$detectorPos2D[2L, ] <- c(1, 0.009)
  S4Vectors::metadata(x) <- meta
  map <- identifyShortChannels(x)
  nearest <- PhysioNIRS:::.nirs_nearest_short(map, c(1L, 3L), 5L)
  expect_identical(
    nearest$tied,
    c(1L, 3L)
  )
  expect_identical(
    nearest$index,
    1L
  )
})

test_that("pair-collapsed regression preserves assay compatibility", {
  hb <- make_short_hb()
  out <- shortSeparationRegress(hb, assay_name = "HbO")
  expect_true("HbO_ssr" %in% SummarizedExperiment::assayNames(out))
  expect_identical(short_assay(out, "HbO_ssr")[, 1L],
                   short_assay(hb, "HbO")[, 1L])
  expect_identical(short_assay(out, "HbR"), short_assay(hb, "HbR"))
  expect_error(
    shortSeparationRegress(hb, assay_name = "OD"),
    "does not exist"
  )
})

test_that("regression rejects malformed control and collisions", {
  x <- make_short_nirs()
  expect_error(shortSeparationRegress(x, method = "near"), "exactly")
  expect_error(shortSeparationRegress(x, condition_limit = 1), "condition_limit")
  expect_error(
    shortSeparationRegress(x, output_assay = "OD"),
    "already exists"
  )
  expect_error(shortSeparationRegress(x, short = c(TRUE, FALSE)), "complete")
  only_short <- make_short_nirs(distances = c(0.005, 0.006, 0.007))
  expect_error(shortSeparationRegress(only_short), "at least one short")
  no_short <- make_short_nirs(distances = c(0.02, 0.03, 0.04))
  expect_error(shortSeparationRegress(no_short), "at least one short")
})
