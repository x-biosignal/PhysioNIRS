test_that("all-channel nuisance design is aligned and standardized", {
  x <- make_short_nirs()
  before <- serialize(x, NULL, version = 3)
  design <- shortSeparationDesign(x)

  expect_s3_class(design, "nirs_nuisance_design")
  expect_true(is.matrix(design))
  expect_identical(dim(design), c(nrow(x), 4L))
  expect_identical(
    colnames(design),
    c(
      "short_S1_D1_wl760", "short_S1_D1_wl850",
      "short_S2_D2_wl760", "short_S2_D2_wl850"
    )
  )
  expect_equal(unname(colMeans(design)), rep(0, 4L), tolerance = 1e-14)
  expect_equal(unname(apply(design, 2L, stats::sd)), rep(1, 4L),
               tolerance = 1e-14)
  expect_identical(attr(design, "time_seconds"),
                   SummarizedExperiment::rowData(x)$time_seconds)
  expect_equal(attr(design, "sampling_rate_hz"), 10, tolerance = 1e-12)
  expect_identical(serialize(x, NULL, version = 3), before)
})

test_that("mean aggregation respects exact compatibility classes", {
  x <- make_short_nirs()
  source <- short_assay(x, "OD")
  design <- shortSeparationDesign(
    x, aggregation = "mean", standardize = FALSE
  )
  expect_identical(
    colnames(design),
    c("short_mean_wl760", "short_mean_wl850")
  )
  expect_equal(design[, 1L], rowMeans(source[, c(1L, 3L)]))
  expect_equal(design[, 2L], rowMeans(source[, c(2L, 4L)]))

  hb <- make_short_hb()
  pair <- shortSeparationDesign(
    hb, assay_name = "HbO", aggregation = "mean", standardize = FALSE
  )
  expect_identical(colnames(pair), "short_mean_HbO")
  expect_equal(pair[, 1L], short_assay(hb, "HbO")[, 1L])
})

test_that("band copies follow requested band then source order", {
  x <- make_short_nirs()
  design <- shortSeparationDesign(
    x,
    aggregation = "mean",
    physiology_bands = c("respiration", "mayer"),
    band_ranges = list(
      respiration = c(0.15, 0.40),
      mayer = c(0.07, 0.13)
    ),
    standardize = FALSE
  )
  expect_identical(
    colnames(design),
    c(
      "short_mean_wl760", "short_mean_wl850",
      "respiration_short_mean_wl760",
      "respiration_short_mean_wl850",
      "mayer_short_mean_wl760", "mayer_short_mean_wl850"
    )
  )
  roles <- attr(design, "column_roles")
  expect_identical(roles$column, colnames(design))
  expect_identical(
    roles$role,
    c("short_mean", "short_mean", rep("physiology_band", 4L))
  )
})

test_that("design output is deterministic and source-bound", {
  x <- make_short_nirs()
  map <- identifyShortChannels(x)
  first <- shortSeparationDesign(
    x, short = map, aggregation = "mean", standardize = FALSE
  )
  second <- shortSeparationDesign(
    x, short = map, aggregation = "mean", standardize = FALSE
  )
  expect_identical(
    serialize(first, NULL, version = 3),
    serialize(second, NULL, version = 3)
  )
  changed <- x
  value <- short_assay(changed, "OD")
  value[1L, 1L] <- value[1L, 1L] + 1
  SummarizedExperiment::assay(
    changed, "OD", withDimnames = FALSE
  ) <- value
  expect_error(
    shortSeparationDesign(changed, short = map),
    "stale"
  )
})

test_that("design controls reject ambiguity and degeneracy", {
  x <- make_short_nirs()
  expect_error(
    shortSeparationDesign(x, aggregation = "m"),
    "exactly"
  )
  expect_error(
    shortSeparationDesign(x, physiology_bands = c("mayer", "mayer")),
    "unique"
  )
  expect_error(
    shortSeparationDesign(
      x, physiology_bands = "mayer",
      band_ranges = list(cardiac = c(0.7, 2))
    ),
    "name every"
  )
  expect_error(
    shortSeparationDesign(
      x, physiology_bands = "may",
      band_ranges = list(may = c(0.07, 0.13))
    ),
    "unique exact"
  )
  expect_error(
    shortSeparationDesign(x, standardize = NA),
    "logical"
  )
  constant <- make_short_nirs(
    values = matrix(1, nrow = 100L, ncol = 6L)
  )
  expect_error(
    shortSeparationDesign(constant),
    "zero or near-zero"
  )
})
