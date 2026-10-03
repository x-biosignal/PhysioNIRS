test_that("2-D, 3-D, and dual geometries select exact dimensions", {
  for (geometry in list("2d", "3d", c("2d", "3d"))) {
    x <- make_snirf_experiment(geometry = geometry)
    path <- write_fixture_file(x)
    on.exit(unlink(path), add = TRUE)
    y <- readSNIRF(path)
    expected <- if ("3d" %in% geometry) "3d" else "2d"
    distances <- sourceDetectorDistances(y, dimension = "auto",
                                         unit = "native")
    expect_true(all(distances$geometry_dimension == expected))
    expect_true(all(distances$distance_unit == "cm"))
    expect_equal(distances$distance[c(1L, 3L)], c(4, 5),
                 tolerance = 1e-12)
  }
})

test_that("meter conversion is exact and does not relabel 2-D geometry", {
  path <- write_fixture_file(
    make_snirf_experiment(geometry = "2d")
  )
  on.exit(unlink(path), add = TRUE)
  x <- readSNIRF(path)
  distances <- sourceDetectorDistances(x, unit = "m")
  expect_equal(distances$distance[c(1L, 3L)], c(0.04, 0.05),
               tolerance = 1e-12)
  expect_true(all(distances$geometry_dimension == "2d"))
  expect_error(sourceDetectorDistances(x, dimension = "3d"),
               "not available")
})

test_that("unknown length units are rejected only for meter conversion", {
  x <- make_snirf_experiment()
  S4Vectors::metadata(x)$snirf$probe$LengthUnit <- "furlong"
  S4Vectors::metadata(x)$snirf$metadata_tags$LengthUnit <- "furlong"
  native <- sourceDetectorDistances(x, unit = "native")
  expect_true(all(native$distance_unit == "furlong"))
  expect_error(sourceDetectorDistances(x, unit = "m"),
               "Unsupported SNIRF LengthUnit")
})

test_that("partial and undersized geometry is rejected", {
  x <- make_snirf_experiment(geometry = "2d")
  S4Vectors::metadata(x)$snirf$probe$detectorPos2D <- NULL
  expect_error(sourceDetectorDistances(x), "complete pair")

  x <- make_snirf_experiment(geometry = "3d")
  S4Vectors::metadata(x)$snirf$probe$sourcePos3D <-
    S4Vectors::metadata(x)$snirf$probe$sourcePos3D[1L, , drop = FALSE]
  expect_error(sourceDetectorDistances(x), "does not cover")
})

test_that("landmark geometry and labels have exact governed dimensions", {
  x <- make_snirf_experiment()
  S4Vectors::metadata(x)$snirf$probe$landmarkPos3D <-
    S4Vectors::metadata(x)$snirf$probe$landmarkPos3D[, 1:2, drop = FALSE]
  expect_error(sourceDetectorDistances(x), "landmark probe geometry")

  x <- make_snirf_experiment()
  S4Vectors::metadata(x)$snirf$probe$landmarkLabels <- "Nz"
  expect_error(sourceDetectorDistances(x), "Landmark labels")
})

test_that("package-owned and official fixtures are readable offline", {
  owned <- system.file(
    "extdata", "snirf_reference.snirf", package = "PhysioNIRS"
  )
  official <- system.file(
    "extdata", "snirf_official_simple_probe.snirf",
    package = "PhysioNIRS"
  )
  expect_true(file.exists(owned))
  expect_true(file.exists(official))
  x <- readSNIRF(owned)
  expect_identical(dim(x), c(17L, 6L))
  y <- readSNIRF(official)
  expect_identical(dim(y), c(1200L, 8L))
  expect_identical(sort(unique(measurementList(y)$wavelength_nm)),
                   c(690, 830))
  expect_identical(
    vapply(S4Vectors::metadata(y)$snirf$stim,
           function(entry) nrow(entry$data), integer(1)),
    c(2L, 1L, 1L)
  )
})
