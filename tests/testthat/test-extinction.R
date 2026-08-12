test_that("packaged extinction table has a verified governed schema", {
  path <- system.file(
    "extdata", "haemoglobin-extinction.csv", package = "PhysioNIRS"
  )
  hash_path <- system.file(
    "extdata", "haemoglobin-extinction.sha256", package = "PhysioNIRS"
  )
  table <- utils::read.csv(path, stringsAsFactors = FALSE)
  expected_hash <- strsplit(
    readLines(hash_path, warn = FALSE)[[1L]], "[[:space:]]+"
  )[[1L]][[1L]]
  actual_hash <- digest::digest(path, algo = "sha256", file = TRUE)
  expect_identical(actual_hash, expected_hash)
  expect_identical(
    names(table),
    c(
      "wavelength_nm", "hbo2_cm1_M1", "hbr_cm1_M1",
      "source", "source_version"
    )
  )
  expect_identical(nrow(table), 376L)
  expect_equal(range(table$wavelength_nm), c(250, 1000))
  expect_true(all(diff(table$wavelength_nm) == 2))
  expect_false(anyNA(table))
  expect_true(all(table$hbo2_cm1_M1 >= 0))
  expect_true(all(table$hbr_cm1_M1 >= 0))
  expect_true(file.exists(system.file(
    "extdata", "LICENSE-mne-extinction.txt", package = "PhysioNIRS"
  )))
})

test_that("extinction lookup pins grid, interpolation, units, and order", {
  grid <- extinctionCoefficients(
    c(690, 760, 850), unit = "cm-1 M-1"
  )
  expect_identical(colnames(grid), c("HbO", "HbR"))
  expect_equal(unname(grid[1L, ]), c(276, 2051.96), tolerance = 0)
  expect_equal(unname(grid[2L, ]), c(586, 1548.52), tolerance = 0)
  expect_equal(unname(grid[3L, ]), c(1058, 691.32), tolerance = 0)

  interpolated <- extinctionCoefficients(761, unit = "cm-1 M-1")
  expect_equal(
    unname(interpolated[1L, ]),
    (c(586, 1548.52) + c(598, 1508.44)) / 2,
    tolerance = 1e-12
  )
  natural <- extinctionCoefficients(c(690, 760, 850))
  expect_equal(
    as.vector(natural),
    as.vector(grid) * (100 * log(10)),
    tolerance = 1e-12
  )
  expect_identical(attr(natural, "unit"), "m-1 M-1")
  expect_identical(attr(natural, "wavelength_nm"), c(690, 760, 850))
})

test_that("extinction domain and extrapolation are explicit", {
  expect_error(extinctionCoefficients(249), "outside")
  expect_warning(
    below <- extinctionCoefficients(
      249, unit = "cm-1 M-1", extrapolate = TRUE
    ),
    "extrapolated"
  )
  expected <- c(106112, 112736) -
    (c(105552, 112736) - c(106112, 112736)) / 2
  expect_equal(unname(below[1L, ]), expected, tolerance = 1e-12)
  expect_true(attr(below, "extrapolated"))
  expect_error(extinctionCoefficients(760, unit = "m"), "exactly")
  expect_error(extinctionCoefficients(760, extrapolate = NA), "non-missing")
  expect_error(extinctionCoefficients(c(760, NA)), "finite")
  expect_error(
    suppressWarnings(extinctionCoefficients(1e308, extrapolate = TRUE)),
    "non-finite or negative"
  )
  expect_error(
    extinctionCoefficients(matrix(760, nrow = 1L)),
    "finite numeric"
  )
})
