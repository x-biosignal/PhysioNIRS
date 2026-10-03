test_that("intensityToOD matches explicit hand calculations", {
  intensity <- rbind(
    c(100, 200),
    c(50, 100),
    c(25, 50),
    c(100, 50)
  )
  x <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03,
    values = intensity
  )
  explicit <- intensityToOD(x, reference = c(100, 200))
  expect_equal(
    unname(SummarizedExperiment::assay(explicit, "OD")),
    unname(-log(sweep(intensity, 2L, c(100, 200), "/"))),
    tolerance = 1e-15
  )

  mean_result <- intensityToOD(x, baseline = "mean")
  expect_equal(
    unname(SummarizedExperiment::assay(mean_result, "OD")),
    unname(-log(sweep(intensity, 2L, colMeans(intensity), "/"))),
    tolerance = 1e-15
  )
  median_result <- intensityToOD(
    x, baseline = "median", baseline_interval = c(0.1, 0.2)
  )
  expected_reference <- apply(intensity[2:3, , drop = FALSE], 2, median)
  expect_equal(
    unname(SummarizedExperiment::assay(median_result, "OD")),
    unname(-log(sweep(intensity, 2L, expected_reference, "/"))),
    tolerance = 1e-15
  )
})

test_that("optical density is scale invariant and sign preserving", {
  intensity <- rbind(c(100, 50), c(80, 25), c(120, 75))
  x <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03,
    values = intensity
  )
  result <- intensityToOD(x, reference = c(100, 50))
  scaled <- x
  SummarizedExperiment::assay(scaled, "raw")[, 1L] <-
    SummarizedExperiment::assay(scaled, "raw")[, 1L] * 17
  scaled_result <- intensityToOD(scaled, reference = c(1700, 50))
  expect_equal(
    SummarizedExperiment::assay(result, "OD"),
    SummarizedExperiment::assay(scaled_result, "OD"),
    tolerance = 1e-14
  )
  expect_gt(SummarizedExperiment::assay(result, "OD")[2L, 1L], 0)
  expect_lt(SummarizedExperiment::assay(result, "OD")[3L, 1L], 0)
  expect_equal(
    unname(SummarizedExperiment::assay(result, "OD")[1L, ]),
    c(0, 0)
  )

  extreme <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03,
    values = rbind(c(1e308, 1e-308), c(1e-308, 1e308))
  )
  extreme_result <- intensityToOD(
    extreme, reference = c(1e308, 1e-308)
  )
  expect_true(all(is.finite(
    SummarizedExperiment::assay(extreme_result, "OD")
  )))
  expect_equal(
    unname(SummarizedExperiment::assay(extreme_result, "OD")[1L, ]),
    c(0, 0)
  )
})

test_that("intensityToOD preserves source and governed object state", {
  x <- make_cw_nirs()
  x <- PhysioCore::logStep(x, "prior", params = list(id = 1))
  source_before <- serialize(x, NULL)
  assays_before <- lapply(
    SummarizedExperiment::assays(x),
    function(value) serialize(value, NULL)
  )
  result <- intensityToOD(x, reference = 1000)
  expect_identical(serialize(x, NULL), source_before)
  expect_identical(
    lapply(SummarizedExperiment::assays(x), function(value) serialize(value, NULL)),
    assays_before
  )
  expect_identical(
    S4Vectors::metadata(result)$events,
    S4Vectors::metadata(x)$events
  )
  expect_identical(
    S4Vectors::metadata(result)$snirf,
    S4Vectors::metadata(x)$snirf
  )
  expect_identical(
    tail(PhysioCore::provenance(result)$activity, 1L),
    "intensityToOD"
  )
  expect_identical(
    S4Vectors::metadata(result)$nirs$assays$OD$log_convention,
    "natural"
  )
})

test_that("non-positive replacement is explicit and local", {
  intensity <- rbind(c(100, 200), c(0, -1), c(50, 100))
  x <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03,
    values = intensity
  )
  expect_error(intensityToOD(x, reference = 100), "row 2.*measurement 1")
  expect_error(
    intensityToOD(x, reference = 100, nonpositive = "replace"),
    "replacement"
  )
  expect_warning(
    result <- intensityToOD(
      x, reference = 100, nonpositive = "replace", replacement = 0.25
    ),
    "Replaced 2"
  )
  expected <- intensity
  expected[expected <= 0] <- 0.25
  expect_equal(
    unname(SummarizedExperiment::assay(result, "OD")),
    unname(-log(expected / 100)),
    tolerance = 1e-15
  )
})

test_that("intensityToOD rejects ambiguous and malformed contracts", {
  x <- make_cw_nirs()
  SummarizedExperiment::assays(x)[["second"]] <-
    SummarizedExperiment::assay(x, "raw")
  expect_error(intensityToOD(x), "exactly one eligible")
  expect_error(intensityToOD(x, assay_name = "missing"), "does not exist")
  expect_error(
    intensityToOD(x, assay_name = "raw", reference = c(1, 2, 3)),
    "length"
  )
  expect_error(
    intensityToOD(
      x, assay_name = "raw", reference = matrix(1, nrow = 1L)
    ),
    "finite numeric"
  )
  expect_error(
    intensityToOD(
      x, assay_name = "raw", reference = 1,
      baseline_interval = c(0, 1)
    ),
    "cannot both"
  )
  expect_error(
    intensityToOD(
      x, assay_name = "raw", baseline_interval = c(100, 101)
    ),
    "selects no"
  )
  expect_error(
    intensityToOD(x, assay_name = "raw", baseline = "mea"), "exactly"
  )
  expect_error(
    intensityToOD(x, assay_name = "raw", nonpositive = "rep"), "exactly"
  )
  expect_error(
    intensityToOD(x, assay_name = "raw", output_assay = "raw"),
    "already exists"
  )
  malformed <- make_cw_nirs()
  S4Vectors::metadata(malformed)$nirs <- list(
    assays = list(raw = list(kind = character()))
  )
  expect_error(intensityToOD(malformed), "Malformed assay contract")

  nonfinite <- make_cw_nirs()
  SummarizedExperiment::assay(nonfinite, "raw")[2L, 1L] <- Inf
  expect_error(intensityToOD(nonfinite), "row 2.*measurement 1")

  unsupported <- make_cw_nirs()
  table <- S4Vectors::metadata(unsupported)$snirf$measurement_list
  table$data_type[[2L]] <- 99999L
  S4Vectors::metadata(unsupported)$snirf$measurement_list <- table
  SummarizedExperiment::colData(unsupported)$data_type[[2L]] <- 99999L
  expect_error(intensityToOD(unsupported), "measurement\\(s\\) 2")
})
