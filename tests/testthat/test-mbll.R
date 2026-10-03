test_that("MBLL recovers independent forward models at 2, 3, and 4 wavelengths", {
  time <- seq.int(0, 40) / 10
  concentration <- cbind(
    HbO = 2.5 * sin(2 * pi * 0.1 * time),
    HbR = -0.8 * sin(2 * pi * 0.1 * time + 0.2)
  )
  for (wavelengths in list(
      c(760, 850),
      c(690, 760, 850),
      c(690, 760, 830, 850))) {
    ppf <- seq(5.5, 6.5, length.out = length(wavelengths))
    optical_density <- forward_mbll_od(
      concentration, wavelengths, distance = 0.032, ppf = ppf
    )
    x <- make_cw_nirs(
      wavelengths = wavelengths,
      distances = 0.032,
      assay_name = "OD",
      values = optical_density
    )
    result <- mbll(x, pathlength_factor = ppf)
    expect_equal(
      SummarizedExperiment::assay(result, "HbO")[, 1L],
      concentration[, "HbO"],
      tolerance = 1e-9
    )
    expect_equal(
      SummarizedExperiment::assay(result, "HbR")[, 1L],
      concentration[, "HbR"],
      tolerance = 1e-9
    )
    expect_equal(
      SummarizedExperiment::assay(result, "HbT"),
      SummarizedExperiment::assay(result, "HbO") +
        SummarizedExperiment::assay(result, "HbR"),
      tolerance = 1e-12
    )
  }
})

test_that("MBLL preserves pair order across wavelength permutations", {
  wavelengths <- c(690, 760, 850)
  time <- seq.int(0, 40) / 10
  concentration_one <- cbind(
    HbO = 3 * sin(time), HbR = -sin(time) / 2
  )
  concentration_two <- cbind(
    HbO = 1.5 * cos(time), HbR = -0.75 * cos(time)
  )
  ppf <- c(5.8, 6.1, 6.4)
  first <- forward_mbll_od(concentration_one, wavelengths, 0.02, ppf)
  second <- forward_mbll_od(concentration_two, wavelengths, 0.045, ppf)
  values <- cbind(first, second)
  x <- make_cw_nirs(
    wavelengths = wavelengths,
    distances = c(0.02, 0.045),
    assay_name = "OD",
    values = values
  )
  permutation <- c(3, 1, 2, 6, 4, 5)
  table <- S4Vectors::metadata(x)$snirf$measurement_list[permutation, ]
  x <- x[, permutation]
  table$measurement_index <- seq_len(nrow(table))
  S4Vectors::metadata(x)$snirf$measurement_list <- table
  SummarizedExperiment::colData(x)$measurement_index <- seq_len(nrow(table))
  result <- mbll(
    x,
    pathlength_factor = stats::setNames(ppf, wavelengths)
  )
  expect_identical(
    SummarizedExperiment::colData(result)$channel_id,
    c("S1_D1", "S2_D2")
  )
  expect_equal(
    SummarizedExperiment::assay(result, "HbO"),
    cbind(S1_D1 = concentration_one[, "HbO"],
          S2_D2 = concentration_two[, "HbO"]),
    tolerance = 1e-9
  )
  expect_equal(
    SummarizedExperiment::assay(result, "HbR"),
    cbind(S1_D1 = concentration_one[, "HbR"],
          S2_D2 = concentration_two[, "HbR"]),
    tolerance = 1e-9
  )
})

test_that("MBLL handles heterogeneous wavelength sets across pairs", {
  wavelengths <- c(690, 760, 830, 850)
  ppf <- c(5.4, 5.8, 6.2, 6.6)
  distances <- c(0.025, 0.035, 0.045)
  time <- seq.int(0, 40) / 10
  concentrations <- lapply(seq_along(distances), function(i) {
    cbind(
      HbO = i * sin(time),
      HbR = -i * cos(time) / 3
    )
  })
  full_values <- do.call(cbind, lapply(seq_along(distances), function(i) {
    forward_mbll_od(
      concentrations[[i]], wavelengths, distances[[i]], ppf
    )
  }))
  full <- make_cw_nirs(
    wavelengths = wavelengths,
    distances = distances,
    assay_name = "OD",
    values = full_values
  )
  keep <- c(1L, 2L, 5L, 6L, 7L, 9L, 10L, 11L, 12L)
  measurement <- S4Vectors::metadata(full)$snirf$measurement_list[keep, ]
  result_input <- full[, keep]
  measurement$measurement_index <- seq_len(nrow(measurement))
  S4Vectors::metadata(result_input)$snirf$measurement_list <- measurement
  SummarizedExperiment::colData(result_input)$measurement_index <-
    seq_len(nrow(measurement))
  result <- mbll(
    result_input,
    pathlength_factor = stats::setNames(ppf, wavelengths)
  )
  expect_equal(
    unname(SummarizedExperiment::assay(result, "HbO")),
    unname(do.call(
      cbind, lapply(concentrations, function(value) value[, "HbO"])
    )),
    tolerance = 1e-9
  )
  expect_identical(
    lengths(SummarizedExperiment::colData(result)$wavelengths_nm),
    c(2L, 3L, 4L)
  )
})

test_that("MBLL can derive wavelength-specific DPF", {
  wavelengths <- c(690, 760, 850)
  age <- 30
  ppf <- as.numeric(dpf(wavelengths, age))
  concentration <- cbind(
    HbO = seq(-2, 2, length.out = 41),
    HbR = seq(1, -1, length.out = 41)
  )
  optical_density <- forward_mbll_od(
    concentration, wavelengths, 0.03, ppf
  )
  x <- make_cw_nirs(
    wavelengths = wavelengths, distances = 0.03,
    assay_name = "OD", values = optical_density
  )
  result <- mbll(
    x,
    pathlength_factor = NULL,
    age_years = age,
    dpf_model = "scholkmann2013"
  )
  expect_equal(
    SummarizedExperiment::assay(result, "HbO")[, 1L],
    concentration[, 1L],
    tolerance = 1e-9
  )
  expect_identical(
    S4Vectors::metadata(result)$nirs$mbll$dpf_source,
    "scholkmann2013"
  )
})

test_that("pathlength factors align by wavelength when lengths are ambiguous", {
  wavelengths <- c(690, 760)
  concentration <- cbind(
    HbO = seq(-1, 1, length.out = 41),
    HbR = seq(0.5, -0.5, length.out = 41)
  )
  ppf <- c(5.5, 6.5)
  first <- forward_mbll_od(concentration, wavelengths, 0.03, ppf)
  values <- cbind(first, first)
  x <- make_cw_nirs(
    wavelengths = wavelengths,
    distances = c(0.03, 0.03),
    assay_name = "OD",
    values = values
  )
  result <- mbll(x, pathlength_factor = ppf)
  expect_equal(
    SummarizedExperiment::assay(result, "HbO"),
    cbind(S1_D1 = concentration[, "HbO"],
          S2_D2 = concentration[, "HbO"]),
    tolerance = 1e-9
  )

  one_pair <- make_cw_nirs(
    wavelengths = rev(wavelengths),
    distances = 0.03,
    assay_name = "OD",
    values = first[, 2:1, drop = FALSE]
  )
  named_result <- mbll(
    one_pair,
    pathlength_factor = c(`690` = 5.5, `760` = 6.5)
  )
  expect_equal(
    SummarizedExperiment::assay(named_result, "HbO")[, 1L],
    concentration[, "HbO"],
    tolerance = 1e-9
  )
  expect_error(
    mbll(
      one_pair,
      pathlength_factor = stats::setNames(c(5.5, 6.5), c("690", ""))
    ),
    "name each unique"
  )
})

test_that("MBLL returns a valid governed pair-collapsed object", {
  x <- make_cw_nirs(wavelengths = c(760, 850), distances = c(0.03, 0.04))
  x <- intensityToOD(x, reference = 1000)
  source_before <- serialize(x, NULL)
  expect_warning(
    result <- mbll(x, distance_m = c(0.03, 0.04)),
    "overrides"
  )
  expect_s4_class(result, "PhysioExperiment")
  expect_true(methods::validObject(result))
  expect_identical(dim(result), c(nrow(x), 2L))
  expect_identical(
    SummarizedExperiment::assayNames(result),
    c("HbO", "HbR", "HbT")
  )
  expect_identical(
    SummarizedExperiment::colData(result)$channel_id,
    c("S1_D1", "S2_D2")
  )
  expect_equal(
    SummarizedExperiment::colData(result)$source_detector_distance_m,
    c(0.03, 0.04)
  )
  expect_identical(
    S4Vectors::metadata(result)$events,
    S4Vectors::metadata(x)$events
  )
  expect_identical(
    SummarizedExperiment::rowData(result),
    SummarizedExperiment::rowData(x)
  )
  expect_identical(PhysioCore::samplingRate(result), 10)
  expect_identical(serialize(x, NULL), source_before)
  expect_identical(
    tail(PhysioCore::provenance(result)$activity, 1L),
    "mbll"
  )
  expect_identical(
    S4Vectors::metadata(result)$nirs$mbll$output_unit,
    "uM"
  )
  expect_identical(
    names(S4Vectors::metadata(result)$nirs$assays),
    c("HbO", "HbR", "HbT")
  )
  expect_identical(
    S4Vectors::metadata(result)$nirs$assays$HbR,
    list(
      kind = "haemoglobin_concentration",
      chromophore = "HbR",
      unit = "uM",
      source_assay = "OD",
      method = "mbll"
    )
  )
  expect_identical(
    S4Vectors::metadata(result)$nirs$mbll$source_assay_contract,
    S4Vectors::metadata(x)$nirs$assays$OD
  )
  expect_match(
    S4Vectors::metadata(result)$nirs$mbll$table_sha256,
    "^[0-9a-f]{64}$"
  )
})

test_that("MBLL rejects malformed measurement and numeric contracts", {
  x <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03, assay_name = "OD"
  )
  expect_error(
    mbll(x, pathlength_factor = 6, age_years = 20),
    "cannot both"
  )
  expect_error(
    mbll(x, pathlength_factor = NULL),
    "age_years"
  )
  expect_error(mbll(x, output_unit = "M"), "exactly")
  expect_error(mbll(x, dpf_model = "schol"), "exactly")
  expect_error(mbll(x, pathlength_factor = -1), "positive")
  expect_error(
    mbll(x, pathlength_factor = matrix(6, nrow = 1L)),
    "finite numeric"
  )
  expect_error(mbll(x, distance_m = c(0.03, 0.03, 0.03)), "length")
  expect_warning(
    expect_error(
      mbll(x, distance_m = c(0.03, 0.031)),
      "distance disagreement"
    ),
    "overrides"
  )

  ungoverned <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03, assay_name = "raw"
  )
  expect_error(
    mbll(ungoverned, assay_name = "raw"),
    "governed natural-log"
  )

  duplicate <- make_cw_nirs(
    wavelengths = c(760, 760), assay_name = "OD"
  )
  expect_error(mbll(duplicate), "duplicate wavelength")

  one <- make_cw_nirs(wavelengths = c(760, 850), assay_name = "OD")
  table <- S4Vectors::metadata(one)$snirf$measurement_list[1L, ]
  one <- one[, 1L]
  table$measurement_index <- 1L
  S4Vectors::metadata(one)$snirf$measurement_list <- table
  SummarizedExperiment::colData(one)$measurement_index <- 1L
  expect_error(mbll(one), "at least two")

  nonfinite <- make_cw_nirs(assay_name = "OD")
  SummarizedExperiment::assay(nonfinite, "OD")[1L, 1L] <- NA_real_
  expect_error(mbll(nonfinite), "row 1.*measurement 1")

  overflow <- make_cw_nirs(
    wavelengths = c(760, 850), distances = 0.03, assay_name = "OD",
    values = matrix(1e308, nrow = 41L, ncol = 2L)
  )
  expect_error(mbll(overflow), "non-finite concentration")

  unsupported <- make_cw_nirs(assay_name = "OD")
  table <- S4Vectors::metadata(unsupported)$snirf$measurement_list
  table$data_type[[1L]] <- 99999L
  S4Vectors::metadata(unsupported)$snirf$measurement_list <- table
  SummarizedExperiment::colData(unsupported)$data_type[[1L]] <- 99999L
  expect_error(mbll(unsupported), "unsupported data_type")
})

test_that("MBLL rank and conditioning gates are pair-specific", {
  x <- make_cw_nirs(
    wavelengths = c(800, 800 + 1e-9),
    assay_name = "OD"
  )
  expect_error(mbll(x), "ill-conditioned")
})
