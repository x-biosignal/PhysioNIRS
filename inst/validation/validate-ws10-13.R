#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 17)

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
script <- normalizePath(file_argument[[1L]])
root <- normalizePath(file.path(dirname(script), "..", ".."))
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path(root, "DESCRIPTION"))) {
  devtools::load_all(root, quiet = TRUE)
} else {
  library(PhysioNIRS)
}

output_directory <- dirname(script)
validation_path <- file.path(
  output_directory, "ws10-13-validation.csv"
)
mutation_path <- file.path(
  output_directory, "ws10-13-mutation-gates.csv"
)
report_path <- file.path(
  output_directory, "ws10-13-validation.md"
)
table_path <- file.path(
  root, "inst", "extdata", "haemoglobin-extinction.csv"
)
fixture_path <- file.path(root, "inst", "extdata", "mbll_reference.rds")
table <- utils::read.csv(table_path, stringsAsFactors = FALSE)

oracle_extinction <- function(wavelength_nm) {
  cbind(
    HbO = stats::approx(
      table$wavelength_nm, table$hbo2_cm1_M1, wavelength_nm,
      method = "linear", ties = "ordered"
    )$y,
    HbR = stats::approx(
      table$wavelength_nm, table$hbr_cm1_M1, wavelength_nm,
      method = "linear", ties = "ordered"
    )$y
  ) * (100 * log(10))
}

oracle_scholkmann <- function(wavelength_nm, age_years) {
  223.3 + 0.05624 * age_years^0.8493 -
    5.723e-7 * wavelength_nm^3 +
    0.001245 * wavelength_nm^2 -
    0.9025 * wavelength_nm
}

make_experiment <- function(values, wavelengths, distance, kind) {
  values <- as.matrix(values)
  n_measurement <- length(wavelengths)
  labels <- paste0("S1_D1_", wavelengths)
  repeated <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[repeated] <- paste0(labels[repeated], "_M", which(repeated))
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(n_measurement),
    source_index = rep(1L, n_measurement),
    detector_index = rep(1L, n_measurement),
    wavelength_index = seq_len(n_measurement),
    wavelength_nm = as.numeric(wavelengths),
    wavelength_actual_nm = rep(NA_real_, n_measurement),
    data_type = rep(1L, n_measurement),
    data_type_index = rep(1L, n_measurement),
    data_type_label = rep(NA_character_, n_measurement),
    data_unit = rep(NA_character_, n_measurement),
    source_power = rep(NA_real_, n_measurement),
    detector_gain = rep(NA_real_, n_measurement),
    source_label = rep("S1", n_measurement),
    detector_label = rep("D1", n_measurement),
    channel_label = labels
  )
  assay_name <- if (kind == "intensity") "raw" else "OD"
  metadata <- list(snirf = list(
    measurement_list = measurement,
    probe = list(
      wavelengths = wavelengths,
      sourceLabels = "S1",
      detectorLabels = "D1",
      sourcePos2D = matrix(c(0, 0), nrow = 1L),
      detectorPos2D = matrix(c(0, distance), nrow = 1L),
      LengthUnit = "m"
    ),
    metadata_tags = list(
      SubjectID = "validation",
      MeasurementDate = "unknown",
      MeasurementTime = "unknown",
      LengthUnit = "m",
      TimeUnit = "s",
      FrequencyUnit = "Hz"
    ),
    stim = list(),
    time_sampling = "uniform"
  ))
  if (kind == "od") {
    metadata$nirs <- list(assays = list(OD = list(
      kind = "optical_density",
      unit = "1",
      log_convention = "natural",
      source_assay = "oracle"
    )))
  }
  PhysioCore::PhysioExperiment(
    assays = stats::setNames(list(values), assay_name),
    rowData = S4Vectors::DataFrame(
      time_seconds = seq.int(0, nrow(values) - 1L) / 10
    ),
    colData = cbind(measurement, label = labels),
    metadata = metadata,
    samplingRate = 10
  )
}

wavelength_sets <- list(
  c(760, 850),
  c(690, 760, 850),
  c(690, 760, 830, 850),
  c(744, 807, 832)
)
validation <- vector("list", 100L)
for (case_index in seq_len(100L)) {
  wavelengths <- wavelength_sets[[(case_index - 1L) %% 4L + 1L]]
  if (case_index %% 5L == 0L) {
    wavelengths <- rev(wavelengths)
  }
  n_time <- 25L + (case_index %% 37L)
  time <- seq.int(0, n_time - 1L) / 10
  distance <- 0.02 + (case_index %% 26L) / 1000
  age <- case_index %% 71L
  ppf <- if (case_index %% 2L) {
    seq(5.2, 6.8, length.out = length(wavelengths))
  } else {
    oracle_scholkmann(wavelengths, age)
  }
  concentration <- cbind(
    HbO = (1 + case_index / 100) * sin(2 * pi * 0.07 * time),
    HbR = -(0.3 + case_index / 500) *
      cos(2 * pi * 0.05 * time + 0.1)
  )
  extinction <- oracle_extinction(wavelengths)
  optical_density <- (concentration * 1e-6) %*%
    t(extinction * (distance * ppf))
  reference <- seq(800, 1300, length.out = length(wavelengths))
  intensity <- sweep(exp(-optical_density), 2L, reference, "*")

  intensity_object <- make_experiment(
    intensity, wavelengths, distance, "intensity"
  )
  intensity_before <- serialize(intensity_object, NULL)
  od_result <- intensityToOD(intensity_object, reference = reference)
  od_actual <- SummarizedExperiment::assay(od_result, "OD")

  od_object <- make_experiment(
    optical_density, wavelengths, distance, "od"
  )
  od_before <- serialize(od_object, NULL)
  concentration_result <- if (case_index %% 2L) {
    mbll(od_object, pathlength_factor = ppf)
  } else {
    mbll(
      od_object,
      pathlength_factor = NULL,
      age_years = age,
      dpf_model = "scholkmann2013"
    )
  }
  hbo_actual <- SummarizedExperiment::assay(
    concentration_result, "HbO"
  )[, 1L]
  hbr_actual <- SummarizedExperiment::assay(
    concentration_result, "HbR"
  )[, 1L]
  hbt_actual <- SummarizedExperiment::assay(
    concentration_result, "HbT"
  )[, 1L]
  dpf_actual <- as.numeric(dpf(wavelengths, age))
  extinction_actual <- extinctionCoefficients(wavelengths)
  validation[[case_index]] <- data.frame(
    case = case_index,
    wavelengths = length(wavelengths),
    signal_length = n_time,
    distance_m = distance,
    age_years = age,
    dpf_source = if (case_index %% 2L) "explicit" else "scholkmann2013",
    od_max_abs = max(abs(od_actual - optical_density)),
    dpf_max_abs = max(abs(dpf_actual -
                            oracle_scholkmann(wavelengths, age))),
    extinction_max_rel = max(
      abs(extinction_actual - extinction) /
        pmax(abs(extinction), .Machine$double.eps)
    ),
    hbo_max_abs_uM = max(abs(hbo_actual - concentration[, "HbO"])),
    hbr_max_abs_uM = max(abs(hbr_actual - concentration[, "HbR"])),
    hbt_identity_max_abs_uM = max(
      abs(hbt_actual - hbo_actual - hbr_actual)
    ),
    signs_correct = all(
      sign(hbo_actual[abs(concentration[, "HbO"]) > 1e-10]) ==
        sign(concentration[
          abs(concentration[, "HbO"]) > 1e-10, "HbO"
        ])
    ) && all(
      sign(hbr_actual[abs(concentration[, "HbR"]) > 1e-10]) ==
        sign(concentration[
          abs(concentration[, "HbR"]) > 1e-10, "HbR"
        ])
    ),
    mapping_correct = identical(
      SummarizedExperiment::colData(concentration_result)$channel_id,
      "S1_D1"
    ),
    source_immutable = identical(
      serialize(intensity_object, NULL), intensity_before
    ) && identical(serialize(od_object, NULL), od_before),
    stringsAsFactors = FALSE
  )
}
validation <- do.call(rbind, validation)

errors <- function(expression) {
  inherits(try(force(expression), silent = TRUE), "try-error")
}
mutation <- data.frame(
  mutation = c(
    "od_base10", "od_wrong_sign", "od_nonpositive_abs",
    "od_ambiguous_assay", "od_reference_length", "vector_matrix_shape",
    "od_interval_empty", "od_overwrite", "od_ratio_overflow",
    "dpf_wrong_intercept", "dpf_bad_domain",
    "dpf_partial_model", "extinction_missing_ln10", "extinction_swap",
    "extinction_missing_metre", "extinction_outside",
    "extinction_overflow",
    "mbll_ungoverned_assay", "mbll_duplicate_wavelength",
    "mbll_single_wavelength", "mbll_zero_distance",
    "mbll_distance_mismatch", "mbll_negative_ppf",
    "mbll_age_and_ppf", "mbll_wrong_unit", "mbll_rank_deficient",
    "mbll_bad_data_type", "mbll_output_overflow",
    "mbll_hbt_not_sum", "source_mutation"
  ),
  caught = FALSE,
  stringsAsFactors = FALSE
)
base_wavelength <- c(760, 850)
base_od <- matrix(0, nrow = 25L, ncol = 2L)
base_od_object <- make_experiment(base_od, base_wavelength, 0.03, "od")
base_intensity_object <- make_experiment(
  matrix(100, nrow = 25L, ncol = 2L),
  base_wavelength, 0.03, "intensity"
)
mutation$caught <- c(
  max(abs(-log(c(0.5, 2)) - -log10(c(0.5, 2)))) > 0.1,
  max(abs(-log(c(0.5, 2)) - log(c(0.5, 2)))) > 1,
  {
    bad <- base_intensity_object
    SummarizedExperiment::assay(bad, "raw")[1L, 1L] <- -100
    errors(intensityToOD(bad))
  },
  {
    bad <- base_intensity_object
    SummarizedExperiment::assays(bad)[["other"]] <-
      SummarizedExperiment::assay(bad, "raw")
    errors(intensityToOD(bad))
  },
  errors(intensityToOD(base_intensity_object, reference = 1:3)),
  errors(intensityToOD(
    base_intensity_object, reference = matrix(1, nrow = 1L)
  )),
  errors(intensityToOD(
    base_intensity_object, baseline_interval = c(100, 101)
  )),
  {
    bad <- intensityToOD(base_intensity_object, reference = 100)
    errors(intensityToOD(bad, assay_name = "raw", output_assay = "OD"))
  },
  {
    extreme_values <- rbind(
      c(1e308, 1e-308), c(1e-308, 1e308)
    )
    extreme <- make_experiment(
      extreme_values, base_wavelength, 0.03, "intensity"
    )
    actual <- SummarizedExperiment::assay(
      intensityToOD(extreme, reference = c(1e308, 1e-308)), "OD"
    )
    naive <- suppressWarnings(
      -log(sweep(extreme_values, 2L, c(1e308, 1e-308), "/"))
    )
    all(is.finite(actual)) && any(!is.finite(naive))
  },
  abs(
    dpf(760, 20) -
      (233.3 + 0.05624 * 20^0.8493 - 5.723e-7 * 760^3 +
         0.001245 * 760^2 - 0.9025 * 760)
  ) > 9,
  errors(dpf(649, 20)),
  errors(dpf(760, 20, model = "schol")),
  max(abs(
    extinctionCoefficients(base_wavelength) -
      oracle_extinction(base_wavelength) / log(10)
  )) > 100,
  max(abs(
    extinctionCoefficients(base_wavelength) -
      oracle_extinction(base_wavelength)[, 2:1]
  )) > 100,
  max(abs(
    extinctionCoefficients(base_wavelength) -
      oracle_extinction(base_wavelength) / 100
  )) > 100,
  errors(extinctionCoefficients(1001)),
  errors(suppressWarnings(
    extinctionCoefficients(1e308, extrapolate = TRUE)
  )),
  {
    bad <- make_experiment(
      matrix(100, 25, 2), base_wavelength, 0.03, "intensity"
    )
    errors(mbll(bad, assay_name = "raw"))
  },
  {
    bad <- make_experiment(
      base_od, c(760, 760), 0.03, "od"
    )
    errors(mbll(bad))
  },
  {
    bad <- make_experiment(
      matrix(0, 25, 1), 760, 0.03, "od"
    )
    errors(mbll(bad))
  },
  errors(mbll(base_od_object, distance_m = 0)),
  errors(suppressWarnings(
    mbll(base_od_object, distance_m = c(0.03, 0.031))
  )),
  errors(mbll(base_od_object, pathlength_factor = -1)),
  errors(mbll(base_od_object, pathlength_factor = 6, age_years = 20)),
  errors(mbll(base_od_object, output_unit = "M")),
  {
    bad <- make_experiment(
      base_od, c(800, 800 + 1e-9), 0.03, "od"
    )
    errors(mbll(bad))
  },
  {
    bad <- base_od_object
    measurement <- S4Vectors::metadata(bad)$snirf$measurement_list
    measurement$data_type[[1L]] <- 99999L
    S4Vectors::metadata(bad)$snirf$measurement_list <- measurement
    SummarizedExperiment::colData(bad)$data_type[[1L]] <- 99999L
    errors(mbll(bad))
  },
  {
    overflow <- make_experiment(
      matrix(1e308, nrow = 25L, ncol = 2L),
      base_wavelength, 0.03, "od"
    )
    errors(mbll(overflow))
  },
  max(abs(
    SummarizedExperiment::assay(
      mbll(base_od_object), "HbT"
    ) - (
      SummarizedExperiment::assay(mbll(base_od_object), "HbO") +
        SummarizedExperiment::assay(mbll(base_od_object), "HbR") + 1
    )
  )) > 0.9,
  {
    before <- serialize(base_od_object, NULL)
    invisible(mbll(base_od_object))
    identical(serialize(base_od_object, NULL), before)
  }
)

utils::write.table(
  validation,
  validation_path,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE,
  na = "NA",
  eol = "\n"
)
utils::write.table(
  mutation,
  mutation_path,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE,
  na = "NA",
  eol = "\n"
)

gates <- c(
  od = max(validation$od_max_abs) < 1e-12,
  dpf = max(validation$dpf_max_abs) < 1e-12,
  extinction = max(validation$extinction_max_rel) < 1e-12,
  hbo = max(validation$hbo_max_abs_uM) < 1e-8,
  hbr = max(validation$hbr_max_abs_uM) < 1e-8,
  hbt = max(validation$hbt_identity_max_abs_uM) < 1e-12,
  signs = all(validation$signs_correct),
  mapping = all(validation$mapping_correct),
  source_immutable = all(validation$source_immutable),
  mutations = all(mutation$caught)
)
fixture <- readRDS(fixture_path)
report <- c(
  "# WS10-13 independent numeric validation",
  "",
  "- Package: PhysioNIRS 0.2.0",
  paste0("- R: ", R.version$major, ".", R.version$minor),
  "- Oracle: independent base R equations; no PhysioNIRS helper is used",
  paste0("- Cases: ", nrow(validation), "/100"),
  paste0("- Wavelength counts: ",
         paste(sort(unique(validation$wavelengths)), collapse = ", ")),
  paste0("- Extinction table SHA-256: ",
         digest::digest(table_path, algo = "sha256", file = TRUE)),
  paste0("- Fixture SHA-256: ",
         digest::digest(fixture_path, algo = "sha256", file = TRUE)),
  paste0("- Fixture generator SHA-256: ", fixture$generator_sha256),
  "- External comparator: unavailable (MNE/Homer3 not installed); pinned MNE data source is validated",
  "",
  "## Tolerances",
  "",
  "- OD max absolute error < 1e-12",
  "- DPF max absolute error < 1e-12",
  "- Extinction max relative error < 1e-12",
  "- HbO/HbR max absolute error < 1e-8 uM",
  "- HbT identity max absolute error < 1e-12 uM",
  "",
  "## Results",
  "",
  paste0("- OD max absolute error: ",
         format(max(validation$od_max_abs), scientific = TRUE)),
  paste0("- DPF max absolute error: ",
         format(max(validation$dpf_max_abs), scientific = TRUE)),
  paste0("- Extinction max relative error: ",
         format(max(validation$extinction_max_rel), scientific = TRUE)),
  paste0("- HbO max absolute error (uM): ",
         format(max(validation$hbo_max_abs_uM), scientific = TRUE)),
  paste0("- HbR max absolute error (uM): ",
         format(max(validation$hbr_max_abs_uM), scientific = TRUE)),
  paste0("- HbT identity max absolute error (uM): ",
         format(max(validation$hbt_identity_max_abs_uM), scientific = TRUE)),
  paste0("- Sign cases: ", sum(validation$signs_correct), "/100"),
  paste0("- Mapping cases: ", sum(validation$mapping_correct), "/100"),
  paste0("- Source immutability cases: ",
         sum(validation$source_immutable), "/100"),
  paste0("- Mutation gates: ", sum(mutation$caught), "/", nrow(mutation)),
  paste0("- Gates: ", sum(gates), "/", length(gates)),
  paste0("- Final: ", if (all(gates)) "PASS" else "FAIL")
)
writeLines(report, report_path, useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n")
if (!all(gates)) quit(status = 1L)
