#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PhysioNIRS)
  library(PhysioCore)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- normalizePath(file_arg[[1L]], mustWork = TRUE)
validation_dir <- dirname(script)
package_root <- normalizePath(
  file.path(validation_dir, "..", ".."), mustWork = TRUE
)

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

caught <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

make_experiment <- function(seed, n_seconds = 20) {
  set.seed(seed)
  fs <- c(10, 12, 20)[(seed %% 3L) + 1L]
  n_wavelength <- c(2L, 3L, 4L)[(seed %% 3L) + 1L]
  wavelength <- c(730, 760, 810, 850)[seq_len(n_wavelength)]
  time <- seq.int(0, n_seconds * fs - 1L) / fs
  shared <- sin(2 * pi * (0.9 + runif(1, -0.04, 0.04)) * time) +
    0.05 * sin(2 * pi * 0.2 * time)
  coupled <- vapply(seq_len(n_wavelength), function(j) {
    (0.8 + 0.1 * j) * shared +
      0.01 * sin(2 * pi * (2 + j / 20) * time + seed / 100)
  }, numeric(length(time)))
  decoupled <- vapply(seq_len(n_wavelength), function(j) {
    sin(2 * pi * (0.75 + 0.17 * j) * time + j + seed / 1000) +
      0.03 * sin(2 * pi * 0.2 * time)
  }, numeric(length(time)))
  values <- cbind(coupled, decoupled)
  source_index <- rep(1:2, each = n_wavelength)
  detector_index <- source_index
  wavelength_nm <- rep(wavelength, 2L)
  labels <- paste0(
    "S", source_index, "_D", detector_index, "_wl", wavelength_nm
  )
  colnames(values) <- labels
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_along(labels),
    source_index = as.integer(source_index),
    detector_index = as.integer(detector_index),
    wavelength_index = rep(seq_len(n_wavelength), 2L),
    wavelength_nm = as.numeric(wavelength_nm),
    wavelength_actual_nm = rep(NA_real_, length(labels)),
    data_type = rep(1L, length(labels)),
    data_type_index = rep(1L, length(labels)),
    data_type_label = rep(NA_character_, length(labels)),
    data_unit = rep(NA_character_, length(labels)),
    source_power = rep(NA_real_, length(labels)),
    detector_gain = rep(NA_real_, length(labels)),
    source_label = paste0("S", source_index),
    detector_label = paste0("D", detector_index),
    channel_label = labels
  )
  x <- PhysioCore::PhysioExperiment(
    assays = list(OD = values),
    rowData = S4Vectors::DataFrame(time_seconds = time),
    colData = cbind(measurement, label = labels),
    metadata = list(
      snirf = list(
        measurement_list = measurement,
        probe = list(
          wavelengths = wavelength,
          sourceLabels = c("S1", "S2"),
          detectorLabels = c("D1", "D2"),
          sourcePos2D = rbind(c(0, 0), c(0.05, 0)),
          detectorPos2D = rbind(c(0, 0.03), c(0.05, 0.03)),
          LengthUnit = "m"
        ),
        metadata_tags = list(
          SubjectID = "synthetic-validation",
          MeasurementDate = "unknown",
          MeasurementTime = "unknown",
          LengthUnit = "m",
          TimeUnit = "s",
          FrequencyUnit = "Hz"
        ),
        stim = list(),
        time_sampling = "uniform"
      ),
      nirs = list(assays = list(OD = list(
        kind = "optical_density",
        unit = "1",
        log_convention = "natural",
        source_assay = "synthetic"
      )))
    ),
    samplingRate = fs
  )
  list(
    x = x,
    fs = fs,
    time = time,
    n_wavelength = n_wavelength,
    expected = rep(c(TRUE, FALSE), each = n_wavelength)
  )
}

oracle_filter <- function(x, fs, lower = 0.7, upper = 1.5, order = 3L) {
  filter <- signal::butter(order, c(lower, upper) / (fs / 2), "pass")
  b <- as.numeric(filter$b)
  a <- as.numeric(filter$a)
  gain <- sum(b) / sum(a)
  vapply(seq_len(ncol(x)), function(j) {
    forward <- signal::filter(
      filter, x[, j],
      init.x = rep(x[[1L, j]], length(b) - 1L),
      init.y = rep(x[[1L, j]] * gain, length(a) - 1L)
    )
    backward_input <- rev(as.numeric(forward))
    backward <- signal::filter(
      filter, backward_input,
      init.x = rep(backward_input[[1L]], length(b) - 1L),
      init.y = rep(backward_input[[1L]] * gain, length(a) - 1L)
    )
    rev(as.numeric(backward))
  }, numeric(nrow(x)))
}

oracle_sci <- function(case) {
  values <- SummarizedExperiment::assay(case$x, "OD")
  filtered <- oracle_filter(values, case$fs)
  measurement <- measurementList(case$x)
  key <- paste(measurement$source_index, measurement$detector_index)
  score <- numeric(ncol(values))
  for (index in split(seq_along(key), factor(key, levels = unique(key)))) {
    pairs <- utils::combn(index, 2L, simplify = FALSE)
    pair_score <- vapply(pairs, function(pair) {
      stats::cor(filtered[, pair[[1L]]], filtered[, pair[[2L]]])
    }, numeric(1))
    score[index] <- min(pair_score)
  }
  names(score) <- colnames(values)
  score
}

reorder_experiment <- function(x, index) {
  out <- x[, index, drop = FALSE]
  metadata <- S4Vectors::metadata(out)
  measurement <- metadata$snirf$measurement_list[index, , drop = FALSE]
  measurement$measurement_index <- seq_len(nrow(measurement))
  metadata$snirf$measurement_list <- measurement
  S4Vectors::metadata(out) <- metadata
  columns <- SummarizedExperiment::colData(out)
  columns$measurement_index <- seq_len(nrow(columns))
  SummarizedExperiment::colData(out) <- columns
  out
}

cases <- vector("list", 120L)
for (i in seq_len(120L)) {
  case <- make_experiment(16000L + i)
  input_before <- serialize(case$x, NULL, version = 3L)
  expected <- oracle_sci(case)
  actual <- scalpCouplingIndex(case$x, window_seconds = NULL)
  actual_score <- actual$score[, 1L]
  set.seed(17000L + i)
  permutation <- sample(seq_along(actual_score))
  permuted <- reorder_experiment(case$x, permutation)
  permuted_score <- scalpCouplingIndex(
    permuted, window_seconds = NULL
  )$score[, 1L]
  permuted_score <- permuted_score[names(actual_score)]

  offset <- case$x
  offset_value <- SummarizedExperiment::assay(offset, "OD") + 1000
  SummarizedExperiment::assay(
    offset, "OD", withDimnames = FALSE
  ) <- offset_value
  scaled <- case$x
  scaled_value <- SummarizedExperiment::assay(scaled, "OD") * 7.5
  SummarizedExperiment::assay(
    scaled, "OD", withDimnames = FALSE
  ) <- scaled_value
  offset_score <- scalpCouplingIndex(
    offset, window_seconds = NULL
  )$score[, 1L]
  scale_score <- scalpCouplingIndex(
    scaled, window_seconds = NULL
  )$score[, 1L]

  cases[[i]] <- data.frame(
    case = i,
    seed = 16000L + i,
    sampling_rate_hz = case$fs,
    wavelength_count = case$n_wavelength,
    channel_count = length(actual_score),
    sci_max_abs_error = max(abs(actual_score - expected)),
    classification_correct = sum(
      (actual_score >= 0.8) == case$expected
    ),
    classification_total = length(case$expected),
    permutation_max_abs = max(abs(
      actual_score - permuted_score
    )),
    offset_max_abs = max(abs(actual_score - offset_score)),
    scale_max_abs = max(abs(actual_score - scale_score)),
    pair_copy_exact = all(vapply(
      split(
        seq_along(actual_score),
        rep(1:2, each = case$n_wavelength)
      ),
      function(index) length(unique(actual_score[index])) == 1L,
      logical(1)
    )),
    input_immutable = identical(
      serialize(case$x, NULL, version = 3L), input_before
    ),
    finite = all(is.finite(actual_score)),
    stringsAsFactors = FALSE
  )
}
case_table <- do.call(rbind, cases)

peak_cases <- lapply(seq_len(20L), function(i) {
  case <- make_experiment(18000L + i)
  result <- signalQualityIndex(
    case$x, method = "peak_power", window_seconds = 20
  )
  data.frame(
    case = i,
    agreement = sum((result$score[, 1L] >= 0.1) == case$expected),
    total = length(case$expected),
    pair_copy_exact = all(vapply(
      split(
        seq_len(nrow(result$score)),
        rep(1:2, each = case$n_wavelength)
      ),
      function(index) length(unique(result$score[index, 1L])) == 1L,
      logical(1)
    )),
    finite = all(is.finite(result$score)),
    stringsAsFactors = FALSE
  )
})
peak_table <- do.call(rbind, peak_cases)

fixture_path <- file.path(
  package_root, "inst", "extdata", "quality_reference.rds"
)
fixture <- readRDS(fixture_path)
snr <- fixture$analytic_snr
snr_case <- make_experiment(16001L)
snr_values <- cbind(snr$values, snr$values)
snr_case$x <- snr_case$x[, 1:2, drop = FALSE]
snr_meta <- S4Vectors::metadata(snr_case$x)
snr_meta$snirf$measurement_list <-
  snr_meta$snirf$measurement_list[1:2, , drop = FALSE]
snr_meta$snirf$measurement_list$measurement_index <- 1:2
S4Vectors::metadata(snr_case$x) <- snr_meta
SummarizedExperiment::assay(
  snr_case$x, "OD", withDimnames = FALSE
) <- snr_values
SummarizedExperiment::rowData(snr_case$x) <-
  S4Vectors::DataFrame(time_seconds = snr$time_seconds)
methods::slot(snr_case$x, "samplingRate") <- snr$sampling_rate_hz
snr_result <- signalQualityIndex(
  snr_case$x,
  method = "snr",
  window_seconds = 20,
  cardiac_range_hz = snr$signal_band_hz,
  noise_range_hz = snr$noise_band_hz
)
snr_error <- max(abs(snr_result$score - snr$expected_snr_db))

prune_case <- make_experiment(19001L)
prune_quality <- scalpCouplingIndex(
  prune_case$x, window_seconds = NULL
)
marked <- pruneChannels(prune_case$x, prune_quality)
dropped <- pruneChannels(
  prune_case$x, prune_quality, action = "drop"
)
mark_exact <- identical(
  as.logical(SummarizedExperiment::colData(marked)$nirs_quality_pass),
  prune_case$expected
)
drop_exact <- ncol(dropped) == prune_case$n_wavelength &&
  all(measurementList(dropped)$source_index == 1L) &&
  identical(
    SummarizedExperiment::assay(dropped, "OD"),
    SummarizedExperiment::assay(prune_case$x, "OD")[
      , seq_len(prune_case$n_wavelength), drop = FALSE
    ]
  )

make_live_source <- function(live) {
  info <- PhysioStream::streamInfo(
    "ws10-16-validation",
    type = "NIRS",
    channel_names = live$channel_id,
    nominal_srate = live$sampling_rate_hz,
    dtype = "float64",
    source_id = "ws10-16-validation",
    clock_domain = "ws10-16-clock",
    channel_units = rep("uM", length(live$channel_id)),
    metadata = list(nirs = list(
      assay_name = "HbO",
      assay_kind = "haemoglobin_concentration",
      identity_kind = "pair",
      channel_id = live$channel_id,
      source_index = live$source_index,
      detector_index = live$detector_index
    ))
  )
  PhysioStream::streamOpen(PhysioStream::loopbackSource(info, 4096L))
}

oracle_weighted <- function(time, values, end, duration) {
  start <- end - duration
  interior <- time > start & time < end
  grid <- c(start, time[interior], end)
  vapply(seq_len(ncol(values)), function(j) {
    ordinate <- stats::approx(
      time, values[, j], xout = grid, ties = "ordered"
    )$y
    sum(diff(grid) * (head(ordinate, -1L) + tail(ordinate, -1L)) / 2) /
      duration
  }, numeric(1))
}

live <- fixture$live
live_source <- make_live_source(live)
live_controller <- nirsNeurofeedback(
  live_source, live$regions, live$contrast,
  baseline_seconds = live$baseline_seconds,
  update_seconds = live$update_seconds,
  smoothing_seconds = live$smoothing_seconds
)
nirsNeurofeedbackStart(live_controller)
PhysioStream::loopbackFeed(
  live_source, live$values, live$time_seconds
)
repeat {
  step <- nirsNeurofeedbackStep(live_controller, 4L)
  if (!step$updated) break
}
live_state <- nirsNeurofeedbackState(live_controller)
region_values <- cbind(
  left = rowMeans(live$values[, live$regions$left, drop = FALSE]),
  right = live$values[, live$regions$right]
)
expected_live <- vapply(
  live_state$update_timestamps,
  function(timestamp) {
    smoothed <- oracle_weighted(
      live$time_seconds, region_values, timestamp,
      live$smoothing_seconds
    )
    sum(live$contrast * (smoothed - live$expected_baseline))
  },
  numeric(1)
)
actual_live <- vapply(
  live_state$target_values,
  function(value) unname(value[[1L]]),
  numeric(1)
)
live_error <- max(abs(expected_live - actual_live))
live_sequence <- vapply(
  live_state$receipts,
  function(value) value$receipt$sequence,
  numeric(1)
)
live_atomic <- all(vapply(
  live_state$receipts,
  function(value) identical(value$receipt$names, "feedback"),
  logical(1)
))
nirsNeurofeedbackStop(live_controller)

mutation_case <- make_experiment(20001L)
mutation_quality <- scalpCouplingIndex(
  mutation_case$x, window_seconds = NULL
)
tampered_score <- mutation_quality
tampered_score$score[[1L]] <- tampered_score$score[[1L]] - 0.1
tampered_pass <- mutation_quality
tampered_pass$pass[[1L]] <- !tampered_pass$pass[[1L]]
nonuniform <- mutation_case$x
SummarizedExperiment::rowData(nonuniform)$time_seconds[[5L]] <-
  SummarizedExperiment::rowData(nonuniform)$time_seconds[[5L]] + 0.01
nonfinite <- mutation_case$x
nonfinite_value <- SummarizedExperiment::assay(nonfinite, "OD")
nonfinite_value[[1L, 1L]] <- Inf
SummarizedExperiment::assay(
  nonfinite, "OD", withDimnames = FALSE
) <- nonfinite_value
stale <- mutation_case$x
stale_value <- SummarizedExperiment::assay(stale, "OD")
stale_value[[1L, 1L]] <- stale_value[[1L, 1L]] + 1e-6
SummarizedExperiment::assay(stale, "OD", withDimnames = FALSE) <- stale_value
all_bad <- make_experiment(20002L)
all_bad_value <- SummarizedExperiment::assay(all_bad$x, "OD")
for (j in seq_len(ncol(all_bad_value))) {
  all_bad_value[, j] <- sin(
    2 * pi * (0.75 + 0.16 * j) * all_bad$time + j
  )
}
SummarizedExperiment::assay(
  all_bad$x, "OD", withDimnames = FALSE
) <- all_bad_value
all_bad_quality <- scalpCouplingIndex(
  all_bad$x, window_seconds = NULL
)
zero <- mutation_case$x
zero_value <- SummarizedExperiment::assay(zero, "OD")
zero_value[, 1L] <- 1
SummarizedExperiment::assay(zero, "OD", withDimnames = FALSE) <- zero_value
zero_quality <- scalpCouplingIndex(zero, window_seconds = NULL)
equality_quality <- scalpCouplingIndex(
  mutation_case$x, window_seconds = NULL,
  threshold = mutation_quality$score[[1L]]
)

bad_live_source <- make_live_source(live)
clock_controller <- nirsNeurofeedback(
  bad_live_source, live$regions, live$contrast,
  baseline_seconds = 1, update_seconds = 0.2, smoothing_seconds = 0.4
)
nirsNeurofeedbackStart(clock_controller)
bad_time <- c(0, 0.1, 0.2, 0.5)
PhysioStream::loopbackFeed(
  bad_live_source, live$values[1:4, , drop = FALSE], bad_time
)
clock_caught <- caught(nirsNeurofeedbackStep(clock_controller, 4L))

mutations <- data.frame(
  mutation = c(
    "partial_method", "invalid_low_band", "nyquist_endpoint",
    "nonexact_window", "partial_tail_excluded", "nonuniform_time",
    "nonfinite_assay", "stale_source", "tampered_score", "tampered_pass",
    "zero_variance_score", "threshold_equality", "noise_vector",
    "noise_overlap", "zero_noise_power", "mark_default", "drop_pair_exact",
    "drop_all_rejected", "duplicate_metric", "require_all_na",
    "partial_action", "region_missing", "region_overlap",
    "positional_contrast", "baseline_zero", "update_too_fast",
    "clock_discontinuity", "no_prebaseline_update",
    "sequence_contiguous", "atomic_target_names",
    "sample_count_schedule_rejected", "duplicate_delivery_absent"
  ),
  detected = c(
    caught(signalQualityIndex(
      mutation_case$x, method = "peak", window_seconds = 10
    )),
    caught(scalpCouplingIndex(mutation_case$x, l_freq = 0)),
    caught(scalpCouplingIndex(
      mutation_case$x, h_freq = mutation_case$fs / 2
    )),
    caught(scalpCouplingIndex(
      mutation_case$x, window_seconds = 1.05
    )),
    {
      tail_case <- make_experiment(20003L)
      tail_result <- scalpCouplingIndex(
        tail_case$x, window_seconds = 7, step_seconds = 7
      )
      all(tail_result$window_end_sample <= nrow(tail_case$x)) &&
        tail(tail_result$window_end_sample, 1L) <= nrow(tail_case$x)
    },
    caught(scalpCouplingIndex(nonuniform)),
    caught(scalpCouplingIndex(nonfinite)),
    caught(pruneChannels(stale, mutation_quality)),
    caught(pruneChannels(mutation_case$x, tampered_score)),
    caught(pruneChannels(mutation_case$x, tampered_pass)),
    all(zero_quality$score[1:mutation_case$n_wavelength, ] == 0),
    all(equality_quality$pass[1:mutation_case$n_wavelength, ]),
    caught(signalQualityIndex(
      mutation_case$x, method = "snr", window_seconds = 10,
      noise_range_hz = c(2, 3)
    )),
    caught(signalQualityIndex(
      mutation_case$x, method = "snr", window_seconds = 10,
      noise_range_hz = matrix(c(1, 2), nrow = 1L)
    )),
    caught(signalQualityIndex(
      zero, method = "snr", window_seconds = 10,
      noise_range_hz = matrix(c(2, 4), nrow = 1L)
    )),
    all(c(
      "nirs_quality_pass", "nirs_quality_fingerprint"
    ) %in% names(SummarizedExperiment::colData(marked))),
    drop_exact,
    caught(pruneChannels(
      all_bad$x, all_bad_quality, action = "drop"
    )),
    caught(pruneChannels(
      mutation_case$x,
      list(a = mutation_quality, b = mutation_quality)
    )),
    caught(pruneChannels(
      mutation_case$x, mutation_quality, require_all = NA
    )),
    caught(pruneChannels(
      mutation_case$x, mutation_quality, action = "dro"
    )),
    caught(nirsNeurofeedback(
      make_live_source(live),
      list(left = "missing", right = "S3_D3"),
      live$contrast, baseline_seconds = 1
    )),
    caught(nirsNeurofeedback(
      make_live_source(live),
      list(left = c("S1_D1", "S2_D2"), right = "S2_D2"),
      live$contrast, baseline_seconds = 1
    )),
    caught(nirsNeurofeedback(
      make_live_source(live), live$regions, c(1, -1),
      baseline_seconds = 1
    )),
    caught(nirsNeurofeedback(
      make_live_source(live), live$regions, live$contrast,
      baseline_seconds = 0
    )),
    caught(nirsNeurofeedback(
      make_live_source(live), live$regions, live$contrast,
      baseline_seconds = 1, update_seconds = 0.001
    )),
    clock_caught,
    min(live_state$update_timestamps) >= live$baseline_seconds,
    identical(live_sequence, seq(0, length.out = length(live_sequence))),
    live_atomic,
    all(abs(diff(live_state$update_timestamps) -
      live$update_seconds) < 1e-10),
    length(unique(live_sequence)) == length(live_sequence)
  ),
  stringsAsFactors = FALSE
)

classification_correct <- sum(case_table$classification_correct)
classification_total <- sum(case_table$classification_total)
coupled_correct <- sum(vapply(seq_len(120L), function(i) {
  case <- make_experiment(16000L + i)
  result <- scalpCouplingIndex(case$x, window_seconds = NULL)
  sum(result$score[seq_len(case$n_wavelength), 1L] >= 0.8)
}, integer(1)))
coupled_total <- sum(case_table$wavelength_count)
decoupled_correct <- classification_correct - coupled_correct
decoupled_total <- classification_total - coupled_total

manifest_path <- file.path(
  package_root, "inst", "extdata", "quality_reference.sha256"
)
manifest <- readLines(manifest_path, warn = FALSE)
fixture_valid <- length(manifest) >= 2L &&
  identical(
    strsplit(manifest[[1L]], " ", fixed = TRUE)[[1L]][[1L]],
    sha256_file(fixture_path)
  )

gates <- data.frame(
  gate = c(
    "quality_cases_120", "classification_agreement_90pct",
    "classification_sensitivity_90pct",
    "classification_specificity_90pct", "sci_error_1e_10",
    "peak_classification_95pct", "pair_copy_exact",
    "snr_error_0_1_db", "permutation_equivariance",
    "offset_invariance", "scale_invariance", "input_immutable",
    "mark_exact", "drop_exact", "live_target_error_1e_10",
    "live_delivery_exact", "all_finite", "mutations_25",
    "fixture_manifest_valid"
  ),
  pass = c(
    nrow(case_table) == 120L,
    classification_correct / classification_total > 0.90,
    coupled_correct / coupled_total >= 0.90,
    decoupled_correct / decoupled_total >= 0.90,
    max(case_table$sci_max_abs_error) <= 1e-10,
    sum(peak_table$agreement) / sum(peak_table$total) >= 0.95,
    all(case_table$pair_copy_exact) && all(peak_table$pair_copy_exact),
    snr_error <= 0.1,
    max(case_table$permutation_max_abs) <= 1e-10,
    max(case_table$offset_max_abs) <= 1e-9,
    max(case_table$scale_max_abs) <= 1e-9,
    all(case_table$input_immutable),
    mark_exact,
    drop_exact,
    live_error <= 1e-10,
    length(live_sequence) > 0L &&
      identical(live_sequence, seq(0, length.out = length(live_sequence))),
    all(case_table$finite) && all(peak_table$finite) &&
      all(is.finite(actual_live)),
    nrow(mutations) >= 25L && all(mutations$detected),
    fixture_valid
  ),
  stringsAsFactors = FALSE
)

options(digits = 17)
utils::write.csv(
  case_table,
  file.path(validation_dir, "ws10-16-validation.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  mutations,
  file.path(validation_dir, "ws10-16-mutations.csv"),
  row.names = FALSE,
  na = ""
)
report <- c(
  "# WS10-16 independent validation",
  "",
  paste("- Quality cases:", nrow(case_table)),
  paste(
    "- SCI classification:",
    classification_correct, "/", classification_total,
    "agreement =", format(
      classification_correct / classification_total, digits = 17
    )
  ),
  paste(
    "- Sensitivity/specificity:",
    format(coupled_correct / coupled_total, digits = 17), "/",
    format(decoupled_correct / decoupled_total, digits = 17)
  ),
  paste("- SCI maximum absolute error:", format(
    max(case_table$sci_max_abs_error), digits = 17
  )),
  paste("- Peak-power classification:", sum(peak_table$agreement), "/",
        sum(peak_table$total)),
  paste("- Analytic SNR error dB:", format(snr_error, digits = 17)),
  paste("- Permutation maximum error:", format(
    max(case_table$permutation_max_abs), digits = 17
  )),
  paste("- Offset/scale maximum error:", format(
    max(case_table$offset_max_abs), digits = 17
  ), "/", format(max(case_table$scale_max_abs), digits = 17)),
  paste("- Live updates:", length(live_sequence),
        "; target maximum error:", format(live_error, digits = 17)),
  paste("- Mutations:", sum(mutations$detected), "/", nrow(mutations)),
  paste("- Gates:", sum(gates$pass), "/", nrow(gates)),
  paste("- Fixture SHA-256:", sha256_file(fixture_path)),
  paste("- Validator SHA-256:", sha256_file(script)),
  "",
  "Runtime equations are independent base-R/signal implementations; no",
  "PhysioNIRS private helper is called.",
  "",
  "## Gate table",
  "",
  "| gate | pass |",
  "|---|---|",
  paste0("| ", gates$gate, " | ", gates$pass, " |")
)
writeLines(
  report,
  file.path(validation_dir, "ws10-16-validation.md"),
  useBytes = TRUE
)

if (!all(gates$pass)) {
  print(gates)
  stop("WS10-16 validation gates failed", call. = FALSE)
}
cat(
  "WS10-16 validation:",
  nrow(case_table), "quality cases;",
  nrow(mutations), "mutations;",
  sum(gates$pass), "/", nrow(gates), "gates PASS\n"
)
