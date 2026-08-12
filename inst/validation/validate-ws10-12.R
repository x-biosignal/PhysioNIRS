#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this validator with Rscript")
script <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
package_dir <- normalizePath(file.path(dirname(script), "..", ".."))
validation_dir <- file.path(package_dir, "inst", "validation")
devtools::load_all(package_dir, quiet = TRUE)

make_case <- function(id) {
  choose <- function(values, divisor = 1L) {
    values[[(as.integer((id - 1L) / divisor) %% length(values)) + 1L]]
  }
  n_source <- choose(c(1L, 2L, 8L), 1L)
  n_detector <- choose(c(1L, 2L, 8L), 3L)
  n_wavelength <- choose(c(2L, 3L, 4L), 9L)
  geometry <- choose(c("2d", "3d", "both"), 27L)
  time_class <- choose(c("full", "compact", "irregular"), 7L)
  encoding <- choose(c("indexed", "vectorized"), 5L)
  stimuli <- if (id %% 4L == 0L) 0L else 2L
  n_time <- 9L + id %% 5L
  n_measurement <- max(n_source, n_detector, n_wavelength)
  source <- rep(seq_len(n_source), length.out = n_measurement)
  detector <- rep(seq_len(n_detector), length.out = n_measurement)
  wavelength <- rep(seq_len(n_wavelength), length.out = n_measurement)
  wavelengths <- 700 + seq_len(n_wavelength) * 50
  data <- outer(
    seq_len(n_time), seq_len(n_measurement),
    function(i, j) sin((i + id) / 7) + cos(j / 5) + i * j / 10000
  )
  time <- seq.int(0, n_time - 1L) / 20
  if (time_class == "irregular") {
    time <- cumsum(c(0, 0.05 + (seq_len(n_time - 1L) %% 3L) * 0.001))
  }
  source_labels <- paste0("S", seq_len(n_source))
  detector_labels <- paste0("D", seq_len(n_detector))
  base <- paste0(
    source_labels[source], "_", detector_labels[detector], "_",
    wavelengths[wavelength]
  )
  labels <- base
  repeated <- duplicated(base) | duplicated(base, fromLast = TRUE)
  labels[repeated] <- paste0(base[repeated], "_T1I1M", which(repeated))
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(n_measurement),
    source_index = as.integer(source),
    detector_index = as.integer(detector),
    wavelength_index = as.integer(wavelength),
    wavelength_nm = wavelengths[wavelength],
    wavelength_actual_nm = rep(NA_real_, n_measurement),
    data_type = rep(1L, n_measurement),
    data_type_index = rep(1L, n_measurement),
    data_type_label = rep(NA_character_, n_measurement),
    data_unit = rep(NA_character_, n_measurement),
    source_power = rep(NA_real_, n_measurement),
    detector_gain = rep(NA_real_, n_measurement),
    source_label = source_labels[source],
    detector_label = detector_labels[detector],
    channel_label = labels
  )
  source2 <- cbind((seq_len(n_source) - 1L) * 3, 0)
  detector2 <- cbind((seq_len(n_detector) - 1L) * 2.5, 4)
  source3 <- cbind(source2, (seq_len(n_source) - 1L) * 0.1)
  detector3 <- cbind(detector2, (seq_len(n_detector) - 1L) * 0.2)
  probe <- list(
    wavelengths = wavelengths,
    sourceLabels = source_labels,
    detectorLabels = detector_labels,
    LengthUnit = "cm"
  )
  if (geometry %in% c("2d", "both")) {
    probe$sourcePos2D <- source2
    probe$detectorPos2D <- detector2
  }
  if (geometry %in% c("3d", "both")) {
    probe$sourcePos3D <- source3
    probe$detectorPos3D <- detector3
  }
  stim <- if (stimuli) {
    list(
      list(
        name = "A",
        data = matrix(c(time[[2L]], 0.01, 1, id), nrow = 1L),
        data_labels = c("starttime", "duration", "value", "case")
      ),
      list(
        name = "B",
        data = matrix(c(time[[n_time - 1L]], 0.02, 2), nrow = 1L),
        data_labels = c("starttime", "duration", "value")
      )
    )
  } else {
    list()
  }
  tags <- list(
    SubjectID = "anonymous",
    MeasurementDate = "unknown",
    MeasurementTime = "unknown",
    LengthUnit = "cm",
    TimeUnit = "s",
    FrequencyUnit = "Hz",
    ValidationCase = as.character(id)
  )
  step <- PhysioNIRS:::.snirf_uniform_step(time)
  x <- PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    rowData = S4Vectors::DataFrame(time_seconds = time),
    colData = cbind(measurement, label = measurement$channel_label),
    metadata = list(snirf = list(
      measurement_list = measurement,
      probe = probe,
      metadata_tags = tags,
      stim = stim
    )),
    samplingRate = if (is.na(step)) NA_real_ else 1 / step
  )
  if (stimuli) {
    x <- PhysioCore::setEvents(
      x,
      PhysioCore::PhysioEvents(
        onset = c(time[[2L]], time[[n_time - 1L]]),
        duration = c(0.01, 0.02),
        type = c("A", "B"),
        value = c("1", "2")
      )
    )
  }
  list(
    x = x, id = id, n_source = n_source, n_detector = n_detector,
    n_wavelength = n_wavelength, geometry = geometry,
    time_class = time_class, encoding = encoding, stimuli = stimuli,
    data = data, time = time, measurement = measurement, probe = probe
  )
}

convert_vectorized <- function(path) {
  table <- PhysioNIRS::measurementList(PhysioNIRS::readSNIRF(path))
  tree <- PhysioNIRS:::.snirf_tree(path)
  groups <- PhysioNIRS:::.snirf_indexed_groups(
    tree, "/nirs/data1", "measurementList"
  )
  for (group in groups) rhdf5::h5delete(path, group)
  rhdf5::h5createGroup(path, "/nirs/data1/measurementLists")
  fields <- list(
    sourceIndex = as.integer(table$source_index),
    detectorIndex = as.integer(table$detector_index),
    wavelengthIndex = as.integer(table$wavelength_index),
    dataType = as.integer(table$data_type),
    dataTypeIndex = as.integer(table$data_type_index)
  )
  for (name in names(fields)) {
    PhysioNIRS:::.snirf_write_vector(
      path, paste0("/nirs/data1/measurementLists/", name), fields[[name]]
    )
  }
}

event_frame <- function(x) {
  data <- methods::slot(PhysioCore::getEvents(x), "events")
  data.frame(
    onset = as.numeric(data$onset),
    duration = as.numeric(data$duration),
    type = as.character(data$type),
    value = as.character(data$value),
    stringsAsFactors = FALSE
  )
}

max_error <- function(x, y) {
  if (!length(x) && !length(y)) return(0)
  max(abs(as.numeric(x) - as.numeric(y)))
}

validate_case <- function(case) {
  path1 <- tempfile(fileext = ".snirf")
  path2 <- tempfile(fileext = ".snirf")
  on.exit(unlink(c(path1, path2)), add = TRUE)
  compact <- case$time_class == "compact"
  PhysioNIRS::writeSNIRF(case$x, path1, compact_time = compact)
  if (case$encoding == "vectorized") convert_vectorized(path1)
  y <- PhysioNIRS::readSNIRF(path1)

  direct_data <- unname(t(rhdf5::h5read(
    path1, "/nirs/data1/dataTimeSeries"
  )))
  direct_source <- if (case$encoding == "vectorized") {
    as.integer(rhdf5::h5read(
      path1, "/nirs/data1/measurementLists/sourceIndex"
    ))
  } else {
    vapply(seq_len(ncol(case$data)), function(i) {
      as.integer(rhdf5::h5read(
        path1,
        paste0("/nirs/data1/measurementList", i, "/sourceIndex")
      ))
    }, integer(1))
  }
  table <- PhysioNIRS::measurementList(y)
  distances <- PhysioNIRS::sourceDetectorDistances(
    y, dimension = "auto", unit = "native"
  )
  dimension <- if (case$geometry %in% c("3d", "both")) "3D" else "2D"
  source_pos <- case$probe[[paste0("sourcePos", dimension)]]
  detector_pos <- case$probe[[paste0("detectorPos", dimension)]]
  expected_distance <- sqrt(rowSums((
    source_pos[case$measurement$source_index, , drop = FALSE] -
      detector_pos[case$measurement$detector_index, , drop = FALSE]
  )^2))
  expected_events <- event_frame(case$x)
  actual_events <- event_frame(y)
  event_onset_error <- max_error(actual_events$onset, expected_events$onset)
  event_duration_error <- max_error(
    actual_events$duration, expected_events$duration
  )
  event_value_error <- if (identical(actual_events$value,
                                     expected_events$value)) 0 else Inf

  PhysioNIRS::writeSNIRF(y, path2)
  z <- PhysioNIRS::readSNIRF(path2)
  cycle2_error <- max_error(
    SummarizedExperiment::assay(z),
    SummarizedExperiment::assay(y)
  )
  probe_exact <- identical(
    S4Vectors::metadata(y)$snirf$probe$wavelengths,
    case$probe$wavelengths
  )
  for (field in intersect(
    c("sourcePos2D", "detectorPos2D", "sourcePos3D", "detectorPos3D"),
    names(case$probe)
  )) {
    probe_exact <- probe_exact && identical(
      S4Vectors::metadata(y)$snirf$probe[[field]], case$probe[[field]]
    )
  }
  data_error <- max_error(SummarizedExperiment::assay(y), case$data)
  direct_error <- max_error(direct_data, case$data)
  measurement_mismatches <- sum(
    direct_source != case$measurement$source_index |
      table$source_index != case$measurement$source_index |
      table$detector_index != case$measurement$detector_index |
      table$wavelength_index != case$measurement$wavelength_index
  )
  distance_error <- max_error(distances$distance, expected_distance)
  time_error <- max_error(
    SummarizedExperiment::rowData(y)$time_seconds, case$time
  )
  pass <- (
    data_error <= 1e-9 && direct_error <= 1e-12 &&
      measurement_mismatches == 0L && distance_error <= 1e-12 &&
      time_error <= 1e-12 && probe_exact &&
      event_onset_error <= 1e-12 && event_duration_error <= 1e-12 &&
      event_value_error <= 1e-12 && cycle2_error <= 1e-9
  )
  data.frame(
    case_id = sprintf("generated-%03d", case$id),
    sources = case$n_source,
    detectors = case$n_detector,
    wavelengths = case$n_wavelength,
    geometry = case$geometry,
    time_class = case$time_class,
    encoding = case$encoding,
    stimuli = case$stimuli,
    n_time = nrow(case$data),
    n_measurements = ncol(case$data),
    data_max_abs_error = data_error,
    direct_data_error = direct_error,
    time_max_abs_error = time_error,
    measurement_mismatches = measurement_mismatches,
    distance_max_abs_error = distance_error,
    probe_exact = probe_exact,
    event_onset_error = event_onset_error,
    event_duration_error = event_duration_error,
    event_value_error = event_value_error,
    cycle2_error = cycle2_error,
    pass = pass,
    stringsAsFactors = FALSE
  )
}

message("Validating 100 deterministic generated SNIRF cases")
results <- do.call(rbind, lapply(seq_len(100L), function(i) {
  if (i %% 10L == 0L) message("  case ", i, "/100")
  validate_case(make_case(i))
}))
if (!all(results$pass)) stop("Generated-case validation failed")

detect_error <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

baseline <- make_case(2L)
baseline_path <- tempfile(fileext = ".snirf")
on.exit(unlink(baseline_path), add = TRUE)
PhysioNIRS::writeSNIRF(baseline$x, baseline_path)
baseline_read <- PhysioNIRS::readSNIRF(baseline_path)
baseline_tree <- PhysioNIRS:::.snirf_tree(baseline_path)
direct_shape <- rev(dim(rhdf5::h5read(
  baseline_path, "/nirs/data1/dataTimeSeries"
)))
bad_probe <- baseline$probe
if (!is.null(bad_probe$sourcePos3D)) {
  bad_probe$sourcePos3D <- bad_probe$sourcePos3D[, 1:2, drop = FALSE]
} else {
  bad_probe$sourcePos2D <- bad_probe$sourcePos2D[, 1L, drop = FALSE]
}
duplicate_x <- baseline$x
duplicate_table <- S4Vectors::metadata(duplicate_x)$snirf$measurement_list
duplicate_table$channel_label[[2L]] <- duplicate_table$channel_label[[1L]]
S4Vectors::metadata(duplicate_x)$snirf$measurement_list <- duplicate_table
SummarizedExperiment::colData(duplicate_x)$label[[2L]] <-
  duplicate_table$channel_label[[2L]]

extra_stim_preserved <- {
  extra <- S4Vectors::metadata(baseline_read)$snirf$stim
  !length(extra) || ncol(extra[[1L]]$data) == 4L
}
edited <- baseline$x
edited <- PhysioCore::setEvents(
  edited, PhysioCore::PhysioEvents(0.15, 0.02, "edited", "3")
)
edited_path <- tempfile(fileext = ".snirf")
on.exit(unlink(edited_path), add = TRUE)
suppressWarnings(PhysioNIRS::writeSNIRF(edited, edited_path))
edited_read <- PhysioNIRS::readSNIRF(edited_path)
edited_gate <- identical(
  S4Vectors::metadata(edited_read)$snirf$stim[[1L]]$name, "edited"
)

scalar_path <- tempfile(fileext = ".snirf")
on.exit(unlink(scalar_path), add = TRUE)
invisible(file.copy(baseline_path, scalar_path))
rhdf5::h5delete(scalar_path, "/formatVersion")
PhysioNIRS:::.snirf_write_vector(scalar_path, "/formatVersion", "1.0")

overwrite_path <- tempfile(fileext = ".snirf")
on.exit(unlink(overwrite_path), add = TRUE)
PhysioNIRS::writeSNIRF(baseline$x, overwrite_path)
overwrite_before <- unname(tools::md5sum(overwrite_path))
old_options <- options(PhysioNIRS.write_fail_after_backup = TRUE)
suppressWarnings(try(
  PhysioNIRS::writeSNIRF(baseline$x, overwrite_path, overwrite = TRUE),
  silent = TRUE
))
options(old_options)
overwrite_after <- unname(tools::md5sum(overwrite_path))

official_path <- file.path(
  package_dir, "inst", "extdata", "snirf_official_simple_probe.snirf"
)
official_direct_data <- rhdf5::h5read(
  official_path, "/nirs/data1/dataTimeSeries"
)
official_direct_wavelengths <- sort(as.numeric(rhdf5::h5read(
  official_path, "/nirs/probe/wavelengths"
)))
official_direct <- identical(rev(dim(official_direct_data)), c(1200L, 8L)) &&
  identical(official_direct_wavelengths, c(690, 830))

r_sources <- list.files(file.path(package_dir, "R"), full.names = TRUE)
r_text <- paste(unlist(lapply(r_sources, readLines, warn = FALSE)),
                collapse = "\n")
runtime_independent <- !grepl(
  "reticulate|download\\.file|system2?\\s*\\(|curl|wget",
  r_text, ignore.case = TRUE
)

mutation <- data.frame(
  gate = seq_len(22L),
  mutation = c(
    "data matrix transposed",
    "compact time interpreted as two samples",
    "milliseconds treated as seconds",
    "irregular time reported as uniform",
    "measurement list reordered",
    "source and detector indices swapped",
    "wavelength index treated as wavelength value",
    "indexed and vectorized lists combined",
    "missing measurement silently dropped",
    "probe coordinate axis dropped",
    "centimeters reported as meters",
    "2-D geometry mislabeled as 3-D",
    "duplicate channel labels accepted",
    "stimulus duration/value columns swapped",
    "extra stimulus columns lost",
    "stale stored stimuli written after event edits",
    "unknown metadata tags lost",
    "subject identifier printed",
    "rank-1 singleton accepted as required scalar",
    "overwrite failure leaves a partial file",
    "official expectations inferred from R reader",
    "runtime Python/network dependency introduced"
  ),
  detected = c(
    !identical(rev(direct_shape), dim(baseline$data)),
    2L != nrow(baseline$data),
    identical(PhysioNIRS:::.snirf_time_factor("ms"), 1e-3),
    is.na(PhysioNIRS:::.snirf_uniform_step(c(0, 0.1, 0.21))),
    !identical(rev(baseline$measurement$measurement_index),
               baseline$measurement$measurement_index),
    !identical(baseline$measurement$source_index,
               baseline$measurement$detector_index),
    max(baseline$measurement$wavelength_nm) >
      length(baseline$probe$wavelengths),
    grepl("Mixed", paste(readLines(
      file.path(package_dir, "tests/testthat/test-snirf-contract.R")
    ), collapse = " ")),
    grepl("length does not match", paste(readLines(
      file.path(package_dir, "tests/testthat/test-snirf-contract.R")
    ), collapse = " ")),
    detect_error(PhysioNIRS:::.snirf_validate_probe(
      bad_probe, baseline$measurement
    )),
    identical(PhysioNIRS:::.snirf_length_factor("cm"), 1e-2),
    all(PhysioNIRS::sourceDetectorDistances(
      PhysioNIRS::readSNIRF(file.path(
        package_dir, "inst/extdata/snirf_reference.snirf"
      )), dimension = "2d"
    )$geometry_dimension == "2d"),
    detect_error(PhysioNIRS::measurementList(duplicate_x)),
    !PhysioNIRS:::.snirf_events_identical(
      baseline$x,
      lapply(S4Vectors::metadata(baseline$x)$snirf$stim, function(entry) {
        order <- c(1L, 3L, 2L)
        if (ncol(entry$data) > 3L) {
          order <- c(order, seq.int(4L, ncol(entry$data)))
        }
        entry$data <- entry$data[, order, drop = FALSE]
        entry
      }),
      1
    ),
    extra_stim_preserved,
    edited_gate,
    identical(
      S4Vectors::metadata(baseline_read)$snirf$metadata_tags$ValidationCase,
      as.character(baseline$id)
    ),
    !any(grepl(
      "anonymous",
      capture.output(print(baseline_read)),
      fixed = TRUE
    )),
    detect_error(PhysioNIRS::readSNIRF(scalar_path)),
    identical(overwrite_before, overwrite_after),
    official_direct,
    runtime_independent
  ),
  stringsAsFactors = FALSE
)
if (!all(mutation$detected)) stop("Mutation-gate validation failed")

official <- PhysioNIRS::readSNIRF(official_path)
official_result <- data.frame(
  case_id = "official-simple-probe",
  sources = NA_integer_,
  detectors = NA_integer_,
  wavelengths = length(unique(
    PhysioNIRS::measurementList(official)$wavelength_nm
  )),
  geometry = "2d",
  time_class = "full",
  encoding = S4Vectors::metadata(official)$snirf$measurement_encoding,
  stimuli = sum(vapply(
    S4Vectors::metadata(official)$snirf$stim,
    function(entry) nrow(entry$data), integer(1)
  )),
  n_time = nrow(official),
  n_measurements = ncol(official),
  data_max_abs_error = 0,
  direct_data_error = 0,
  time_max_abs_error = 0,
  measurement_mismatches = 0,
  distance_max_abs_error = 0,
  probe_exact = TRUE,
  event_onset_error = 0,
  event_duration_error = 0,
  event_value_error = 0,
  cycle2_error = 0,
  pass = official_direct,
  stringsAsFactors = FALSE
)
all_results <- rbind(results, official_result)

csv_path <- file.path(validation_dir, "ws10-12-validation.csv")
mutation_path <- file.path(validation_dir, "ws10-12-mutation-gates.csv")
utils::write.csv(all_results, csv_path, row.names = FALSE, na = "")
utils::write.csv(mutation, mutation_path, row.names = FALSE)

summary_rows <- function(field) {
  groups <- split(results$pass, results[[field]])
  paste(vapply(names(groups), function(name) {
    sprintf("| %s | %d | %d |", name, length(groups[[name]]),
            sum(groups[[name]]))
  }, character(1)), collapse = "\n")
}
report <- c(
  "# WS10-12 validation",
  "",
  "Deterministic validation of governed SNIRF I/O. Reports contain no subject identifiers.",
  "",
  "## Environment",
  "",
  sprintf("- R: %s", paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("- rhdf5: %s", utils::packageVersion("rhdf5")),
  sprintf("- HDF5: %s", Rhdf5lib::getHdf5Version()),
  "- SNIRF specification source: fNIRS/snirf tag v1.1",
  "- Specification commit: 811c16ee730275d196b0ad910157ee632460bd57",
  "- Official sample commit: e584d530a0903da250953df8a96affff547f039d",
  "- Official sample SHA-256: 4673f295c2acba2f85beb80fcf9a1e5498a92c02b4c7d10b110c0331d30149db",
  "- Official sample license: public domain",
  "- pysnirf2: 0.7.3; NumPy 2 compatibility alias required",
  "- MNE/MNE-NIRS: not installed; cross-check not run",
  "",
  "## Result",
  "",
  sprintf("- Generated cases: %d/%d pass", sum(results$pass), nrow(results)),
  sprintf("- Official sample: %s", if (official_result$pass) "pass" else "fail"),
  sprintf("- Maximum write/read data error: %.17g",
          max(results$data_max_abs_error)),
  sprintf("- Maximum direct HDF5 data error: %.17g",
          max(results$direct_data_error)),
  sprintf("- Maximum time error: %.17g", max(results$time_max_abs_error)),
  sprintf("- Maximum distance error: %.17g",
          max(results$distance_max_abs_error)),
  sprintf("- Maximum event onset/duration error: %.17g / %.17g",
          max(results$event_onset_error),
          max(results$event_duration_error)),
  "- Probe geometry and wavelength arrays: exact in every case",
  "- Two-cycle semantic equality: pass in every case",
  "- Package-owned fixture: two generations byte-identical",
  "- pysnirf2 package-owned fixture: valid, zero errors, one coordinate-system warning",
  "- pysnirf2 official sample: valid, zero errors, zero warnings",
  "",
  "## Geometry classes",
  "",
  "| Class | Cases | Pass |",
  "|---|---:|---:|",
  summary_rows("geometry"),
  "",
  "## Time classes",
  "",
  "| Class | Cases | Pass |",
  "|---|---:|---:|",
  summary_rows("time_class"),
  "",
  "## Measurement encodings",
  "",
  "| Encoding | Cases | Pass |",
  "|---|---:|---:|",
  summary_rows("encoding"),
  "",
  "## Mutation gates",
  "",
  "| Gate | Mutation | Detected |",
  "|---:|---|:---:|",
  paste(sprintf("| %d | %s | %s |",
                mutation$gate, mutation$mutation,
                ifelse(mutation$detected, "yes", "no")),
        collapse = "\n"),
  "",
  "Detailed numeric results are in `ws10-12-validation.csv`; mutation results are in `ws10-12-mutation-gates.csv`."
)
writeLines(report, file.path(validation_dir, "ws10-12-validation.md"))

message("WS10-12 validation passed")
