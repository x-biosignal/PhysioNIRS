make_cw_nirs <- function(
    wavelengths = c(760, 850),
    distances = c(0.03, 0.04),
    n_time = 41L,
    assay_name = "raw",
    values = NULL) {
  n_pair <- length(distances)
  stopifnot(n_pair >= 1L, length(wavelengths) >= 2L)
  source_index <- rep(seq_len(n_pair), each = length(wavelengths))
  detector_index <- source_index
  wavelength_index <- rep(seq_along(wavelengths), times = n_pair)
  n_measurement <- length(source_index)
  source_labels <- paste0("S", seq_len(n_pair))
  detector_labels <- paste0("D", seq_len(n_pair))
  channel_label <- paste0(
    source_labels[source_index], "_",
    detector_labels[detector_index], "_",
    wavelengths[wavelength_index]
  )
  repeated <- duplicated(channel_label) |
    duplicated(channel_label, fromLast = TRUE)
  channel_label[repeated] <- paste0(
    channel_label[repeated], "_M", which(repeated)
  )
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(n_measurement),
    source_index = as.integer(source_index),
    detector_index = as.integer(detector_index),
    wavelength_index = as.integer(wavelength_index),
    wavelength_nm = as.numeric(wavelengths[wavelength_index]),
    wavelength_actual_nm = rep(NA_real_, n_measurement),
    data_type = rep(1L, n_measurement),
    data_type_index = rep(1L, n_measurement),
    data_type_label = rep(NA_character_, n_measurement),
    data_unit = rep(NA_character_, n_measurement),
    source_power = rep(NA_real_, n_measurement),
    detector_gain = rep(NA_real_, n_measurement),
    source_label = source_labels[source_index],
    detector_label = detector_labels[detector_index],
    channel_label = channel_label
  )
  source_position <- cbind(seq_len(n_pair) - 1, 0)
  detector_position <- source_position
  detector_position[, 2L] <- distances
  probe <- list(
    wavelengths = as.numeric(wavelengths),
    sourceLabels = source_labels,
    detectorLabels = detector_labels,
    sourcePos2D = source_position,
    detectorPos2D = detector_position,
    LengthUnit = "m"
  )
  if (is.null(values)) {
    time <- seq.int(0, n_time - 1L) / 10
    values <- outer(
      time,
      seq_len(n_measurement),
      function(t, j) 1000 + 20 * sin(2 * pi * 0.1 * t + j / 5)
    )
  } else {
    values <- as.matrix(values)
    n_time <- nrow(values)
    stopifnot(ncol(values) == n_measurement)
  }
  time <- seq.int(0, n_time - 1L) / 10
  dimnames(values) <- list(NULL, channel_label)
  metadata <- list(snirf = list(
    measurement_list = measurement,
    probe = probe,
    metadata_tags = list(
      SubjectID = "private-fixture",
      MeasurementDate = "unknown",
      MeasurementTime = "unknown",
      LengthUnit = "m",
      TimeUnit = "s",
      FrequencyUnit = "Hz"
    ),
    stim = list(),
    time_sampling = "uniform"
  ))
  if (assay_name == "OD") {
    metadata$nirs <- list(assays = list(OD = list(
      kind = "optical_density",
      unit = "1",
      log_convention = "natural",
      source_assay = "synthetic"
    )))
  }
  x <- PhysioCore::PhysioExperiment(
    assays = stats::setNames(list(values), assay_name),
    rowData = S4Vectors::DataFrame(time_seconds = time),
    colData = cbind(measurement, label = channel_label),
    metadata = metadata,
    samplingRate = 10
  )
  PhysioCore::setEvents(
    x,
    PhysioCore::PhysioEvents(
      onset = 1, duration = 0.5, type = "task", value = "1"
    )
  )
}

forward_mbll_od <- function(concentration_uM, wavelengths, distance, ppf) {
  table <- utils::read.csv(
    system.file(
      "extdata", "haemoglobin-extinction.csv", package = "PhysioNIRS"
    ),
    stringsAsFactors = FALSE
  )
  interpolate <- function(column) {
    stats::approx(
      table$wavelength_nm, table[[column]], wavelengths,
      method = "linear", rule = 1, ties = "ordered"
    )$y
  }
  extinction <- cbind(
    HbO = interpolate("hbo2_cm1_M1"),
    HbR = interpolate("hbr_cm1_M1")
  ) * (100 * log(10))
  concentration_M <- concentration_uM * 1e-6
  concentration_M %*% t(extinction * (distance * ppf))
}
