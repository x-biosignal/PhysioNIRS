make_snirf_experiment <- function(
    n_time = 17L,
    time = seq.int(0, n_time - 1L) / 10,
    geometry = c("2d", "3d"),
    duplicate = FALSE,
    stim = TRUE) {
  wavelengths <- c(760, 850)
  source <- c(1L, 1L, 1L, 2L, 2L, 2L)
  detector <- c(1L, 1L, 2L, 1L, 2L, 2L)
  wavelength <- c(1L, 2L, 1L, 2L, 1L, 2L)
  if (duplicate) {
    source[[6L]] <- source[[5L]]
    detector[[6L]] <- detector[[5L]]
    wavelength[[6L]] <- wavelength[[5L]]
  }
  values <- outer(
    seq_len(n_time), seq_len(6L),
    function(i, j) sin(i / 3) + cos(j / 4) + i * j / 1000
  )
  source_labels <- c("S1", "S2")
  detector_labels <- c("D1", "D2")
  base <- paste0(
    source_labels[source], "_", detector_labels[detector], "_",
    wavelengths[wavelength]
  )
  labels <- base
  repeated <- duplicated(base) | duplicated(base, fromLast = TRUE)
  labels[repeated] <- paste0(
    base[repeated], "_T1I1M", which(repeated)
  )
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(6L),
    source_index = source,
    detector_index = detector,
    wavelength_index = wavelength,
    wavelength_nm = wavelengths[wavelength],
    wavelength_actual_nm = rep(NA_real_, 6L),
    data_type = c(1L, 1L, 1L, 1L, 99999L, 101L),
    data_type_index = rep(1L, 6L),
    data_type_label = c(rep(NA_character_, 4L), "processed", "moments"),
    data_unit = rep(NA_character_, 6L),
    source_power = rep(NA_real_, 6L),
    detector_gain = rep(NA_real_, 6L),
    source_label = source_labels[source],
    detector_label = detector_labels[detector],
    channel_label = labels
  )
  probe <- list(
    wavelengths = wavelengths,
    wavelengthsEmission = c(755, 845),
    sourceLabels = source_labels,
    detectorLabels = detector_labels,
    landmarkLabels = c("Nz", "Cz"),
    coordinateSystem = "Other",
    coordinateSystemDescription = "Synthetic Cartesian fixture",
    LengthUnit = "cm"
  )
  if ("2d" %in% geometry) {
    probe$sourcePos2D <- matrix(c(0, 0, 3, 0), ncol = 2L, byrow = TRUE)
    probe$detectorPos2D <- matrix(c(0, 4, 3, 4), ncol = 2L, byrow = TRUE)
    probe$landmarkPos2D <- matrix(c(0, 5, 3, 5), ncol = 2L, byrow = TRUE)
  }
  if ("3d" %in% geometry) {
    probe$sourcePos3D <- matrix(
      c(0, 0, 0, 3, 0, 0), ncol = 3L, byrow = TRUE
    )
    probe$detectorPos3D <- matrix(
      c(0, 4, 0, 3, 4, 0), ncol = 3L, byrow = TRUE
    )
    probe$landmarkPos3D <- matrix(
      c(0, 5, 0, 3, 5, 0), ncol = 3L, byrow = TRUE
    )
  }
  tags <- list(
    SubjectID = "fixture-subject",
    MeasurementDate = "unknown",
    MeasurementTime = "unknown",
    LengthUnit = "cm",
    TimeUnit = "s",
    FrequencyUnit = "Hz",
    CustomTag = "preserved"
  )
  stimuli <- if (stim) {
    list(
      list(
        name = "left",
        data = matrix(c(0.2, 0.1, 1, 11), nrow = 1L),
        data_labels = c("starttime", "duration", "value", "difficulty")
      ),
      list(
        name = "right",
        data = matrix(c(0.8, 0.2, 2), nrow = 1L),
        data_labels = c("starttime", "duration", "value")
      )
    )
  } else {
    list()
  }
  x <- PhysioCore::PhysioExperiment(
    assays = list(raw = values),
    rowData = S4Vectors::DataFrame(time_seconds = time),
    colData = cbind(measurement, label = measurement$channel_label),
    metadata = list(snirf = list(
      measurement_list = measurement,
      probe = probe,
      metadata_tags = tags,
      stim = stimuli
    )),
    samplingRate = {
      step <- PhysioNIRS:::.snirf_uniform_step(time)
      if (is.na(step)) NA_real_ else 1 / step
    }
  )
  if (stim) {
    x <- PhysioCore::setEvents(
      x,
      PhysioCore::PhysioEvents(
        onset = c(0.2, 0.8), duration = c(0.1, 0.2),
        type = c("left", "right"), value = c("1", "2")
      )
    )
  }
  x
}

write_fixture_file <- function(x = make_snirf_experiment(), compact = FALSE) {
  path <- tempfile(fileext = ".snirf")
  writeSNIRF(x, path, compact_time = compact)
  path
}

convert_to_vectorized <- function(path) {
  table <- measurementList(readSNIRF(path))
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
    dataTypeIndex = as.integer(table$data_type_index),
    dataTypeLabel = ifelse(is.na(table$data_type_label), "",
                           table$data_type_label),
    dataUnit = ifelse(is.na(table$data_unit), "", table$data_unit),
    wavelengthActual = as.numeric(table$wavelength_actual_nm),
    sourcePower = as.numeric(table$source_power),
    detectorGain = as.numeric(table$detector_gain)
  )
  for (name in names(fields)) {
    PhysioNIRS:::.snirf_write_vector(
      path, paste0("/nirs/data1/measurementLists/", name), fields[[name]]
    )
  }
  invisible(path)
}
