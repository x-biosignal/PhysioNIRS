.snirf_read_tags <- function(path, tree, root) {
  group <- paste0(root, "/metaDataTags")
  .snirf_tree_row(tree, group, "H5I_GROUP")
  entries <- .snirf_direct(tree, group)
  groups <- entries[entries$otype == "H5I_GROUP", , drop = FALSE]
  if (nrow(groups)) {
    stop("Nested groups are not supported under metaDataTags: ",
         paste(groups$.path, collapse = ", "), call. = FALSE)
  }
  datasets <- entries[entries$otype == "H5I_DATASET", , drop = FALSE]
  tags <- stats::setNames(vector("list", nrow(datasets)), datasets$name)
  for (i in seq_len(nrow(datasets))) {
    dataset <- datasets$.path[[i]]
    rank <- .snirf_rank(tree, dataset)
    if (rank == 0L) {
      tags[[i]] <- .snirf_read_scalar(path, dataset, tree)
    } else if (rank == 1L) {
      value <- rhdf5::h5read(path, dataset)
      if (is.list(value)) {
        stop("Unsupported metadata dataset type: ", dataset, call. = FALSE)
      }
      tags[[i]] <- as.vector(value)
    } else {
      stop("Metadata tags must be scalar or rank-1: ", dataset, call. = FALSE)
    }
  }
  missing <- setdiff(.snirf_required_tags, names(tags))
  if (length(missing)) {
    stop("Missing required SNIRF metadata tags: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  for (name in .snirf_required_tags) {
    dataset <- paste0(group, "/", name)
    tags[[name]] <- .snirf_read_scalar(path, dataset, tree, "string")
  }
  .snirf_time_factor(tags$TimeUnit)
  tags
}

.snirf_read_probe <- function(path, tree, root, tags) {
  group <- paste0(root, "/probe")
  .snirf_tree_row(tree, group, "H5I_GROUP")
  probe <- list()
  matrix_fields <- c(
    "sourcePos2D", "detectorPos2D", "sourcePos3D", "detectorPos3D",
    "landmarkPos2D", "landmarkPos3D"
  )
  vector_numeric <- c("wavelengths", "wavelengthsEmission")
  vector_string <- c("sourceLabels", "detectorLabels", "landmarkLabels")
  scalar_string <- c("coordinateSystem", "coordinateSystemDescription")
  for (field in .snirf_probe_fields) {
    dataset <- paste0(group, "/", field)
    if (!.snirf_has(tree, dataset, "H5I_DATASET")) next
    probe[[field]] <- if (field %in% matrix_fields) {
      .snirf_read_matrix(path, dataset, tree)
    } else if (field %in% vector_numeric) {
      .snirf_read_vector(path, dataset, tree, "numeric")
    } else if (field %in% vector_string) {
      .snirf_read_vector(path, dataset, tree, "string")
    } else if (field %in% scalar_string) {
      .snirf_read_scalar(path, dataset, tree, "string")
    }
  }
  probe$LengthUnit <- tags$LengthUnit
  .snirf_validate_probe(probe)
}

.snirf_validate_probe <- function(probe, measurement = NULL) {
  if (!is.list(probe)) stop("SNIRF probe metadata must be a list", call. = FALSE)
  complete2 <- all(c("sourcePos2D", "detectorPos2D") %in% names(probe))
  complete3 <- all(c("sourcePos3D", "detectorPos3D") %in% names(probe))
  partial2 <- any(c("sourcePos2D", "detectorPos2D") %in% names(probe))
  partial3 <- any(c("sourcePos3D", "detectorPos3D") %in% names(probe))
  if ((partial2 && !complete2) || (partial3 && !complete3)) {
    stop("Source and detector geometry must be present as a complete pair",
         call. = FALSE)
  }
  if (!complete2 && !complete3) {
    stop("SNIRF probe requires complete 2-D or 3-D geometry", call. = FALSE)
  }
  for (field in c("sourcePos2D", "detectorPos2D")) {
    if (field %in% names(probe) &&
        (!is.matrix(probe[[field]]) || ncol(probe[[field]]) != 2L ||
         anyNA(probe[[field]]) || any(!is.finite(probe[[field]])))) {
      stop("Invalid 2-D probe geometry: ", field, call. = FALSE)
    }
  }
  for (field in c("sourcePos3D", "detectorPos3D")) {
    if (field %in% names(probe) &&
        (!is.matrix(probe[[field]]) || ncol(probe[[field]]) != 3L ||
         anyNA(probe[[field]]) || any(!is.finite(probe[[field]])))) {
      stop("Invalid 3-D probe geometry: ", field, call. = FALSE)
    }
  }
  for (field in c("landmarkPos2D", "landmarkPos3D")) {
    if (!field %in% names(probe)) next
    dimension <- if (field == "landmarkPos2D") 2L else 3L
    if (!is.matrix(probe[[field]]) || ncol(probe[[field]]) != dimension ||
        anyNA(probe[[field]]) || any(!is.finite(probe[[field]]))) {
      stop("Invalid landmark probe geometry: ", field, call. = FALSE)
    }
  }
  for (field in c("wavelengths", "wavelengthsEmission")) {
    if (field %in% names(probe) &&
        (!is.numeric(probe[[field]]) || anyNA(probe[[field]]) ||
         any(!is.finite(probe[[field]])) || any(probe[[field]] <= 0))) {
      stop("Probe wavelength arrays must contain finite positive values: ",
           field, call. = FALSE)
    }
  }
  if ("landmarkLabels" %in% names(probe)) {
    landmark_rows <- unique(vapply(
      c("landmarkPos2D", "landmarkPos3D"),
      function(field) {
        if (is.null(probe[[field]])) NA_integer_ else nrow(probe[[field]])
      },
      integer(1)
    ))
    landmark_rows <- landmark_rows[!is.na(landmark_rows)]
    if (length(landmark_rows) > 1L ||
        (length(landmark_rows) == 1L &&
         length(probe$landmarkLabels) != landmark_rows[[1L]])) {
      stop("Landmark labels must match every landmark geometry row",
           call. = FALSE)
    }
  }
  if (!is.null(measurement)) {
    source_max <- max(measurement$source_index)
    detector_max <- max(measurement$detector_index)
    for (dimension in c("2D", "3D")) {
      source <- probe[[paste0("sourcePos", dimension)]]
      detector <- probe[[paste0("detectorPos", dimension)]]
      if (!is.null(source) && (nrow(source) < source_max ||
                              nrow(detector) < detector_max)) {
        stop("Probe geometry does not cover referenced measurement indices",
             call. = FALSE)
      }
    }
  }
  probe
}

.snirf_measurement_value <- function(path, tree, group, field, required) {
  dataset <- paste0(group, "/", field)
  if (!.snirf_has(tree, dataset, "H5I_DATASET")) {
    if (required) stop("Missing required measurement field: ", dataset,
                       call. = FALSE)
    return(NULL)
  }
  .snirf_read_scalar(
    path, dataset, tree,
    if (field %in% c("dataTypeLabel", "dataUnit")) "string" else "numeric"
  )
}

.snirf_read_measurements <- function(path, tree, data_root, n_measurement,
                                     probe) {
  indexed_rows <- .snirf_direct(tree, data_root, "H5I_GROUP")
  indexed <- indexed_rows$.path[
    grepl("/measurementList[0-9]+$", indexed_rows$.path)
  ]
  vector_root <- paste0(data_root, "/measurementLists")
  vectorized <- .snirf_has(tree, vector_root, "H5I_GROUP")
  if (length(indexed) && vectorized) {
    stop("Mixed indexed and vectorized measurement-list encodings: ",
         paste(c(indexed, vector_root), collapse = ", "), call. = FALSE)
  }
  required <- c(
    "sourceIndex", "detectorIndex", "wavelengthIndex",
    "dataType", "dataTypeIndex"
  )
  rows <- vector("list", n_measurement)
  if (length(indexed)) {
    groups <- .snirf_indexed_groups(tree, data_root, "measurementList")
    if (length(groups) != n_measurement) {
      stop("Measurement-list length does not match data columns", call. = FALSE)
    }
    for (i in seq_along(groups)) {
      values <- stats::setNames(
        vector("list", length(.snirf_measurement_fields)),
        .snirf_measurement_fields
      )
      for (field in .snirf_measurement_fields) {
        values[[field]] <- .snirf_measurement_value(
          path, tree, groups[[i]], field, field %in% required
        )
      }
      rows[[i]] <- values
    }
    encoding <- "indexed"
  } else if (vectorized) {
    for (field in required) {
      dataset <- paste0(vector_root, "/", field)
      if (!.snirf_has(tree, dataset, "H5I_DATASET")) {
        stop("Missing required vectorized measurement field: ", dataset,
             call. = FALSE)
      }
    }
    arrays <- stats::setNames(
      vector("list", length(.snirf_measurement_fields)),
      .snirf_measurement_fields
    )
    for (field in .snirf_measurement_fields) {
      dataset <- paste0(vector_root, "/", field)
      if (!.snirf_has(tree, dataset, "H5I_DATASET")) next
      arrays[[field]] <- .snirf_read_vector(
        path, dataset, tree,
        if (field %in% c("dataTypeLabel", "dataUnit")) "string" else "numeric",
        allow_missing = !field %in% required
      )
      if (length(arrays[[field]]) != n_measurement) {
        stop("Vectorized measurement field has the wrong length: ", dataset,
             call. = FALSE)
      }
    }
    for (i in seq_len(n_measurement)) {
      rows[[i]] <- lapply(arrays, function(value) {
        if (is.null(value)) NULL else value[[i]]
      })
    }
    encoding <- "vectorized"
  } else {
    stop("SNIRF data block has no measurement-list encoding", call. = FALSE)
  }

  get_integer <- function(field) {
    vapply(seq_along(rows), function(i) {
      .snirf_positive_whole(rows[[i]][[field]],
                            paste0(data_root, "/measurement[", i, "]/", field))
    }, integer(1))
  }
  source <- get_integer("sourceIndex")
  detector <- get_integer("detectorIndex")
  wavelength <- get_integer("wavelengthIndex")
  data_type <- get_integer("dataType")
  data_type_index <- get_integer("dataTypeIndex")
  wavelengths <- probe$wavelengths
  if (is.null(wavelengths) || !length(wavelengths)) {
    stop("Probe wavelengths are required for continuous-wave measurements",
         call. = FALSE)
  }
  if (any(wavelength > length(wavelengths))) {
    stop("Measurement wavelengthIndex is outside probe wavelengths",
         call. = FALSE)
  }
  source_count <- max(vapply(
    c("sourcePos2D", "sourcePos3D"),
    function(field) if (is.null(probe[[field]])) 0L else nrow(probe[[field]]),
    integer(1)
  ))
  detector_count <- max(vapply(
    c("detectorPos2D", "detectorPos3D"),
    function(field) if (is.null(probe[[field]])) 0L else nrow(probe[[field]]),
    integer(1)
  ))
  if (any(source > source_count) || any(detector > detector_count)) {
    stop("Measurement source/detector index is outside probe geometry",
         call. = FALSE)
  }
  source_labels <- probe$sourceLabels %||%
    paste0("S", seq_len(source_count))
  detector_labels <- probe$detectorLabels %||%
    paste0("D", seq_len(detector_count))
  if (length(source_labels) < source_count ||
      length(detector_labels) < detector_count) {
    stop("Probe labels do not cover the probe geometry", call. = FALSE)
  }
  nominal <- as.numeric(wavelengths[wavelength])
  base <- paste0(
    source_labels[source], "_", detector_labels[detector], "_",
    format(nominal, trim = TRUE, scientific = FALSE)
  )
  labels <- base
  duplicated_base <- duplicated(base) | duplicated(base, fromLast = TRUE)
  labels[duplicated_base] <- paste0(
    base[duplicated_base], "_T", data_type[duplicated_base],
    "I", data_type_index[duplicated_base],
    "M", which(duplicated_base)
  )
  if (any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("SNIRF channel labels must be complete and unique", call. = FALSE)
  }
  optional_numeric <- function(field) {
    vapply(rows, function(row) {
      value <- row[[field]]
      if (is.null(value)) NA_real_ else {
        value <- as.numeric(value)
        if (length(value) == 1L && is.na(value)) return(NA_real_)
        if (length(value) != 1L || !is.finite(value)) {
          stop("Optional measurement numeric field must be finite: ", field,
               call. = FALSE)
        }
        value
      }
    }, numeric(1))
  }
  optional_string <- function(field) {
    vapply(rows, function(row) {
      value <- row[[field]]
      if (is.null(value)) NA_character_ else as.character(value)
    }, character(1))
  }
  table <- S4Vectors::DataFrame(
    measurement_index = seq_len(n_measurement),
    source_index = source,
    detector_index = detector,
    wavelength_index = wavelength,
    wavelength_nm = nominal,
    wavelength_actual_nm = optional_numeric("wavelengthActual"),
    data_type = data_type,
    data_type_index = data_type_index,
    data_type_label = optional_string("dataTypeLabel"),
    data_unit = optional_string("dataUnit"),
    source_power = optional_numeric("sourcePower"),
    detector_gain = optional_numeric("detectorGain"),
    source_label = as.character(source_labels[source]),
    detector_label = as.character(detector_labels[detector]),
    channel_label = labels
  )
  attr(table, "snirf_encoding") <- encoding
  table
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.snirf_read_stim <- function(path, tree, root, time_factor) {
  groups <- .snirf_indexed_groups(tree, root, "stim", allow_empty = TRUE)
  stim <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    name <- .snirf_read_scalar(path, paste0(group, "/name"), tree, "string")
    data <- .snirf_read_matrix(path, paste0(group, "/data"), tree)
    if (ncol(data) < 3L || any(data[, 2L] < 0)) {
      stop("Stimulus data requires start, nonnegative duration, and value: ",
           group, call. = FALSE)
    }
    labels_path <- paste0(group, "/dataLabels")
    labels <- if (.snirf_has(tree, labels_path, "H5I_DATASET")) {
      .snirf_read_vector(path, labels_path, tree, "string")
    } else {
      c("starttime", "duration", "value",
        if (ncol(data) > 3L) paste0("extra", seq_len(ncol(data) - 3L)))
    }
    if (length(labels) != ncol(data)) {
      stop("Stimulus dataLabels length does not match data columns: ", group,
           call. = FALSE)
    }
    stim[[i]] <- list(name = name, data = data, data_labels = labels)
  }
  stim
}

.snirf_flatten_stim <- function(stim, time_factor) {
  if (!length(stim)) {
    return(data.frame(
      onset = numeric(), duration = numeric(),
      type = character(), value = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, lapply(stim, function(entry) {
    data.frame(
      onset = entry$data[, 1L] * time_factor,
      duration = entry$data[, 2L] * time_factor,
      type = rep(entry$name, nrow(entry$data)),
      value = vapply(entry$data[, 3L], format,
                     character(1), digits = 17, trim = TRUE,
                     scientific = FALSE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$onset), , drop = FALSE]
}

.snirf_event_frame <- function(x) {
  events <- PhysioCore::getEvents(x)
  data <- methods::slot(events, "events")
  frame <- data.frame(
    onset = as.numeric(data$onset),
    duration = as.numeric(data$duration),
    type = as.character(data$type),
    value = as.character(data$value),
    stringsAsFactors = FALSE
  )
  rownames(frame) <- NULL
  frame
}

.snirf_events_identical <- function(x, stim, time_factor) {
  current <- .snirf_event_frame(x)
  stored <- .snirf_flatten_stim(stim, time_factor)
  identical(current, stored)
}

#' Read a Shared Near Infrared Spectroscopy Format file
#'
#' @param path Path to one readable `.snirf` file.
#' @param nirs_index Positive index of the NIRS root.
#' @param data_index Positive index of the data block.
#'
#' @return A [PhysioCore::PhysioExperiment()] with a time-by-measurement
#'   `raw` assay.
#' @export
readSNIRF <- function(path, nirs_index = 1L, data_index = 1L) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !grepl("\\.snirf$", path) ||
      !file.exists(path) || dir.exists(path)) {
    stop("`path` must name one readable regular .snirf file", call. = FALSE)
  }
  if (file.access(path, 4L) != 0L) {
    stop("SNIRF file is not readable: ", path, call. = FALSE)
  }
  nirs_index <- .snirf_index(nirs_index, "nirs_index")
  data_index <- .snirf_index(data_index, "data_index")
  normalized <- normalizePath(path, mustWork = TRUE)
  tree <- .snirf_tree(normalized)
  format_version <- .snirf_read_scalar(
    normalized, "/formatVersion", tree, "string"
  )
  root <- .snirf_select_root(tree, nirs_index)
  data_root <- .snirf_select_data(tree, root, data_index)
  tags <- .snirf_read_tags(normalized, tree, root)
  probe <- .snirf_read_probe(normalized, tree, root, tags)
  data <- .snirf_read_matrix(
    normalized, paste0(data_root, "/dataTimeSeries"), tree
  )
  time_raw <- .snirf_read_vector(
    normalized, paste0(data_root, "/time"), tree, "numeric"
  )
  if (length(time_raw) == 2L && nrow(data) != 2L) {
    if (time_raw[[2L]] <= 0) {
      stop("Compact SNIRF time spacing must be positive", call. = FALSE)
    }
    time_raw <- time_raw[[1L]] +
      seq.int(0, nrow(data) - 1L) * time_raw[[2L]]
  } else {
    if (length(time_raw) != nrow(data)) {
      stop("SNIRF time length does not match data rows", call. = FALSE)
    }
    if (length(time_raw) > 1L && any(diff(time_raw) <= 0)) {
      stop("Full SNIRF time values must be strictly increasing", call. = FALSE)
    }
  }
  time_factor <- .snirf_time_factor(tags$TimeUnit)
  time_seconds <- time_raw * time_factor
  step <- .snirf_uniform_step(time_seconds)
  sampling_rate <- if (is.na(step)) NA_real_ else 1 / step
  measurement <- .snirf_read_measurements(
    normalized, tree, data_root, ncol(data), probe
  )
  probe <- .snirf_validate_probe(probe, measurement)
  stim <- .snirf_read_stim(normalized, tree, root, time_factor)
  col_data <- measurement
  col_data$label <- measurement$channel_label
  metadata <- list(snirf = list(
    format_version = format_version,
    nirs_index = nirs_index,
    data_index = data_index,
    measurement_encoding = attr(measurement, "snirf_encoding"),
    measurement_list = measurement,
    probe = probe,
    metadata_tags = tags,
    stim = stim,
    time_sampling = if (is.na(step)) "irregular" else "uniform"
  ))
  x <- PhysioCore::PhysioExperiment(
    assays = list(raw = data),
    rowData = S4Vectors::DataFrame(time_seconds = time_seconds),
    colData = col_data,
    metadata = metadata,
    samplingRate = sampling_rate,
    provenance = normalized
  )
  flat <- .snirf_flatten_stim(stim, time_factor)
  if (nrow(flat)) {
    x <- PhysioCore::setEvents(
      x,
      PhysioCore::PhysioEvents(
        onset = flat$onset, duration = flat$duration,
        type = flat$type, value = flat$value
      )
    )
  }
  x
}

#' Return the governed SNIRF measurement list
#'
#' @param x A `PhysioExperiment` imported from SNIRF or carrying the same
#'   governed metadata contract.
#'
#' @return A defensive copy of the normalized measurement `DataFrame`.
#' @export
measurementList <- function(x) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  metadata <- S4Vectors::metadata(x)
  table <- metadata$snirf$measurement_list
  required <- c(
    "measurement_index", "source_index", "detector_index",
    "wavelength_index", "wavelength_nm", "wavelength_actual_nm",
    "data_type", "data_type_index", "data_type_label", "data_unit",
    "source_power", "detector_gain", "source_label", "detector_label",
    "channel_label"
  )
  if (!inherits(table, "DataFrame") ||
      !all(required %in% names(table)) ||
      nrow(table) != ncol(x)) {
    stop("`x` has no complete governed SNIRF measurement list",
         call. = FALSE)
  }
  positive_fields <- c(
    "measurement_index", "source_index", "detector_index",
    "wavelength_index", "data_type", "data_type_index"
  )
  for (field in positive_fields) {
    value <- table[[field]]
    if (!is.numeric(value) || length(value) != nrow(table) ||
        anyNA(value) || any(!is.finite(value)) || any(value < 1) ||
        any(value > .Machine$integer.max) || any(value != floor(value))) {
      stop("SNIRF measurement field must contain positive whole numbers: ",
           field, call. = FALSE)
    }
  }
  if (!identical(as.integer(table$measurement_index), seq_len(nrow(table)))) {
    stop("SNIRF measurement_index must preserve exact column order",
         call. = FALSE)
  }
  for (field in c("wavelength_nm")) {
    value <- table[[field]]
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
        any(value <= 0)) {
      stop("SNIRF nominal wavelengths must be finite positive values",
           call. = FALSE)
    }
  }
  for (field in c(
    "wavelength_actual_nm", "source_power", "detector_gain"
  )) {
    value <- table[[field]]
    if (!is.numeric(value) || any(!is.na(value) & !is.finite(value))) {
      stop("Optional SNIRF measurement numeric field is invalid: ", field,
           call. = FALSE)
    }
  }
  for (field in c("source_label", "detector_label")) {
    value <- as.character(table[[field]])
    if (anyNA(value) || any(!nzchar(value))) {
      stop("SNIRF measurement identity labels must be complete: ", field,
           call. = FALSE)
    }
  }
  for (field in c("data_type_label", "data_unit")) {
    value <- as.character(table[[field]])
    if (any(!is.na(value) & !nzchar(value))) {
      stop("Optional SNIRF measurement labels cannot be empty: ", field,
           call. = FALSE)
    }
  }
  labels <- as.character(table$channel_label)
  if (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("SNIRF channel labels must be complete and unique", call. = FALSE)
  }
  col_data <- SummarizedExperiment::colData(x)
  if (!"label" %in% names(col_data) ||
      !identical(labels, as.character(col_data$label))) {
    stop("SNIRF measurement identities disagree with colData labels",
         call. = FALSE)
  }
  S4Vectors::DataFrame(table, check.names = FALSE)
}

.snirf_prepare_time <- function(x, n_time, compact_time) {
  row_data <- SummarizedExperiment::rowData(x)
  if ("time_seconds" %in% names(row_data)) {
    time <- as.numeric(row_data$time_seconds)
    if (length(time) != n_time || anyNA(time) || any(!is.finite(time)) ||
        (length(time) > 1L && any(diff(time) <= 0))) {
      stop("rowData time_seconds must be finite, strictly increasing, and ",
           "match assay rows", call. = FALSE)
    }
  } else {
    rate <- PhysioCore::samplingRate(x)
    if (length(rate) != 1L || is.na(rate) || !is.finite(rate) || rate <= 0) {
      stop("A valid time_seconds column or positive sampling rate is required",
           call. = FALSE)
    }
    time <- seq.int(0, n_time - 1L) / rate
  }
  if (compact_time) {
    step <- .snirf_uniform_step(time)
    if (is.na(step)) {
      stop("`compact_time=TRUE` requires uniformly spaced time values",
           call. = FALSE)
    }
    return(c(time[[1L]], step))
  }
  time
}

.snirf_prepare_tags <- function(snirf, probe) {
  tags <- snirf$metadata_tags %||% list()
  if (!is.list(tags) || is.null(names(tags))) {
    stop("SNIRF metadata_tags must be a named list", call. = FALSE)
  }
  defaults <- list(
    SubjectID = "unknown",
    MeasurementDate = "unknown",
    MeasurementTime = "unknown",
    LengthUnit = probe$LengthUnit,
    TimeUnit = "s",
    FrequencyUnit = "Hz"
  )
  for (name in names(defaults)) {
    if (is.null(tags[[name]]) || !length(tags[[name]])) {
      tags[[name]] <- defaults[[name]]
    }
  }
  tags$LengthUnit <- probe$LengthUnit
  tags$TimeUnit <- "s"
  for (name in names(tags)) {
    value <- tags[[name]]
    if (!is.atomic(value) || is.object(value) || is.null(value) ||
        !length(value) || length(dim(value)) ||
        (!is.numeric(value) && !is.character(value) && !is.logical(value)) ||
        anyNA(value) ||
        (is.numeric(value) && any(!is.finite(value))) ||
        (is.character(value) && any(!nzchar(value)))) {
      stop("SNIRF metadata tag must be one finite atomic scalar or vector: ",
           name, call. = FALSE)
    }
  }
  for (name in .snirf_required_tags) {
    value <- tags[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
      stop("Required SNIRF metadata tag must be a non-empty string: ", name,
           call. = FALSE)
    }
  }
  tags
}

.snirf_prepare_stim <- function(x, snirf) {
  stored <- snirf$stim %||% list()
  old_time_factor <- .snirf_time_factor(
    (snirf$metadata_tags %||% list(TimeUnit = "s"))$TimeUnit %||% "s"
  )
  if (.snirf_events_identical(x, stored, old_time_factor)) {
    return(lapply(stored, function(entry) {
      data <- entry$data
      data[, 1:2] <- data[, 1:2, drop = FALSE] * old_time_factor
      list(name = entry$name, data = data, data_labels = entry$data_labels)
    }))
  }
  current <- .snirf_event_frame(x)
  if (length(stored) && any(vapply(stored, function(entry) {
    ncol(entry$data) > 3L
  }, logical(1)))) {
    warning("PhysioEvents were edited; obsolete SNIRF stimulus extra ",
            "columns were dropped", call. = FALSE)
  } else if (length(stored)) {
    warning("PhysioEvents were edited; stored SNIRF stimulus tables were ",
            "replaced", call. = FALSE)
  }
  if (!nrow(current)) return(list())
  if (anyNA(current) || any(!is.finite(current$onset)) ||
      any(!is.finite(current$duration)) || any(current$duration < 0) ||
      any(!nzchar(current$type))) {
    stop("PhysioEvents cannot be represented as SNIRF stimuli", call. = FALSE)
  }
  amplitude <- suppressWarnings(as.numeric(current$value))
  if (anyNA(amplitude) || any(!is.finite(amplitude))) {
    stop("SNIRF event values must be finite numeric amplitudes",
         call. = FALSE)
  }
  types <- unique(as.character(current$type))
  lapply(types, function(type) {
    keep <- current$type == type
    list(
      name = type,
      data = cbind(
        starttime = current$onset[keep],
        duration = current$duration[keep],
        value = amplitude[keep]
      ),
      data_labels = c("starttime", "duration", "value")
    )
  })
}

.snirf_write_contents <- function(x, path, assay_name, compact_time) {
  data <- SummarizedExperiment::assay(x, assay_name)
  data <- as.matrix(data)
  if (!is.numeric(data) || length(dim(data)) != 2L ||
      any(dim(data) == 0L) || anyNA(data) || any(!is.finite(data))) {
    stop("Selected assay must be a non-empty finite numeric matrix",
         call. = FALSE)
  }
  measurement <- measurementList(x)
  if (nrow(measurement) != ncol(data)) {
    stop("Measurement-list rows must match assay columns", call. = FALSE)
  }
  snirf <- S4Vectors::metadata(x)$snirf
  probe <- .snirf_validate_probe(snirf$probe, measurement)
  wavelength_index <- as.integer(measurement$wavelength_index)
  if (is.null(probe$wavelengths) ||
      any(wavelength_index > length(probe$wavelengths)) ||
      !identical(
        as.numeric(measurement$wavelength_nm),
        as.numeric(probe$wavelengths[wavelength_index])
      )) {
    stop("Measurement wavelengths disagree with the governed probe lookup",
         call. = FALSE)
  }
  if (!is.null(probe$sourceLabels) &&
      !identical(
        as.character(measurement$source_label),
        as.character(probe$sourceLabels[measurement$source_index])
      )) {
    stop("Measurement source labels disagree with the governed probe lookup",
         call. = FALSE)
  }
  if (!is.null(probe$detectorLabels) &&
      !identical(
        as.character(measurement$detector_label),
        as.character(probe$detectorLabels[measurement$detector_index])
      )) {
    stop("Measurement detector labels disagree with the governed probe lookup",
         call. = FALSE)
  }
  tags <- .snirf_prepare_tags(snirf, probe)
  time <- .snirf_prepare_time(x, nrow(data), compact_time)
  stim <- .snirf_prepare_stim(x, snirf)

  if (!rhdf5::h5createFile(path)) {
    stop("Unable to create temporary SNIRF file", call. = FALSE)
  }
  .snirf_write_scalar(path, "/formatVersion", "1.0")
  for (group in c(
    "/nirs", "/nirs/metaDataTags", "/nirs/data1", "/nirs/probe"
  )) {
    rhdf5::h5createGroup(path, group)
  }
  .snirf_write_matrix(path, "/nirs/data1/dataTimeSeries", data)
  .snirf_write_vector(path, "/nirs/data1/time", as.numeric(time))
  for (i in seq_len(nrow(measurement))) {
    group <- paste0("/nirs/data1/measurementList", i)
    rhdf5::h5createGroup(path, group)
    values <- list(
      sourceIndex = as.integer(measurement$source_index[[i]]),
      detectorIndex = as.integer(measurement$detector_index[[i]]),
      wavelengthIndex = as.integer(measurement$wavelength_index[[i]]),
      dataType = as.integer(measurement$data_type[[i]]),
      dataTypeIndex = as.integer(measurement$data_type_index[[i]]),
      dataTypeLabel = measurement$data_type_label[[i]],
      dataUnit = measurement$data_unit[[i]],
      wavelengthActual = measurement$wavelength_actual_nm[[i]],
      sourcePower = measurement$source_power[[i]],
      detectorGain = measurement$detector_gain[[i]]
    )
    for (field in names(values)) {
      value <- values[[field]]
      if (length(value) && !is.na(value)) {
        .snirf_write_scalar(path, paste0(group, "/", field), value)
      }
    }
  }
  for (field in .snirf_probe_fields) {
    value <- probe[[field]]
    if (is.null(value) || !length(value)) next
    dataset <- paste0("/nirs/probe/", field)
    if (is.matrix(value)) {
      .snirf_write_matrix(path, dataset, value)
    } else if (length(value) == 1L &&
               field %in% c("coordinateSystem",
                            "coordinateSystemDescription")) {
      .snirf_write_scalar(path, dataset, as.character(value))
    } else {
      .snirf_write_vector(path, dataset, value)
    }
  }
  for (name in names(tags)) {
    value <- tags[[name]]
    dataset <- paste0("/nirs/metaDataTags/", name)
    if (length(value) == 1L) {
      .snirf_write_scalar(path, dataset, value)
    } else {
      .snirf_write_vector(path, dataset, value)
    }
  }
  for (i in seq_along(stim)) {
    group <- paste0("/nirs/stim", i)
    rhdf5::h5createGroup(path, group)
    .snirf_write_scalar(path, paste0(group, "/name"), stim[[i]]$name)
    .snirf_write_matrix(path, paste0(group, "/data"), stim[[i]]$data)
    .snirf_write_vector(
      path, paste0(group, "/dataLabels"), stim[[i]]$data_labels
    )
  }
  invisible(list(data = data, measurement = measurement))
}

.snirf_low_level_validate <- function(path, expected) {
  tree <- .snirf_tree(path)
  scalar_paths <- c(
    "/formatVersion",
    paste0("/nirs/metaDataTags/", .snirf_required_tags)
  )
  for (dataset in scalar_paths) {
    if (.snirf_rank(tree, dataset) != 0L) {
      stop("Temporary SNIRF contains a non-scalar required field: ", dataset,
           call. = FALSE)
    }
  }
  data <- .snirf_read_matrix(
    path, "/nirs/data1/dataTimeSeries", tree
  )
  if (!identical(dim(data), dim(expected$data)) ||
      max(abs(data - expected$data)) > 1e-12) {
    stop("Temporary SNIRF data failed independent orientation validation",
         call. = FALSE)
  }
  groups <- .snirf_indexed_groups(
    tree, "/nirs/data1", "measurementList"
  )
  if (length(groups) != nrow(expected$measurement)) {
    stop("Temporary SNIRF measurement count failed validation",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Write a Shared Near Infrared Spectroscopy Format file
#'
#' @param x A governed `PhysioExperiment`.
#' @param path Destination ending exactly in `.snirf`.
#' @param assay_name Assay to write. `NULL` uses the default assay.
#' @param overwrite Whether to replace an existing destination.
#' @param compact_time Whether to use `[start, spacing]` time encoding.
#'
#' @return `path`, invisibly.
#' @export
writeSNIRF <- function(x, path, assay_name = NULL, overwrite = FALSE,
                       compact_time = FALSE) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !grepl("\\.snirf$", path)) {
    stop("`path` must end exactly in .snirf", call. = FALSE)
  }
  overwrite <- .snirf_flag(overwrite, "overwrite")
  compact_time <- .snirf_flag(compact_time, "compact_time")
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop("Destination parent directory does not exist: ", parent,
         call. = FALSE)
  }
  if (file.exists(path) && !overwrite) {
    stop("Destination exists; use `overwrite=TRUE` to replace it",
         call. = FALSE)
  }
  assay_names <- SummarizedExperiment::assayNames(x)
  if (is.null(assay_name)) {
    assay_name <- PhysioCore::defaultAssay(x)
  } else if (!is.character(assay_name) || length(assay_name) != 1L ||
             is.na(assay_name) || !assay_name %in% assay_names) {
    stop("`assay_name` must exactly identify one assay", call. = FALSE)
  }
  if (length(assay_name) != 1L || is.na(assay_name) ||
      !assay_name %in% assay_names) {
    stop("No writable default assay is available", call. = FALSE)
  }
  destination <- normalizePath(parent, mustWork = TRUE)
  destination <- file.path(destination, basename(path))
  temp <- tempfile(pattern = ".physionirs-", tmpdir = parent,
                   fileext = ".snirf")
  backup <- NULL
  installed <- FALSE
  on.exit({
    if (file.exists(temp)) unlink(temp)
    if (!installed && !is.null(backup) && file.exists(backup) &&
        !file.exists(destination)) {
      file.rename(backup, destination)
    }
    if (installed && !is.null(backup) && file.exists(backup)) unlink(backup)
  }, add = TRUE)
  expected <- .snirf_write_contents(
    x, temp, assay_name = assay_name, compact_time = compact_time
  )
  .snirf_low_level_validate(temp, expected)
  if (isTRUE(getOption("PhysioNIRS.write_fail_after_validation", FALSE))) {
    stop("Injected SNIRF write failure after validation", call. = FALSE)
  }
  if (file.exists(destination)) {
    backup <- tempfile(pattern = ".physionirs-backup-", tmpdir = parent,
                       fileext = ".snirf")
    if (!file.rename(destination, backup)) {
      stop("Unable to move existing SNIRF destination to backup",
           call. = FALSE)
    }
  }
  if (isTRUE(getOption("PhysioNIRS.write_fail_after_backup", FALSE))) {
    stop("Injected SNIRF write failure after backup", call. = FALSE)
  }
  if (!file.rename(temp, destination)) {
    stop("Unable to atomically install SNIRF destination", call. = FALSE)
  }
  installed <- TRUE
  invisible(destination)
}
