.nirs_short_reference <- list(
  mne_nirs_commit = "0a5081735144b902a3953e81d010420e1210c556",
  short_source_sha256 =
    "c9c0f621301e6c0aafca29bb9e987e2d455cf06c59ae1f8e3418e6c39f88991f",
  correction_source_sha256 =
    "0d8c34bd559e49839ff41a8b3699a556bb4efb0b17bf6ecc2c0867b875c406c4"
)

.nirs_short_number_id <- function(x) {
  sprintf("%.15g", as.numeric(x))
}

.nirs_short_identity <- function(x) {
  metadata <- S4Vectors::metadata(x)
  is_mbll <- is.list(metadata$nirs$mbll)
  if (!is_mbll) {
    measurement <- measurementList(x)
    channel_id <- paste0(
      "S", as.integer(measurement$source_index),
      "_D", as.integer(measurement$detector_index),
      "_wl", .nirs_short_number_id(measurement$wavelength_nm)
    )
    if (anyDuplicated(channel_id)) {
      stop(
        "Governed measurements contain duplicate source-detector-wavelength ",
        "identity",
        call. = FALSE
      )
    }
    return(data.frame(
      channel_index = seq_len(ncol(x)),
      channel_id = channel_id,
      source_index = as.integer(measurement$source_index),
      detector_index = as.integer(measurement$detector_index),
      wavelength_nm = as.numeric(measurement$wavelength_nm),
      identity_kind = rep("measurement", ncol(x)),
      stringsAsFactors = FALSE
    ))
  }

  columns <- SummarizedExperiment::colData(x)
  required <- c(
    "channel_id", "source_index", "detector_index",
    "source_detector_distance_m"
  )
  if (!all(required %in% names(columns)) || nrow(columns) != ncol(x)) {
    stop(
      "`x` must carry governed SNIRF measurement or MBLL pair identity",
      call. = FALSE
    )
  }
  channel_id <- as.character(columns$channel_id)
  source <- as.numeric(columns$source_index)
  detector <- as.numeric(columns$detector_index)
  stored_distance <- as.numeric(columns$source_detector_distance_m)
  if (anyNA(channel_id) || any(!nzchar(channel_id)) ||
      anyDuplicated(channel_id) ||
      anyNA(source) || any(!is.finite(source)) ||
      any(source < 1 | source != floor(source)) ||
      anyNA(detector) || any(!is.finite(detector)) ||
      any(detector < 1 | detector != floor(detector)) ||
      anyNA(stored_distance) || any(!is.finite(stored_distance)) ||
      any(stored_distance <= 0)) {
    stop("Governed MBLL pair identity is malformed", call. = FALSE)
  }
  pair_key <- paste(source, detector, sep = ":")
  if (anyDuplicated(pair_key)) {
    stop("Governed MBLL source-detector pairs must be unique", call. = FALSE)
  }
  data.frame(
    channel_index = seq_len(ncol(x)),
    channel_id = channel_id,
    source_index = as.integer(source),
    detector_index = as.integer(detector),
    wavelength_nm = rep(NA_real_, ncol(x)),
    identity_kind = rep("pair", ncol(x)),
    stored_distance_m = stored_distance,
    stringsAsFactors = FALSE
  )
}

.nirs_short_probe_geometry <- function(x, identity) {
  metadata <- S4Vectors::metadata(x)
  probe <- .snirf_validate_probe(metadata$snirf$probe, identity)
  has_3d <- all(c("sourcePos3D", "detectorPos3D") %in% names(probe))
  has_2d <- all(c("sourcePos2D", "detectorPos2D") %in% names(probe))
  dimension <- if (has_3d) 3L else if (has_2d) 2L else {
    stop("Governed probe has no complete source-detector geometry",
         call. = FALSE)
  }
  suffix <- if (dimension == 3L) "3D" else "2D"
  factor <- .snirf_length_factor(probe$LengthUnit)
  source_all <- probe[[paste0("sourcePos", suffix)]] * factor
  detector_all <- probe[[paste0("detectorPos", suffix)]] * factor
  source <- source_all[identity$source_index, , drop = FALSE]
  detector <- detector_all[identity$detector_index, , drop = FALSE]
  distance <- sqrt(rowSums((source - detector)^2))
  if (anyNA(distance) || any(!is.finite(distance)) || any(distance <= 0)) {
    stop("Probe geometry produced non-positive or non-finite distances",
         call. = FALSE)
  }
  if ("stored_distance_m" %in% names(identity)) {
    stored <- identity$stored_distance_m
    tolerance <- pmax(1e-9, 1e-6 * pmax(distance, 1))
    bad <- which(abs(stored - distance) > tolerance)
    if (length(bad)) {
      stop(
        "Stored source-detector distance disagrees with probe coordinates for ",
        identity$channel_id[[bad[[1L]]]],
        call. = FALSE
      )
    }
  }
  midpoint <- (source + detector) / 2
  if (dimension == 2L) {
    midpoint <- cbind(midpoint, 0)
  }
  list(
    distance_m = as.numeric(distance),
    midpoint_m = unname(midpoint),
    dimension = dimension,
    fingerprint = .nirs_sha256(list(
      dimension = dimension,
      length_unit = probe$LengthUnit,
      source_position_m = unname(source_all),
      detector_position_m = unname(detector_all)
    ))
  )
}

.nirs_short_all_assays <- function(x) {
  names <- SummarizedExperiment::assayNames(x)
  values <- lapply(names, function(name) {
    value <- .nirs_assay_matrix(x, name)
    if (anyNA(value) || any(!is.finite(value))) {
      stop("All source assays must contain only finite values", call. = FALSE)
    }
    list(value = unname(value), dimnames = dimnames(value))
  })
  stats::setNames(values, names)
}

.nirs_short_context <- function(
    x, assay_name = NULL, min_samples = 2L, require_contract = FALSE) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  identity <- .nirs_short_identity(x)
  geometry <- .nirs_short_probe_geometry(x, identity)
  time <- .nirs_time_seconds(x)
  if (length(time) < min_samples) {
    stop("`x` has too few time samples", call. = FALSE)
  }
  delta <- diff(time)
  step <- stats::median(delta)
  tolerance <- max(1e-9, 1e-6 * step)
  if (!is.finite(step) || step <= 0 ||
      any(abs(delta - step) > tolerance)) {
    stop("`x` time must be uniformly sampled", call. = FALSE)
  }
  fs <- 1 / step
  recorded_fs <- suppressWarnings(as.numeric(PhysioCore::samplingRate(x)))
  if (length(recorded_fs) != 1L || !is.finite(recorded_fs) ||
      recorded_fs <= 0 ||
      abs(recorded_fs - fs) > max(1e-9, 1e-6 * fs)) {
    stop("Recorded sampling rate disagrees with the governed time base",
         call. = FALSE)
  }
  fs <- as.numeric(recorded_fs)
  assays <- .nirs_short_all_assays(x)
  data <- NULL
  contract <- NULL
  if (!is.null(assay_name)) {
    assay_name <- .nirs_scalar_name(assay_name, "assay_name")
    data <- .nirs_assay_matrix(x, assay_name)
    contract <- S4Vectors::metadata(x)$nirs$assays[[assay_name]]
    if (require_contract) {
      if (!is.list(contract) || !is.character(contract$kind) ||
          length(contract$kind) != 1L || is.na(contract$kind)) {
        stop("Assay `", assay_name, "` has no complete governed contract",
             call. = FALSE)
      }
      if (identity$identity_kind[[1L]] == "measurement" &&
          contract$kind != "optical_density") {
        stop("Measurement-level short separation requires optical density",
             call. = FALSE)
      }
      if (identity$identity_kind[[1L]] == "pair" &&
          contract$kind != "haemoglobin_concentration") {
        stop("Pair-level short separation requires haemoglobin concentration",
             call. = FALSE)
      }
    }
  }
  identity_fingerprint <- .nirs_sha256(identity)
  source_fingerprint <- .nirs_sha256(list(
    assays = assays,
    assay_contracts = S4Vectors::metadata(x)$nirs$assays,
    time = time,
    identity = identity,
    probe_fingerprint = geometry$fingerprint
  ))
  list(
    assay_name = assay_name,
    data = data,
    contract = contract,
    time = time,
    fs = fs,
    identity = identity,
    geometry = geometry,
    identity_fingerprint = identity_fingerprint,
    source_fingerprint = source_fingerprint
  )
}

.nirs_make_short_map <- function(context, threshold_m) {
  midpoint <- context$geometry$midpoint_m
  out <- data.frame(
    channel_index = as.integer(context$identity$channel_index),
    channel_id = as.character(context$identity$channel_id),
    source_index = as.integer(context$identity$source_index),
    detector_index = as.integer(context$identity$detector_index),
    distance_m = as.numeric(context$geometry$distance_m),
    midpoint_x_m = as.numeric(midpoint[, 1L]),
    midpoint_y_m = as.numeric(midpoint[, 2L]),
    midpoint_z_m = as.numeric(midpoint[, 3L]),
    is_short = as.logical(context$geometry$distance_m < threshold_m),
    identity_kind = as.character(context$identity$identity_kind),
    stringsAsFactors = FALSE
  )
  class(out) <- c("nirs_short_channels", "data.frame")
  attr(out, "schema_version") <- 1L
  attr(out, "threshold_m") <- threshold_m
  attr(out, "source_fingerprint") <- context$source_fingerprint
  attr(out, "identity_fingerprint") <- context$identity_fingerprint
  attr(out, "probe_fingerprint") <- context$geometry$fingerprint
  out
}

.nirs_validate_short_map <- function(short, x, context = NULL) {
  required <- c(
    "channel_index", "channel_id", "source_index", "detector_index",
    "distance_m", "midpoint_x_m", "midpoint_y_m", "midpoint_z_m",
    "is_short", "identity_kind"
  )
  threshold <- attr(short, "threshold_m", exact = TRUE)
  if (!inherits(short, "nirs_short_channels") ||
      !is.data.frame(short) || !identical(names(short), required) ||
      !is.numeric(threshold) || !is.null(dim(threshold)) ||
      length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold) || threshold <= 0) {
    stop("`short` must be a complete `nirs_short_channels` map",
         call. = FALSE)
  }
  if (is.null(context)) context <- .nirs_short_context(x)
  expected <- .nirs_make_short_map(context, as.numeric(threshold))
  columns_equal <- vapply(
    required,
    function(name) identical(short[[name]], expected[[name]]),
    logical(1)
  )
  attrs <- c(
    "schema_version", "threshold_m", "source_fingerprint",
    "identity_fingerprint", "probe_fingerprint"
  )
  attrs_equal <- vapply(
    attrs,
    function(name) identical(
      attr(short, name, exact = TRUE),
      attr(expected, name, exact = TRUE)
    ),
    logical(1)
  )
  expected_attrs <- sort(names(attributes(expected)))
  actual_attrs <- sort(names(attributes(short)))
  if (!all(columns_equal) || !all(attrs_equal) ||
      !identical(nrow(short), nrow(expected)) ||
      !identical(actual_attrs, expected_attrs)) {
    stop("`short` is stale, malformed, or belongs to another object",
         call. = FALSE)
  }
  short
}

#' Identify short source-detector channels
#'
#' Calculates governed source-detector distances and midpoints in metres.
#' Channels strictly below `threshold_m` are classified as short. A returned
#' map is fingerprint-bound to the complete source object.
#'
#' @param x A governed measurement-level or MBLL pair-level
#'   `PhysioExperiment`.
#' @param threshold_m Positive short-separation threshold in metres.
#'
#' @return A source-bound `nirs_short_channels` data frame.
#' @export
identifyShortChannels <- function(x, threshold_m = 0.01) {
  threshold_m <- .nirs_motion_scalar(
    threshold_m, "threshold_m", lower = .Machine$double.xmin
  )
  context <- .nirs_short_context(x)
  .nirs_make_short_map(context, threshold_m)
}
