.nirs_motion_reference <- list(
  homer3_commit = "a2bdfcf65e932478110cd9abdd4f0d1b773c5217",
  tddr_commit = "2b104674fdf39027f5148d7d97f61b60bad9327c",
  csaps_commit = "4c1d003e822a3432cd52cd9e5a6c9662e966d0c9"
)

.nirs_motion_scalar <- function(
    value, arg, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 1L ||
      is.na(value) || !is.finite(value) || value < lower || value > upper ||
      (integer && value != floor(value))) {
    qualifier <- if (integer) "integer " else ""
    stop(
      "`", arg, "` must be one finite ", qualifier,
      "value in [", lower, ", ", upper, "]",
      call. = FALSE
    )
  }
  if (integer) as.integer(value) else as.numeric(value)
}

.nirs_motion_flag <- function(value, arg) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", arg, "` must be TRUE or FALSE", call. = FALSE)
  }
  value
}

.nirs_sha256 <- function(value) {
  digest::digest(
    value,
    algo = "sha256",
    serialize = TRUE,
    serializeVersion = 2
  )
}

.nirs_motion_identity <- function(x) {
  measurement <- tryCatch(measurementList(x), error = function(e) NULL)
  if (!is.null(measurement)) {
    return(data.frame(
      measurement_index = seq_len(ncol(x)),
      measurement_id = as.character(measurement$channel_label),
      source_index = as.integer(measurement$source_index),
      detector_index = as.integer(measurement$detector_index),
      wavelength_nm = as.numeric(measurement$wavelength_nm),
      stringsAsFactors = FALSE
    ))
  }

  col_data <- SummarizedExperiment::colData(x)
  required <- c("channel_id", "source_index", "detector_index")
  if (!all(required %in% names(col_data)) || nrow(col_data) != ncol(x)) {
    stop(
      "`x` must carry governed SNIRF measurement or MBLL pair identity",
      call. = FALSE
    )
  }
  id <- as.character(col_data$channel_id)
  source <- as.numeric(col_data$source_index)
  detector <- as.numeric(col_data$detector_index)
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id) ||
      anyNA(source) || anyNA(detector) ||
      any(!is.finite(source)) || any(!is.finite(detector)) ||
      any(source < 1 | source != floor(source)) ||
      any(detector < 1 | detector != floor(detector))) {
    stop("Governed MBLL pair identity is malformed", call. = FALSE)
  }
  data.frame(
    measurement_index = seq_len(ncol(x)),
    measurement_id = id,
    source_index = as.integer(source),
    detector_index = as.integer(detector),
    wavelength_nm = rep(NA_real_, ncol(x)),
    stringsAsFactors = FALSE
  )
}

.nirs_motion_context <- function(x, assay_name, min_samples = 2L) {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  assay_name <- .nirs_scalar_name(assay_name, "assay_name")
  data <- .nirs_assay_matrix(x, assay_name)
  if (anyNA(data) || any(!is.finite(data))) {
    bad <- which(!is.finite(data), arr.ind = TRUE)[1L, ]
    stop(
      "Assay `", assay_name, "` row ", bad[[1L]], ", channel ",
      bad[[2L]], " is non-finite",
      call. = FALSE
    )
  }
  min_samples <- as.integer(min_samples)
  if (nrow(data) < min_samples) {
    stop(
      "Assay `", assay_name, "` requires at least ", min_samples,
      " samples; found ", nrow(data),
      call. = FALSE
    )
  }
  time <- .nirs_time_seconds(x)
  delta <- diff(time)
  step <- stats::median(delta)
  tolerance <- max(1e-9, 1e-6 * step)
  if (!is.finite(step) || step <= 0 ||
      any(abs(delta - step) > tolerance)) {
    stop(
      "`x` time must be uniformly sampled within tolerance ",
      format(tolerance, scientific = TRUE),
      call. = FALSE
    )
  }
  fs <- 1 / step
  recorded_fs <- suppressWarnings(as.numeric(PhysioCore::samplingRate(x)))
  if (length(recorded_fs) == 1L && is.finite(recorded_fs) &&
      recorded_fs > 0 &&
      abs(recorded_fs - fs) > max(1e-9, 1e-6 * fs)) {
    stop(
      "Recorded sampling rate disagrees with `rowData(x)$time_seconds`",
      call. = FALSE
    )
  }
  identity <- .nirs_motion_identity(x)
  if (nrow(identity) != ncol(data) ||
      !identical(identity$measurement_index, seq_len(ncol(data)))) {
    stop("Governed channel identity is not aligned to assay columns",
         call. = FALSE)
  }
  list(
    assay_name = assay_name,
    data = data,
    time = time,
    fs = fs,
    identity = identity,
    fingerprint = .nirs_sha256(list(
      assay_name = assay_name,
      data = unname(data),
      dimnames = dimnames(data),
      time = time,
      identity = identity
    ))
  )
}

.nirs_motion_intervals <- function(mask, time, identity) {
  empty <- data.frame(
    measurement_index = integer(),
    start_index = integer(),
    end_index = integer(),
    start_time = numeric(),
    end_time = numeric()
  )
  rows <- vector("list", 0L)
  k <- 0L
  for (j in seq_len(ncol(mask))) {
    runs <- rle(mask[, j])
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    bad <- which(runs$values)
    for (i in bad) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        measurement_index = as.integer(identity$measurement_index[[j]]),
        start_index = as.integer(starts[[i]]),
        end_index = as.integer(ends[[i]]),
        start_time = as.numeric(time[[starts[[i]]]]),
        end_time = as.numeric(time[[ends[[i]]]])
      )
    }
  }
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

.nirs_motion_mask_fingerprint <- function(mask) {
  .nirs_sha256(list(
    sample_by_measurement = unname(mask$sample_by_measurement),
    measurement_id = mask$measurement_id,
    assay_name = mask$assay_name,
    sampling_rate_hz = mask$sampling_rate_hz,
    input_fingerprint = mask$input_fingerprint,
    parameters = mask$parameters
  ))
}

.nirs_validate_motion_mask <- function(mask, context) {
  required <- c(
    "sample_by_measurement", "global", "intervals", "parameters",
    "assay_name", "sampling_rate_hz", "measurement_id",
    "input_fingerprint", "fingerprint"
  )
  if (!inherits(mask, "nirs_motion_mask") || !is.list(mask) ||
      !all(required %in% names(mask))) {
    stop("`mask` must be a complete `nirs_motion_mask`", call. = FALSE)
  }
  matrix <- mask$sample_by_measurement
  if (!is.matrix(matrix) || !is.logical(matrix) || anyNA(matrix) ||
      !identical(dim(matrix), dim(context$data))) {
    stop("`mask$sample_by_measurement` is malformed or misaligned",
         call. = FALSE)
  }
  if (!is.logical(mask$global) || anyNA(mask$global) ||
      length(mask$global) != nrow(context$data) ||
      !identical(mask$global, apply(matrix, 1L, any))) {
    stop("`mask$global` is malformed or inconsistent", call. = FALSE)
  }
  expected_ids <- as.character(context$identity$measurement_id)
  if (!identical(mask$measurement_id, expected_ids) ||
      !identical(mask$assay_name, context$assay_name) ||
      !is.numeric(mask$sampling_rate_hz) ||
      length(mask$sampling_rate_hz) != 1L ||
      !isTRUE(all.equal(
        as.numeric(mask$sampling_rate_hz),
        context$fs,
        tolerance = 1e-12
      )) ||
      !identical(mask$input_fingerprint, context$fingerprint)) {
    stop("`mask` is stale or belongs to another assay/object", call. = FALSE)
  }
  expected_intervals <- .nirs_motion_intervals(
    matrix, context$time, context$identity
  )
  if (!identical(mask$intervals, expected_intervals)) {
    stop("`mask$intervals` is inconsistent with its logical matrix",
         call. = FALSE)
  }
  if (!identical(mask$fingerprint, .nirs_motion_mask_fingerprint(mask))) {
    stop("`mask` fingerprint is invalid", call. = FALSE)
  }
  mask
}

.nirs_add_motion_assay <- function(
    x, context, output_assay, values, method, params) {
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  if (output_assay %in% SummarizedExperiment::assayNames(x)) {
    stop("`output_assay` already exists: ", output_assay, call. = FALSE)
  }
  if (!is.matrix(values) || !identical(dim(values), dim(context$data)) ||
      anyNA(values) || any(!is.finite(values))) {
    stop("Motion correction produced a malformed or non-finite assay",
         call. = FALSE)
  }
  dimnames(values) <- dimnames(context$data)
  out <- x
  SummarizedExperiment::assay(out, output_assay, withDimnames = FALSE) <-
    values
  source_contract <- S4Vectors::metadata(x)$nirs$assays[[context$assay_name]]
  out <- .nirs_set_assay_contract(out, output_assay, list(
    kind = if (is.null(source_contract$kind)) {
      "motion_corrected_signal"
    } else {
      as.character(source_contract$kind)
    },
    unit = source_contract$unit,
    source_assay = context$assay_name,
    method = method,
    motion_corrected = TRUE
  ))
  .nirs_append_step(
    out,
    method,
    params = c(
      list(
        implementation_version = "0.3.0",
        source_fingerprint = context$fingerprint,
        sampling_rate_hz = context$fs,
        sample_count = nrow(context$data),
        measurement_count = ncol(context$data),
        reference = .nirs_motion_reference
      ),
      params
    ),
    input_assay = context$assay_name,
    output_assay = output_assay
  )
}
