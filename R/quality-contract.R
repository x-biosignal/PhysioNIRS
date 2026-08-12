.nirs_quality_schema <- "1.0.0"

.nirs_quality_reference <- list(
  pollonini_sci_doi = "10.1117/1.JBO.19.8.086007",
  pollonini_phoebe_doi = "10.1364/BOE.7.005104",
  mne_version = "1.3.1",
  mne_sci_sha256 =
    "403842797280f8e64d181763bfa1a48ba4b0e504108c9c2bb5852425e64e8315",
  mne_nirs_commit = "0a5081735144b902a3953e81d010420e1210c556",
  mne_nirs_segmented_sha256 =
    "6c1b6a01adf01acedff097bf485c8d40116b77b839aaaaa5be66fa7792535e35",
  mne_nirs_peak_power_sha256 =
    "da6ff09d847fb671d6dc922375d17ec23d753568837214c341a5c7a5eec9caf2"
)

.nirs_quality_fields <- c(
  "schema_version", "metric", "assay_name", "identity_kind", "channel_id",
  "source_index", "detector_index", "wavelength_nm",
  "window_start_sample", "window_end_sample", "window_start_time",
  "window_end_time", "score", "pass", "threshold", "parameters",
  "source_fingerprint", "identity_fingerprint", "implementation"
)

.nirs_quality_number <- function(
    value, arg, lower = -Inf, upper = Inf, integer = FALSE,
    lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 1L ||
      is.na(value) || !is.finite(value) ||
      (if (lower_open) value <= lower else value < lower) ||
      (if (upper_open) value >= upper else value > upper) ||
      (integer && value != floor(value))) {
    stop("`", arg, "` must be one finite ",
         if (integer) "integer " else "",
         "value in the governed range", call. = FALSE)
  }
  if (integer) as.integer(value) else as.numeric(value)
}

.nirs_quality_exact_samples <- function(seconds, fs, arg) {
  seconds <- .nirs_quality_number(
    seconds, arg, lower = 0, lower_open = TRUE
  )
  exact <- seconds * fs
  rounded <- round(exact)
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(exact))
  if (!is.finite(exact) || abs(exact - rounded) > tolerance ||
      rounded < 1 || rounded > .Machine$integer.max) {
    stop("`", arg, "` must correspond to an exact positive sample count",
         call. = FALSE)
  }
  as.integer(rounded)
}

.nirs_quality_windows <- function(
    n_time, time, fs, window_seconds, step_seconds) {
  if (is.null(window_seconds)) {
    if (!is.null(step_seconds)) {
      stop("`step_seconds` must be NULL for a complete-record window",
           call. = FALSE)
    }
    start <- 1L
    end <- as.integer(n_time)
  } else {
    window_samples <- .nirs_quality_exact_samples(
      window_seconds, fs, "window_seconds"
    )
    if (is.null(step_seconds)) step_seconds <- window_seconds
    step_samples <- .nirs_quality_exact_samples(
      step_seconds, fs, "step_seconds"
    )
    if (window_samples > n_time) {
      stop("No complete quality window fits in `x`", call. = FALSE)
    }
    start <- seq.int(1L, n_time - window_samples + 1L, by = step_samples)
    end <- start + window_samples - 1L
  }
  ids <- sprintf("w%04d", seq_along(start))
  list(
    id = ids,
    start = as.integer(start),
    end = as.integer(end),
    start_time = as.numeric(time[start]),
    end_time = as.numeric(time[end])
  )
}

.nirs_quality_context <- function(
    x, assay_name, window_seconds, step_seconds, min_samples = 8L) {
  context <- .nirs_short_context(
    x, assay_name, min_samples = min_samples, require_contract = TRUE
  )
  identity_kind <- unique(context$identity$identity_kind)
  if (length(identity_kind) != 1L ||
      !identity_kind %in% c("measurement", "pair")) {
    stop("Governed NIRS identity kind is malformed", call. = FALSE)
  }
  if (anyNA(context$data) || any(!is.finite(context$data))) {
    stop("Quality input assay must contain only finite real values",
         call. = FALSE)
  }
  if (!is.numeric(context$data) || is.complex(context$data)) {
    stop("Quality input assay must be a real numeric matrix",
         call. = FALSE)
  }
  windows <- .nirs_quality_windows(
    nrow(context$data), context$time, context$fs,
    window_seconds, step_seconds
  )
  context$windows <- windows
  context$identity_kind <- identity_kind
  context
}

.nirs_quality_groups <- function(identity) {
  key <- paste(identity$source_index, identity$detector_index, sep = ":")
  split(seq_along(key), factor(key, levels = unique(key)))
}

.nirs_quality_filter <- function(context, l_freq, h_freq, order) {
  l_freq <- .nirs_quality_number(
    l_freq, "l_freq", lower = 0, lower_open = TRUE
  )
  h_freq <- .nirs_quality_number(
    h_freq, "h_freq", lower = l_freq, lower_open = TRUE,
    upper = context$fs / 2, upper_open = TRUE
  )
  order <- .nirs_quality_number(
    order, "order", lower = 1, upper = 20, integer = TRUE
  )
  if (nrow(context$data) < max(8L, 3L * order + 1L)) {
    stop("Quality input is too short for the requested Butterworth order",
         call. = FALSE)
  }
  filter <- signal::butter(
    order, c(l_freq, h_freq) / (context$fs / 2), type = "pass"
  )
  value <- vapply(
    seq_len(ncol(context$data)),
    function(j) .nirs_filtfilt_pad0(filter, context$data[, j]),
    numeric(nrow(context$data))
  )
  value <- as.matrix(value)
  dimnames(value) <- dimnames(context$data)
  if (anyNA(value) || any(!is.finite(value))) {
    stop("Zero-phase quality filter produced non-finite values",
         call. = FALSE)
  }
  list(
    value = value,
    filter = filter,
    l_freq = l_freq,
    h_freq = h_freq,
    order = order
  )
}

.nirs_near_constant <- function(x) {
  scale <- max(1, max(abs(x)))
  stats::sd(x) <= sqrt(.Machine$double.eps) * scale
}

.nirs_quality_fingerprint <- function(result) {
  value <- result
  attr(value, "quality_fingerprint") <- NULL
  .nirs_sha256(value)
}

.nirs_quality_result <- function(
    context, metric, score, threshold, parameters, implementation) {
  windows <- context$windows
  channel_id <- as.character(context$identity$channel_id)
  score <- as.matrix(score)
  if (!identical(dim(score), c(length(channel_id), length(windows$id))) ||
      anyNA(score) || any(!is.finite(score))) {
    stop("Quality metric produced a malformed or non-finite score matrix",
         call. = FALSE)
  }
  dimnames(score) <- list(channel_id, windows$id)
  pass <- score >= threshold
  dimnames(pass) <- dimnames(score)
  result <- list(
    schema_version = .nirs_quality_schema,
    metric = metric,
    assay_name = context$assay_name,
    identity_kind = context$identity_kind,
    channel_id = channel_id,
    source_index = as.integer(context$identity$source_index),
    detector_index = as.integer(context$identity$detector_index),
    wavelength_nm = as.numeric(context$identity$wavelength_nm),
    window_start_sample = windows$start,
    window_end_sample = windows$end,
    window_start_time = windows$start_time,
    window_end_time = windows$end_time,
    score = score,
    pass = pass,
    threshold = as.numeric(threshold),
    parameters = parameters,
    source_fingerprint = context$source_fingerprint,
    identity_fingerprint = context$identity_fingerprint,
    implementation = implementation
  )
  class(result) <- c("nirs_quality", "list")
  attr(result, "quality_fingerprint") <- .nirs_quality_fingerprint(result)
  result
}

.nirs_quality_plain <- function(x) {
  # NULL is a plain, serialization-stable sentinel used for optional
  # parameters (e.g. step_seconds, noise_range_hz, filter). It must count
  # as plain. Handle it explicitly: R >= 4.4.0 made is.atomic(NULL) return
  # FALSE, so relying on the atomic leaf test below would wrongly reject
  # every governed result that carries a NULL-valued parameter.
  if (is.null(x)) {
    return(TRUE)
  }
  if (is.environment(x) || methods::is(x, "externalptr") ||
      inherits(x, "connection") || is.function(x) || isS4(x)) {
    return(FALSE)
  }
  if (is.list(x)) {
    return(all(vapply(x, .nirs_quality_plain, logical(1))))
  }
  is.atomic(x) && !is.object(x)
}

.nirs_validate_quality <- function(quality, x, context = NULL) {
  if (!inherits(quality, "nirs_quality") || !is.list(quality) ||
      !identical(names(quality), .nirs_quality_fields) ||
      !identical(quality$schema_version, .nirs_quality_schema) ||
      !is.character(quality$metric) || length(quality$metric) != 1L ||
      is.na(quality$metric) || !nzchar(quality$metric) ||
      !is.character(quality$assay_name) ||
      length(quality$assay_name) != 1L || is.na(quality$assay_name) ||
      !nzchar(quality$assay_name) ||
      !is.character(quality$identity_kind) ||
      length(quality$identity_kind) != 1L ||
      !quality$identity_kind %in% c("measurement", "pair")) {
    stop("`quality` must be a complete governed `nirs_quality` result",
         call. = FALSE)
  }
  if (is.null(context)) {
    context <- .nirs_quality_context(
      x, quality$assay_name, window_seconds = NULL,
      step_seconds = NULL, min_samples = 2L
    )
  }
  ids <- as.character(context$identity$channel_id)
  n_window <- length(quality$window_start_sample)
  matrix_ok <- is.matrix(quality$score) && is.numeric(quality$score) &&
    !is.complex(quality$score) && !anyNA(quality$score) &&
    all(is.finite(quality$score)) &&
    identical(dim(quality$score), c(length(ids), n_window)) &&
    identical(rownames(quality$score), ids) &&
    is.matrix(quality$pass) && is.logical(quality$pass) &&
    !anyNA(quality$pass) &&
    identical(dim(quality$pass), dim(quality$score)) &&
    identical(dimnames(quality$pass), dimnames(quality$score))
  vectors_ok <- n_window >= 1L &&
    is.integer(quality$window_start_sample) &&
    is.integer(quality$window_end_sample) &&
    length(quality$window_end_sample) == n_window &&
    is.numeric(quality$window_start_time) &&
    length(quality$window_start_time) == n_window &&
    is.numeric(quality$window_end_time) &&
    length(quality$window_end_time) == n_window &&
    all(is.finite(quality$window_start_time)) &&
    all(is.finite(quality$window_end_time)) &&
    all(quality$window_start_sample >= 1L) &&
    all(quality$window_end_sample <= nrow(x)) &&
    all(quality$window_start_sample <= quality$window_end_sample)
  window_ids <- sprintf("w%04d", seq_len(n_window))
  windows_semantic <- vectors_ok &&
    identical(colnames(quality$score), window_ids) &&
    identical(
      quality$window_start_time,
      as.numeric(context$time[quality$window_start_sample])
    ) &&
    identical(
      quality$window_end_time,
      as.numeric(context$time[quality$window_end_sample])
    ) &&
    (n_window == 1L ||
      all(diff(quality$window_start_sample) > 0L))
  identity_ok <- identical(quality$channel_id, ids) &&
    identical(quality$source_index,
              as.integer(context$identity$source_index)) &&
    identical(quality$detector_index,
              as.integer(context$identity$detector_index)) &&
    identical(quality$wavelength_nm,
              as.numeric(context$identity$wavelength_nm)) &&
    identical(quality$identity_kind, context$identity_kind) &&
    identical(quality$source_fingerprint, context$source_fingerprint) &&
    identical(quality$identity_fingerprint,
              context$identity_fingerprint)
  threshold_ok <- is.numeric(quality$threshold) &&
    is.null(dim(quality$threshold)) && length(quality$threshold) == 1L &&
    !is.na(quality$threshold) && is.finite(quality$threshold) &&
    identical(quality$pass, quality$score >= quality$threshold)
  fingerprint <- attr(quality, "quality_fingerprint", exact = TRUE)
  fingerprint_ok <- is.character(fingerprint) &&
    length(fingerprint) == 1L && !is.na(fingerprint) &&
    identical(fingerprint, .nirs_quality_fingerprint(quality))
  if (!matrix_ok || !vectors_ok || !windows_semantic ||
      !identity_ok || !threshold_ok ||
      !.nirs_quality_plain(quality$parameters) ||
      !.nirs_quality_plain(quality$implementation) || !fingerprint_ok) {
    stop("`quality` is stale, malformed, or belongs to another source object",
         call. = FALSE)
  }
  quality
}

#' Display a governed NIRS quality result
#'
#' @param x A `nirs_quality` result.
#' @param ... Unused.
#' @export
print.nirs_quality <- function(x, ...) {
  cat(
    "NIRS quality <", x$metric, ">: ", nrow(x$score), " channels x ",
    ncol(x$score), " windows; pass ",
    sum(x$pass), "/", length(x$pass), "\n", sep = ""
  )
  invisible(x)
}
