.nirs_max_lag_change <- function(x, window) {
  n <- length(x)
  out <- numeric(n - 1L)
  for (lag in seq_len(window)) {
    current <- abs(x[(lag + 1L):n] - x[seq_len(n - lag)])
    if (lag > 1L) current <- c(current, numeric(lag - 1L))
    out <- pmax(out, current)
  }
  out
}

.nirs_expand_diff_mask <- function(bad_diff, buffer, n_sample) {
  index <- which(bad_diff)
  if (!length(index)) return(rep(FALSE, n_sample))
  expanded <- unique(as.vector(outer(
    index, seq.int(-buffer, buffer), `+`
  )))
  expanded <- expanded[expanded >= 1L & expanded <= n_sample - 1L]
  out <- rep(FALSE, n_sample)
  out[expanded + 1L] <- TRUE
  out
}

#' Detect fNIRS motion artifacts by channel
#'
#' Uses strict amplitude and derivative-SD thresholds over a finite lag window.
#' Detections are expanded in time and can be unioned across wavelengths from
#' the same source-detector pair.
#'
#' @param x A governed `PhysioExperiment`.
#' @param assay_name Exact source assay name.
#' @param t_motion Positive lag window in seconds.
#' @param t_mask Non-negative expansion on each side in seconds.
#' @param sd_threshold Positive multiple of channel derivative SD.
#' @param amplitude_threshold Positive absolute-change threshold.
#' @param group_wavelengths Whether to union measurements sharing an exact
#'   source-detector pair.
#'
#' @return A typed `nirs_motion_mask` aligned to assay rows and columns.
#' @export
motionArtifactDetect <- function(
    x,
    assay_name = "OD",
    t_motion = 0.5,
    t_mask = 1,
    sd_threshold = 20,
    amplitude_threshold = 0.5,
    group_wavelengths = TRUE) {
  context <- .nirs_motion_context(x, assay_name, min_samples = 2L)
  t_motion <- .nirs_motion_scalar(
    t_motion, "t_motion", lower = .Machine$double.xmin
  )
  t_mask <- .nirs_motion_scalar(t_mask, "t_mask", lower = 0)
  sd_threshold <- .nirs_motion_scalar(
    sd_threshold, "sd_threshold", lower = .Machine$double.xmin
  )
  amplitude_threshold <- .nirs_motion_scalar(
    amplitude_threshold, "amplitude_threshold",
    lower = .Machine$double.xmin
  )
  group_wavelengths <- .nirs_motion_flag(
    group_wavelengths, "group_wavelengths"
  )

  window <- max(1L, as.integer(round(t_motion * context$fs)))
  if (window >= nrow(context$data)) {
    stop(
      "`t_motion` requires ", window,
      " lag samples but the assay has only ", nrow(context$data),
      " rows",
      call. = FALSE
    )
  }
  buffer <- as.integer(round(t_mask * context$fs))
  raw_mask <- matrix(
    FALSE,
    nrow = nrow(context$data),
    ncol = ncol(context$data),
    dimnames = dimnames(context$data)
  )
  for (j in seq_len(ncol(context$data))) {
    derivative <- diff(context$data[, j])
    derivative_sd <- stats::sd(derivative)
    if (!is.finite(derivative_sd)) {
      stop("Cannot calculate derivative SD for channel ", j,
           call. = FALSE)
    }
    change <- .nirs_max_lag_change(context$data[, j], window)
    bad <- change > sd_threshold * derivative_sd |
      change > amplitude_threshold
    raw_mask[, j] <- .nirs_expand_diff_mask(
      bad, buffer, nrow(context$data)
    )
  }

  result_mask <- raw_mask
  if (group_wavelengths) {
    key <- paste(
      context$identity$source_index,
      context$identity$detector_index,
      sep = ":"
    )
    groups <- split(seq_along(key), factor(key, levels = unique(key)))
    for (index in groups) {
      pair_mask <- apply(
        raw_mask[, index, drop = FALSE], 1L, any
      )
      result_mask[, index] <- pair_mask
    }
  }
  global <- apply(result_mask, 1L, any)
  parameters <- list(
    t_motion = t_motion,
    t_mask = t_mask,
    sd_threshold = sd_threshold,
    amplitude_threshold = amplitude_threshold,
    group_wavelengths = group_wavelengths,
    lag_window_samples = window,
    mask_buffer_samples = buffer,
    comparison = "strict_greater_than",
    diff_alignment = "sample_after_change"
  )
  out <- list(
    sample_by_measurement = result_mask,
    global = global,
    intervals = .nirs_motion_intervals(
      result_mask, context$time, context$identity
    ),
    parameters = parameters,
    assay_name = context$assay_name,
    sampling_rate_hz = context$fs,
    measurement_id = as.character(context$identity$measurement_id),
    input_fingerprint = context$fingerprint
  )
  class(out) <- "nirs_motion_mask"
  out$fingerprint <- .nirs_motion_mask_fingerprint(out)
  out
}
