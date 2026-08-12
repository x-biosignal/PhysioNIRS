.nirs_filtfilt_pad0 <- function(filter, x) {
  # Match the pinned SciPy/TDDR padlen=0 boundary state. signal::filtfilt()
  # applies a different end-padding convention, so the two passes are explicit.
  b <- as.numeric(filter$b)
  a <- as.numeric(filter$a)
  gain <- sum(b) / sum(a)
  if (!is.finite(gain)) {
    stop("Butterworth filter has non-finite steady-state gain",
         call. = FALSE)
  }
  forward <- signal::filter(
    filter,
    x,
    init.x = rep(x[[1L]], length(b) - 1L),
    init.y = rep(x[[1L]] * gain, length(a) - 1L)
  )
  backward_input <- rev(as.numeric(forward))
  backward <- signal::filter(
    filter,
    backward_input,
    init.x = rep(backward_input[[1L]], length(b) - 1L),
    init.y = rep(backward_input[[1L]] * gain, length(a) - 1L)
  )
  rev(as.numeric(backward))
}

.nirs_tddr_robust <- function(derivative, tune, max_iter) {
  mu <- mean(derivative)
  weight <- rep(1, length(derivative))
  if (!is.finite(mu) || max(abs(derivative - mu)) == 0) {
    return(list(identity = TRUE, derivative = derivative, iterations = 0L))
  }
  iterations <- 0L
  for (iteration in seq_len(max_iter)) {
    dev <- abs(derivative - mu)
    sigma <- 1.4826 * stats::median(dev)
    if (!is.finite(sigma) || sigma <= 0) {
      return(list(
        identity = TRUE,
        derivative = derivative,
        iterations = as.integer(iteration)
      ))
    }
    scaled <- dev / (sigma * tune)
    weight <- ((1 - scaled^2) * (scaled < 1))^2
    weight_sum <- sum(weight)
    if (!is.finite(weight_sum) || weight_sum <= 0) {
      return(list(
        identity = TRUE,
        derivative = derivative,
        iterations = as.integer(iteration)
      ))
    }
    updated <- sum(weight * derivative) / weight_sum
    iterations <- as.integer(iteration)
    tolerance <- sqrt(.Machine$double.eps) *
      max(abs(updated), abs(mu))
    if (abs(updated - mu) < tolerance) {
      mu <- updated
      break
    }
    mu <- updated
  }
  final_dev <- abs(derivative - mu)
  final_sigma <- 1.4826 * stats::median(final_dev)
  if (!is.finite(final_sigma) || final_sigma <= 0) {
    return(list(
      identity = TRUE,
      derivative = derivative,
      iterations = iterations
    ))
  }
  final_scaled <- final_dev / (final_sigma * tune)
  weight <- ((1 - final_scaled^2) * (final_scaled < 1))^2
  if (!is.finite(sum(weight)) || sum(weight) <= 0) {
    return(list(
      identity = TRUE,
      derivative = derivative,
      iterations = iterations
    ))
  }
  list(
    identity = FALSE,
    derivative = weight * (derivative - mu),
    iterations = iterations
  )
}

.nirs_tddr_channel <- function(
    x, fs, cutoff_hz, tune, max_iter) {
  source_mean <- mean(x)
  centered <- x - source_mean
  normalized_cutoff <- 2 * cutoff_hz / fs
  branch <- "butterworth_padlen0"
  if (normalized_cutoff >= 1) {
    low <- centered
    branch <- "all_low_frequency"
  } else {
    filter <- signal::butter(
      n = 3L, W = normalized_cutoff, type = "low"
    )
    low <- .nirs_filtfilt_pad0(filter, centered)
  }
  high <- centered - low
  robust <- .nirs_tddr_robust(diff(low), tune, max_iter)
  if (robust$identity) {
    return(list(
      value = x,
      iterations = robust$iterations,
      identity = TRUE,
      filter_branch = branch
    ))
  }
  repaired_low <- c(0, cumsum(robust$derivative))
  repaired_low <- repaired_low - mean(repaired_low)
  value <- repaired_low + high + source_mean
  value <- value + (source_mean - mean(value))
  list(
    value = as.numeric(value),
    iterations = robust$iterations,
    identity = FALSE,
    filter_branch = branch
  )
}

#' Repair fNIRS motion artifacts with TDDR
#'
#' Implements temporal derivative distribution repair with a third-order
#' Butterworth low-pass and Tukey bisquare derivative weights.
#'
#' @param x A governed `PhysioExperiment`.
#' @param assay_name Exact source assay name.
#' @param output_assay Exact new assay name.
#' @param cutoff_hz Positive low-pass cutoff below Nyquist.
#' @param tune Positive Tukey bisquare tuning constant.
#' @param max_iter Positive integer maximum robust iterations.
#'
#' @return A clone of `x` with one motion-corrected assay.
#' @export
tddr <- function(
    x,
    assay_name = "OD",
    output_assay = paste0(assay_name, "_tddr"),
    cutoff_hz = 0.5,
    tune = 4.685,
    max_iter = 50L) {
  context <- .nirs_motion_context(x, assay_name, min_samples = 5L)
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  cutoff_hz <- .nirs_motion_scalar(
    cutoff_hz, "cutoff_hz", lower = .Machine$double.xmin
  )
  tune <- .nirs_motion_scalar(
    tune, "tune", lower = .Machine$double.xmin
  )
  max_iter <- .nirs_motion_scalar(
    max_iter, "max_iter", lower = 1, upper = .Machine$integer.max,
    integer = TRUE
  )
  nyquist <- context$fs / 2
  cutoff_tolerance <- max(1e-12, 1e-12 * nyquist)
  if (cutoff_hz > nyquist + cutoff_tolerance) {
    stop("`cutoff_hz` must not exceed Nyquist (", nyquist, " Hz)",
         call. = FALSE)
  }
  if (cutoff_hz >= nyquist - cutoff_tolerance) {
    cutoff_hz <- nyquist
  }

  corrected <- matrix(
    NA_real_,
    nrow = nrow(context$data),
    ncol = ncol(context$data),
    dimnames = dimnames(context$data)
  )
  diagnostics <- vector("list", ncol(context$data))
  for (j in seq_len(ncol(context$data))) {
    result <- .nirs_tddr_channel(
      context$data[, j], context$fs, cutoff_hz, tune, max_iter
    )
    corrected[, j] <- result$value
    diagnostics[[j]] <- list(
      measurement_id = context$identity$measurement_id[[j]],
      iterations = result$iterations,
      identity_branch = result$identity,
      filter_branch = result$filter_branch
    )
  }
  .nirs_add_motion_assay(
    x,
    context,
    output_assay,
    corrected,
    "tddr",
    params = list(
      cutoff_hz = cutoff_hz,
      filter_order = 3L,
      filter_type = "Butterworth low-pass, forward/reverse, padlen=0",
      tune = tune,
      max_iter = max_iter,
      diagnostics = diagnostics
    )
  )
}
