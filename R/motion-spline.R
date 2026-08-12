.nirs_csaps_values <- function(x, y, p) {
  n <- length(x)
  if (n != length(y) || n < 1L) {
    stop("Cubic smoothing spline inputs are misaligned", call. = FALSE)
  }
  if (n <= 2L || p == 1) return(as.numeric(y))
  dx <- diff(x)
  if (any(!is.finite(dx)) || any(dx <= 0)) {
    stop("Cubic smoothing spline sites must increase", call. = FALSE)
  }
  m <- n - 2L
  r <- matrix(0, nrow = m, ncol = m)
  diag(r) <- 2 * (dx[seq_len(m)] + dx[seq_len(m) + 1L])
  if (m > 1L) {
    off <- dx[2:m]
    r[cbind(seq_len(m - 1L), 2:m)] <- off
    r[cbind(2:m, seq_len(m - 1L))] <- off
  }
  q <- matrix(0, nrow = m, ncol = n)
  for (i in seq_len(m)) {
    q[i, i] <- 1 / dx[[i]]
    q[i, i + 1L] <- -(1 / dx[[i]] + 1 / dx[[i + 1L]])
    q[i, i + 2L] <- 1 / dx[[i + 1L]]
  }
  pp <- 6 * (1 - p)
  a <- pp * tcrossprod(q) + p * r
  slope <- diff(y) / dx
  u <- tryCatch(
    solve(a, diff(slope)),
    error = function(e) {
      stop("Cubic smoothing spline system is singular", call. = FALSE)
    }
  )
  padded_u <- c(0, as.numeric(u), 0)
  d1 <- diff(padded_u) / dx
  d2 <- diff(c(0, d1, 0))
  fitted <- as.numeric(y) - pp * d2
  if (anyNA(fitted) || any(!is.finite(fitted))) {
    stop("Cubic smoothing spline produced non-finite values",
         call. = FALSE)
  }
  fitted
}

.nirs_spline_window <- function(segment_length, fs) {
  if (segment_length < 1L) return(0L)
  if (segment_length < 0.3 * fs) {
    segment_length
  } else if (segment_length < 3 * fs) {
    max(1L, as.integer(floor(0.3 * fs)))
  } else {
    max(1L, as.integer(floor(segment_length / 10)))
  }
}

.nirs_mean_head <- function(x, count) {
  mean(x[seq_len(min(length(x), count))])
}

.nirs_mean_tail <- function(x, count) {
  n <- length(x)
  mean(x[seq.int(max(1L, n - count + 1L), n)])
}

.nirs_spline_channel <- function(y, time, bad, p, fs) {
  runs <- rle(bad)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  artifact <- which(runs$values)
  if (!length(artifact)) {
    return(list(value = y, intervals = 0L, all_clean = TRUE))
  }
  if (all(bad)) {
    stop("Spline correction requires at least one clean sample",
         call. = FALSE)
  }
  artifact_start <- starts[artifact]
  artifact_end <- ends[artifact]
  out <- y

  for (k in seq_along(artifact_start)) {
    index <- seq.int(artifact_start[[k]], artifact_end[[k]])
    fitted <- .nirs_csaps_values(time[index], y[index], p)
    out[index] <- y[index] - fitted
  }

  first_start <- artifact_start[[1L]]
  first_end <- artifact_end[[1L]]
  first_index <- seq.int(first_start, first_end)
  current_window <- .nirs_spline_window(length(first_index), fs)
  if (first_start > 1L) {
    previous <- seq_len(first_start - 1L)
    target <- .nirs_mean_tail(out[previous], .nirs_spline_window(
      length(previous), fs
    ))
    current <- .nirs_mean_head(out[first_index], current_window)
  } else {
    next_end <- if (length(artifact_start) > 1L) {
      artifact_start[[2L]] - 1L
    } else {
      length(y)
    }
    following <- seq.int(first_end + 1L, next_end)
    target <- .nirs_mean_head(y[following], .nirs_spline_window(
      length(following), fs
    ))
    current <- .nirs_mean_tail(out[first_index], current_window)
  }
  out[first_index] <- out[first_index] - current + target

  if (length(artifact_start) > 1L) {
    for (k in seq_len(length(artifact_start) - 1L)) {
      clean <- seq.int(
        artifact_end[[k]] + 1L,
        artifact_start[[k + 1L]] - 1L
      )
      previous_artifact <- seq.int(
        artifact_start[[k]], artifact_end[[k]]
      )
      target <- .nirs_mean_tail(
        out[previous_artifact],
        .nirs_spline_window(length(previous_artifact), fs)
      )
      current <- .nirs_mean_head(
        y[clean], .nirs_spline_window(length(clean), fs)
      )
      out[clean] <- y[clean] - current + target

      next_artifact <- seq.int(
        artifact_start[[k + 1L]], artifact_end[[k + 1L]]
      )
      target <- .nirs_mean_tail(
        out[clean], .nirs_spline_window(length(clean), fs)
      )
      current <- .nirs_mean_head(
        out[next_artifact],
        .nirs_spline_window(length(next_artifact), fs)
      )
      out[next_artifact] <- out[next_artifact] - current + target
    }
  }

  last_end <- artifact_end[[length(artifact_end)]]
  if (last_end < length(y)) {
    tail_index <- seq.int(last_end + 1L, length(y))
    last_artifact <- seq.int(
      artifact_start[[length(artifact_start)]], last_end
    )
    target <- .nirs_mean_tail(
      out[last_artifact],
      .nirs_spline_window(length(last_artifact), fs)
    )
    current <- .nirs_mean_head(
      y[tail_index], .nirs_spline_window(length(tail_index), fs)
    )
    out[tail_index] <- y[tail_index] - current + target
  }
  if (anyNA(out) || any(!is.finite(out))) {
    stop("Spline correction produced non-finite values", call. = FALSE)
  }
  list(
    value = as.numeric(out),
    intervals = length(artifact_start),
    all_clean = FALSE
  )
}

#' Correct detected fNIRS motion intervals with cubic smoothing splines
#'
#' Uses the MATLAB `csaps` smoothing-parameter definition and stitches
#' corrected artifact and clean segments with short local means.
#'
#' @param x A governed `PhysioExperiment`.
#' @param mask A compatible `nirs_motion_mask`, or `NULL` to detect once.
#' @param assay_name Exact source assay name.
#' @param output_assay Exact new assay name.
#' @param p Cubic smoothing parameter in `[0, 1]`.
#' @param t_motion Positive detection lag window in seconds, used only when
#'   `mask = NULL`.
#' @param t_mask Non-negative detection-mask expansion in seconds, used only
#'   when `mask = NULL`.
#' @param sd_threshold Positive derivative-SD detection multiplier, used only
#'   when `mask = NULL`.
#' @param amplitude_threshold Positive absolute-change detection threshold,
#'   used only when `mask = NULL`.
#' @param group_wavelengths Whether detection unions wavelengths from each
#'   source-detector pair, used only when `mask = NULL`.
#'
#' @return A clone of `x` with one motion-corrected assay.
#' @export
splineMotionCorrect <- function(
    x,
    mask = NULL,
    assay_name = "OD",
    output_assay = paste0(assay_name, "_spline"),
    p = 0.99,
    t_motion = 0.5,
    t_mask = 1,
    sd_threshold = 20,
    amplitude_threshold = 0.5,
    group_wavelengths = TRUE) {
  context <- .nirs_motion_context(x, assay_name, min_samples = 2L)
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  p <- .nirs_motion_scalar(p, "p", lower = 0, upper = 1)
  if (is.null(mask)) {
    mask <- motionArtifactDetect(
      x,
      assay_name = context$assay_name,
      t_motion = t_motion,
      t_mask = t_mask,
      sd_threshold = sd_threshold,
      amplitude_threshold = amplitude_threshold,
      group_wavelengths = group_wavelengths
    )
  } else {
    mask <- .nirs_validate_motion_mask(mask, context)
  }

  corrected <- matrix(
    NA_real_,
    nrow = nrow(context$data),
    ncol = ncol(context$data),
    dimnames = dimnames(context$data)
  )
  diagnostics <- vector("list", ncol(context$data))
  for (j in seq_len(ncol(context$data))) {
    result <- .nirs_spline_channel(
      context$data[, j],
      context$time,
      mask$sample_by_measurement[, j],
      p,
      context$fs
    )
    corrected[, j] <- result$value
    diagnostics[[j]] <- list(
      measurement_id = context$identity$measurement_id[[j]],
      artifact_intervals = result$intervals,
      identity_branch = result$all_clean
    )
  }
  .nirs_add_motion_assay(
    x,
    context,
    output_assay,
    corrected,
    "splineMotionCorrect",
    params = list(
      p = p,
      p_semantics = "MATLAB csaps variational smoothing parameter",
      dt_short_seconds = 0.3,
      dt_long_seconds = 3,
      mask_fingerprint = mask$fingerprint,
      artifact_intervals = nrow(mask$intervals),
      artifact_samples = sum(mask$sample_by_measurement),
      detection_parameters = mask$parameters,
      diagnostics = diagnostics
    )
  )
}
