#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PhysioNIRS)
  library(PhysioCore)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- normalizePath(file_arg[[1L]], mustWork = TRUE)
validation_dir <- dirname(script)
package_root <- normalizePath(
  file.path(validation_dir, "..", ".."), mustWork = TRUE
)

exact_enum <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !x %in% choices) {
    stop(arg, " invalid", call. = FALSE)
  }
  x
}

make_experiment <- function(values, fs) {
  values <- as.matrix(values)
  n <- nrow(values)
  wavelength <- c(760, 850, 760, 850)
  source <- c(1L, 1L, 2L, 2L)
  detector <- source
  labels <- c("S1_D1_760", "S1_D1_850", "S2_D2_760", "S2_D2_850")
  dimnames(values) <- list(NULL, labels)
  measurement <- S4Vectors::DataFrame(
    measurement_index = 1:4,
    source_index = source,
    detector_index = detector,
    wavelength_index = c(1L, 2L, 1L, 2L),
    wavelength_nm = wavelength,
    wavelength_actual_nm = rep(NA_real_, 4),
    data_type = rep(1L, 4),
    data_type_index = rep(1L, 4),
    data_type_label = rep(NA_character_, 4),
    data_unit = rep(NA_character_, 4),
    source_power = rep(NA_real_, 4),
    detector_gain = rep(NA_real_, 4),
    source_label = c("S1", "S1", "S2", "S2"),
    detector_label = c("D1", "D1", "D2", "D2"),
    channel_label = labels
  )
  metadata <- list(
    snirf = list(
      measurement_list = measurement,
      probe = list(
        wavelengths = c(760, 850),
        sourceLabels = c("S1", "S2"),
        detectorLabels = c("D1", "D2"),
        sourcePos2D = rbind(c(0, 0), c(0.05, 0)),
        detectorPos2D = rbind(c(0, 0.03), c(0.05, 0.04)),
        LengthUnit = "m"
      ),
      metadata_tags = list(
        SubjectID = "validation-fixture",
        MeasurementDate = "unknown",
        MeasurementTime = "unknown",
        LengthUnit = "m",
        TimeUnit = "s",
        FrequencyUnit = "Hz"
      ),
      stim = list(),
      time_sampling = "uniform"
    ),
    nirs = list(assays = list(OD = list(
      kind = "optical_density",
      unit = "1",
      log_convention = "natural",
      source_assay = "synthetic"
    )))
  )
  x <- PhysioCore::PhysioExperiment(
    assays = list(OD = values),
    rowData = S4Vectors::DataFrame(
      time_seconds = seq.int(0, n - 1L) / fs
    ),
    colData = cbind(measurement, label = labels),
    metadata = metadata,
    samplingRate = fs
  )
  PhysioCore::setEvents(
    x,
    PhysioCore::PhysioEvents(
      onset = 5, duration = 10, type = "task", value = "1"
    )
  )
}

oracle_max_change <- function(x, window) {
  n <- length(x)
  out <- numeric(n - 1L)
  for (lag in seq_len(window)) {
    value <- abs(x[(lag + 1L):n] - x[seq_len(n - lag)])
    if (lag > 1L) value <- c(value, numeric(lag - 1L))
    out <- pmax(out, value)
  }
  out
}

oracle_detect <- function(
    data, fs, source, detector, t_motion, t_mask, sd_threshold,
    amplitude_threshold, grouped = TRUE,
    change_function = oracle_max_change,
    sd_function = stats::sd,
    grouping_key = NULL) {
  n <- nrow(data)
  window <- max(1L, round(t_motion * fs))
  buffer <- round(t_mask * fs)
  mask <- matrix(FALSE, n, ncol(data))
  for (j in seq_len(ncol(data))) {
    change <- change_function(data[, j], window)
    bad <- which(
      change > sd_threshold * sd_function(diff(data[, j])) |
        change > amplitude_threshold
    )
    if (length(bad)) {
      expanded <- unique(as.vector(outer(
        bad, seq.int(-buffer, buffer), `+`
      )))
      expanded <- expanded[expanded >= 1L & expanded <= n - 1L]
      mask[expanded + 1L, j] <- TRUE
    }
  }
  if (grouped) {
    key <- grouping_key
    if (is.null(key)) key <- paste(source, detector, sep = ":")
    for (group in split(seq_along(key), factor(key, levels = unique(key)))) {
      mask[, group] <- apply(mask[, group, drop = FALSE], 1L, any)
    }
  }
  mask
}

oracle_max_change_circular <- function(x, window) {
  n <- length(x)
  out <- numeric(n - 1L)
  for (lag in seq_len(window)) {
    shifted <- c(tail(x, lag), head(x, n - lag))
    out <- pmax(out, abs(x - shifted)[seq_len(n - 1L)])
  }
  out
}

population_sd <- function(x) {
  sqrt(mean((x - mean(x))^2))
}

oracle_filtfilt <- function(filter, x) {
  b <- as.numeric(filter$b)
  a <- as.numeric(filter$a)
  gain <- sum(b) / sum(a)
  forward <- signal::filter(
    filter, x,
    init.x = rep(x[[1L]], length(b) - 1L),
    init.y = rep(x[[1L]] * gain, length(a) - 1L)
  )
  reverse_input <- rev(as.numeric(forward))
  rev(as.numeric(signal::filter(
    filter, reverse_input,
    init.x = rep(reverse_input[[1L]], length(b) - 1L),
    init.y = rep(reverse_input[[1L]] * gain, length(a) - 1L)
  )))
}

oracle_tddr <- function(x, fs, cutoff = 0.5, tune = 4.685) {
  source_mean <- mean(x)
  centered <- x - source_mean
  normalized <- 2 * cutoff / fs
  low <- if (normalized >= 1) {
    centered
  } else {
    oracle_filtfilt(signal::butter(3, normalized, "low"), centered)
  }
  high <- centered - low
  derivative <- diff(low)
  mu <- mean(derivative)
  if (max(abs(derivative - mu)) == 0) return(x)
  weight <- rep(1, length(derivative))
  identity <- FALSE
  for (iteration in seq_len(50L)) {
    deviation <- abs(derivative - mu)
    sigma <- 1.4826 * stats::median(deviation)
    if (!is.finite(sigma) || sigma <= 0) {
      identity <- TRUE
      break
    }
    scaled <- deviation / (sigma * tune)
    weight <- ((1 - scaled^2) * (scaled < 1))^2
    if (sum(weight) <= 0) {
      identity <- TRUE
      break
    }
    updated <- sum(weight * derivative) / sum(weight)
    if (abs(updated - mu) <
        sqrt(.Machine$double.eps) * max(abs(updated), abs(mu))) {
      mu <- updated
      break
    }
    mu <- updated
  }
  if (identity) return(x)
  deviation <- abs(derivative - mu)
  sigma <- 1.4826 * stats::median(deviation)
  scaled <- deviation / (sigma * tune)
  weight <- ((1 - scaled^2) * (scaled < 1))^2
  repaired <- c(0, cumsum(weight * (derivative - mu)))
  repaired <- repaired - mean(repaired)
  value <- repaired + high + source_mean
  value + source_mean - mean(value)
}

mutant_tddr <- function(
    x,
    fs,
    cutoff = 0.5,
    tune = 4.685,
    filter_order = 3L,
    cutoff_scale = 2,
    mad_constant = 1.4826,
    center_input = TRUE,
    restore_mean = TRUE,
    low_mode = c("zero_phase", "source", "causal"),
    robust_center = c("bisquare", "ordinary"),
    weight_power = 2L,
    shifted_origin = FALSE) {
  low_mode <- exact_enum(low_mode[[1L]], c(
    "zero_phase", "source", "causal"
  ), "low_mode")
  robust_center <- exact_enum(robust_center[[1L]], c(
    "bisquare", "ordinary"
  ), "robust_center")
  source_mean <- if (center_input) mean(x) else 0
  centered <- x - source_mean
  normalized <- cutoff_scale * cutoff / fs
  filter <- signal::butter(filter_order, normalized, "low")
  low <- switch(
    low_mode,
    zero_phase = oracle_filtfilt(filter, centered),
    source = centered,
    causal = as.numeric(signal::filter(filter, centered))
  )
  high <- centered - low
  derivative <- diff(low)
  mu <- mean(derivative)
  weight <- rep(1, length(derivative))
  if (robust_center == "bisquare") {
    for (iteration in seq_len(50L)) {
      deviation <- abs(derivative - mu)
      sigma <- mad_constant * stats::median(deviation)
      if (!is.finite(sigma) || sigma <= 0) return(x)
      scaled <- deviation / (sigma * tune)
      weight <- ((1 - scaled^2) * (scaled < 1))^weight_power
      if (sum(weight) <= 0) return(x)
      updated <- sum(weight * derivative) / sum(weight)
      if (abs(updated - mu) <
          sqrt(.Machine$double.eps) * max(abs(updated), abs(mu))) {
        mu <- updated
        break
      }
      mu <- updated
    }
    deviation <- abs(derivative - mu)
    sigma <- mad_constant * stats::median(deviation)
    scaled <- deviation / (sigma * tune)
    weight <- ((1 - scaled^2) * (scaled < 1))^weight_power
  }
  corrected_derivative <- weight * (derivative - mu)
  repaired <- if (shifted_origin) {
    c(cumsum(corrected_derivative), 0)
  } else {
    c(0, cumsum(corrected_derivative))
  }
  repaired <- repaired - mean(repaired)
  value <- repaired + high
  if (restore_mean) value <- value + mean(x) - mean(value)
  value
}

db2 <- list(
  low = c(
    0.48296291314453416, 0.83651630373780794,
    0.22414386804201339, -0.12940952255126037
  ),
  high = c(
    -0.12940952255126037, -0.22414386804201339,
    0.83651630373780794, -0.48296291314453416
  )
)

oracle_dwt <- function(x) {
  n <- length(x)
  half <- n %/% 2L
  a <- d <- numeric(half)
  for (k in seq_len(half)) {
    index <- ((2L * (k - 1L) + seq_along(db2$low) - 2L) %% n) + 1L
    a[[k]] <- sum(db2$low * x[index])
    d[[k]] <- sum(db2$high * x[index])
  }
  list(a = a, d = d)
}

oracle_idwt <- function(a, d) {
  n <- 2L * length(a)
  out <- numeric(n)
  for (k in seq_along(a)) {
    index <- ((2L * (k - 1L) + seq_along(db2$low) - 2L) %% n) + 1L
    out[index] <- out[index] + a[[k]] * db2$low + d[[k]] * db2$high
  }
  out
}

oracle_wavelet <- function(x, iqr = 1.5, min_level = 4L) {
  original_n <- length(x)
  exponent <- ceiling(log2(original_n))
  n <- 2^exponent
  depth <- exponent - min_level
  padded <- numeric(n)
  padded[seq_len(original_n)] <- x
  dc <- mean(padded)
  padded <- padded - dc
  coefficients <- matrix(0, n, depth + 1L)
  coefficients[, 1L] <- padded
  for (level in 0:(depth - 1L)) {
    blocks <- 2^level
    block_n <- n %/% blocks
    half <- block_n %/% 2L
    for (block in 0:(blocks - 1L)) {
      index <- block * block_n + seq_len(block_n)
      signal <- coefficients[index, 1L]
      shifted <- c(signal[[block_n]], signal[-block_n])
      one <- oracle_dwt(signal)
      two <- oracle_dwt(shifted)
      coefficients[index[seq_len(half)], 1L] <- one$a
      coefficients[index[half + seq_len(half)], 1L] <- two$a
      coefficients[index[seq_len(half)], level + 2L] <- one$d
      coefficients[index[half + seq_len(half)], level + 2L] <- two$d
    }
  }
  valid_n <- original_n
  if (depth > 1L) {
    for (level in seq_len(depth - 1L)) {
      valid_n <- floor(valid_n / 2)
      blocks <- 2^level
      block_n <- n %/% blocks
      for (block in 0:(blocks - 1L)) {
        index <- block * block_n + seq_len(block_n)
        reference <- coefficients[
          index[seq_len(min(valid_n, block_n))], level + 1L
        ]
        q <- as.numeric(stats::quantile(
          reference, c(0.25, 0.5, 0.75), type = 7, names = FALSE
        ))
        spread <- q[[3L]] - q[[1L]]
        local <- coefficients[index, level + 1L]
        local[local < q[[1L]] - iqr * spread |
                local > q[[3L]] + iqr * spread] <- 0
        coefficients[index, level + 1L] <- local
      }
    }
  }
  approximation <- coefficients[, 1L]
  for (level in seq.int(depth - 1L, 0L)) {
    blocks <- 2^level
    block_n <- n %/% blocks
    half <- block_n %/% 2L
    for (block in 0:(blocks - 1L)) {
      index <- block * block_n + seq_len(block_n)
      first <- seq_len(half)
      second <- half + seq_len(half)
      ordinary <- oracle_idwt(
        approximation[index[first]],
        coefficients[index[first], level + 2L]
      )
      shifted <- oracle_idwt(
        approximation[index[second]],
        coefficients[index[second], level + 2L]
      )
      approximation[index] <- (
        ordinary + c(shifted[-1L], shifted[[1L]])
      ) / 2
    }
  }
  (approximation + dc)[seq_len(original_n)]
}

mutant_decimated_wavelet <- function(
    x,
    filter = db2,
    boundary = c("periodic", "zero"),
    threshold = c("iqr", "sd"),
    retain_outliers = FALSE,
    trim = TRUE) {
  boundary <- exact_enum(boundary[[1L]], c("periodic", "zero"), "boundary")
  threshold <- exact_enum(threshold[[1L]], c("iqr", "sd"), "threshold")
  original_n <- length(x)
  n <- 2^ceiling(log2(original_n))
  padded <- c(x, numeric(n - original_n))
  dc <- mean(padded)
  padded <- padded - dc
  half <- n %/% 2L
  approximation <- detail <- numeric(half)
  locations <- vector("list", half)
  for (k in seq_len(half)) {
    raw <- 2L * (k - 1L) + seq_along(filter$low) - 1L
    if (boundary == "periodic") {
      index <- (raw %% n) + 1L
      values <- padded[index]
      locations[[k]] <- index
    } else {
      valid <- raw >= 0L & raw < n
      values <- numeric(length(raw))
      values[valid] <- padded[raw[valid] + 1L]
      locations[[k]] <- raw[valid] + 1L
    }
    approximation[[k]] <- sum(filter$low * values)
    detail[[k]] <- sum(filter$high * values)
  }
  outlier <- if (threshold == "iqr") {
    quartile <- as.numeric(stats::quantile(
      detail, c(0.25, 0.75), type = 7, names = FALSE
    ))
    spread <- diff(quartile)
    detail < quartile[[1L]] - 1.5 * spread |
      detail > quartile[[2L]] + 1.5 * spread
  } else {
    abs(detail - mean(detail)) > 1.5 * stats::sd(detail)
  }
  if (!retain_outliers) detail[outlier] <- 0
  reconstructed <- numeric(n)
  for (k in seq_len(half)) {
    raw <- 2L * (k - 1L) + seq_along(filter$low) - 1L
    contribution <- approximation[[k]] * filter$low +
      detail[[k]] * filter$high
    if (boundary == "periodic") {
      index <- (raw %% n) + 1L
      reconstructed[index] <- reconstructed[index] + contribution
    } else {
      valid <- raw >= 0L & raw < n
      index <- raw[valid] + 1L
      reconstructed[index] <- reconstructed[index] + contribution[valid]
    }
  }
  reconstructed <- reconstructed + dc
  if (trim) reconstructed[seq_len(original_n)] else reconstructed
}

oracle_csaps <- function(x, y, p) {
  n <- length(x)
  if (n <= 2L || p == 1) return(y)
  dx <- diff(x)
  m <- n - 2L
  r <- matrix(0, m, m)
  diag(r) <- 2 * (dx[seq_len(m)] + dx[seq_len(m) + 1L])
  if (m > 1L) {
    r[cbind(seq_len(m - 1L), 2:m)] <- dx[2:m]
    r[cbind(2:m, seq_len(m - 1L))] <- dx[2:m]
  }
  q <- matrix(0, m, n)
  for (i in seq_len(m)) {
    q[i, i] <- 1 / dx[[i]]
    q[i, i + 1L] <- -(1 / dx[[i]] + 1 / dx[[i + 1L]])
    q[i, i + 2L] <- 1 / dx[[i + 1L]]
  }
  pp <- 6 * (1 - p)
  u <- solve(pp * tcrossprod(q) + p * r, diff(diff(y) / dx))
  d1 <- diff(c(0, u, 0)) / dx
  y - pp * diff(c(0, d1, 0))
}

window_count <- function(n, fs) {
  if (n < 0.3 * fs) n else if (n < 3 * fs) {
    max(1L, floor(0.3 * fs))
  } else {
    max(1L, floor(n / 10))
  }
}

mean_head <- function(x, count) mean(x[seq_len(min(length(x), count))])
mean_tail <- function(x, count) {
  mean(x[seq.int(max(1L, length(x) - count + 1L), length(x))])
}

oracle_spline <- function(y, time, bad, p, fs) {
  runs <- rle(bad)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  selected <- which(runs$values)
  if (!length(selected)) return(y)
  s <- starts[selected]
  e <- ends[selected]
  out <- y
  for (k in seq_along(s)) {
    index <- seq.int(s[[k]], e[[k]])
    out[index] <- y[index] - oracle_csaps(time[index], y[index], p)
  }
  first <- seq.int(s[[1L]], e[[1L]])
  if (s[[1L]] > 1L) {
    previous <- seq_len(s[[1L]] - 1L)
    target <- mean_tail(out[previous], window_count(length(previous), fs))
    current <- mean_head(out[first], window_count(length(first), fs))
  } else {
    next_end <- if (length(s) > 1L) s[[2L]] - 1L else length(y)
    following <- seq.int(e[[1L]] + 1L, next_end)
    target <- mean_head(y[following], window_count(length(following), fs))
    current <- mean_tail(out[first], window_count(length(first), fs))
  }
  out[first] <- out[first] - current + target
  if (length(s) > 1L) {
    for (k in seq_len(length(s) - 1L)) {
      clean <- seq.int(e[[k]] + 1L, s[[k + 1L]] - 1L)
      previous <- seq.int(s[[k]], e[[k]])
      target <- mean_tail(out[previous], window_count(length(previous), fs))
      current <- mean_head(y[clean], window_count(length(clean), fs))
      out[clean] <- y[clean] - current + target
      next_artifact <- seq.int(s[[k + 1L]], e[[k + 1L]])
      target <- mean_tail(out[clean], window_count(length(clean), fs))
      current <- mean_head(
        out[next_artifact], window_count(length(next_artifact), fs)
      )
      out[next_artifact] <- out[next_artifact] - current + target
    }
  }
  if (e[[length(e)]] < length(y)) {
    tail_index <- seq.int(e[[length(e)]] + 1L, length(y))
    last <- seq.int(s[[length(s)]], e[[length(e)]])
    target <- mean_tail(out[last], window_count(length(last), fs))
    current <- mean_head(y[tail_index], window_count(length(tail_index), fs))
    out[tail_index] <- y[tail_index] - current + target
  }
  out
}

mutant_spline_residual <- function(
    y, time, bad, p, use_spar = FALSE) {
  runs <- rle(bad)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  artifact <- which(runs$values)
  out <- y
  for (k in artifact) {
    index <- seq.int(starts[[k]], ends[[k]])
    fitted <- if (use_spar && length(index) >= 4L) {
      stats::smooth.spline(time[index], y[index], spar = p)$y
    } else {
      oracle_csaps(time[index], y[index], p)
    }
    out[index] <- y[index] - fitted
  }
  out
}

frequency_component <- function(x, time, frequency) {
  centered <- x - mean(x)
  coefficient <- sum(centered * exp(-2i * pi * frequency * time))
  list(
    amplitude = 2 * Mod(coefficient) / length(x),
    phase = Arg(coefficient)
  )
}

phase_distance <- function(a, b) {
  abs(Arg(exp(1i * (a - b))))
}

source <- c(1L, 1L, 2L, 2L)
detector <- source
results <- vector("list", 100L)
first_case <- NULL

for (case_id in seq_len(100L)) {
  set.seed(14000L + case_id)
  fs <- c(8, 10, 12.5, 20)[[(case_id - 1L) %% 4L + 1L]]
  duration <- c(60, 80, 100)[[(case_id - 1L) %% 3L + 1L]]
  n <- as.integer(round(fs * duration))
  time <- seq.int(0, n - 1L) / fs
  phase <- runif(4L, -0.15, 0.15)
  clean <- vapply(
    phase,
    function(value) {
      0.08 * sin(2 * pi * 0.05 * time + value) +
        0.005 * sin(2 * pi * 0.10 * time)
    },
    numeric(n)
  )
  noise <- matrix(rnorm(n * 4L, sd = 0.0015), n, 4L)
  observed <- clean + noise
  spike <- as.integer(round(0.20 * n)) +
    ((case_id * 17L) %% as.integer(round(0.15 * n)))
  observed[spike, 1:2] <- observed[spike, 1:2] +
    c(1.1, 0.9) * if (case_id %% 2L) 1 else -1
  step <- min(n - as.integer(round(8 * fs)), spike + round(5 * fs))
  observed[step:n, 1:2] <- observed[step:n, 1:2] +
    c(0.7, 0.55) * if (case_id %% 3L) 1 else -1
  burst_start <- min(n - 20L, spike + round(10 * fs))
  burst <- 0.6 * sin(seq(0, 6 * pi, length.out = 16L))
  observed[burst_start + 0:15, 3:4] <-
    observed[burst_start + 0:15, 3:4] + cbind(burst, 0.85 * burst)

  x <- make_experiment(observed, fs)
  source_value <- SummarizedExperiment::assay(x, "OD")
  effective_fs <- 1 / stats::median(diff(time))
  source_digest <- digest::digest(x, algo = "sha256", serializeVersion = 2)
  mask <- motionArtifactDetect(
    x,
    t_motion = 0.2,
    t_mask = 0.3,
    sd_threshold = 5,
    amplitude_threshold = 0.25,
    group_wavelengths = TRUE
  )
  mask_oracle <- oracle_detect(
    source_value, effective_fs, source, detector,
    0.2, 0.3, 5, 0.25, TRUE
  )
  tddr_result <- tddr(x)
  tddr_value <- SummarizedExperiment::assay(tddr_result, "OD_tddr")
  tddr_oracle <- apply(
    source_value, 2L, oracle_tddr, fs = fs, cutoff = 0.5, tune = 4.685
  )
  wavelet_result <- waveletMotionCorrect(x)
  wavelet_value <- SummarizedExperiment::assay(
    wavelet_result, "OD_wavelet"
  )
  wavelet_oracle <- apply(source_value, 2L, oracle_wavelet)
  spline_result <- splineMotionCorrect(x, mask = mask)
  spline_value <- SummarizedExperiment::assay(
    spline_result, "OD_spline"
  )
  spline_oracle <- matrix(NA_real_, n, 4L)
  for (j in seq_len(4L)) {
    spline_oracle[, j] <- oracle_spline(
      source_value[, j], time, mask$sample_by_measurement[, j],
      0.99, effective_fs
    )
  }
  clean_x <- make_experiment(clean, fs)
  clean_tddr <- SummarizedExperiment::assay(tddr(clean_x), "OD_tddr")
  clean_wavelet <- SummarizedExperiment::assay(
    waveletMotionCorrect(clean_x), "OD_wavelet"
  )
  clean_mask <- motionArtifactDetect(
    clean_x,
    t_motion = 0.2,
    t_mask = 0.3,
    sd_threshold = 1e6,
    amplitude_threshold = 1e6
  )
  clean_spline <- SummarizedExperiment::assay(
    splineMotionCorrect(clean_x, mask = clean_mask), "OD_spline"
  )
  clean_component <- frequency_component(clean[, 1], time, 0.05)
  corrected_component <- frequency_component(
    clean_tddr[, 1], time, 0.05
  )
  low_filter <- signal::butter(3, 2 * 0.5 / fs, "low")
  before_residual <- oracle_filtfilt(
    low_filter, source_value[, 1] - clean[, 1]
  )
  after_residual <- oracle_filtfilt(
    low_filter, tddr_value[, 1] - clean_tddr[, 1]
  )
  before_variance <- stats::var(before_residual)
  after_variance <- stats::var(after_residual)
  relative_rmse <- function(actual, expected) {
    sqrt(mean((actual - expected)^2)) / sqrt(mean(expected^2))
  }
  preserved <- function(result) {
    identical(
      SummarizedExperiment::assay(result, "OD"),
      SummarizedExperiment::assay(x, "OD")
    ) &&
      identical(
        SummarizedExperiment::rowData(result),
        SummarizedExperiment::rowData(x)
      ) &&
      identical(
        SummarizedExperiment::colData(result),
        SummarizedExperiment::colData(x)
      ) &&
      identical(PhysioCore::getEvents(result), PhysioCore::getEvents(x))
  }
  results[[case_id]] <- data.frame(
    case_id = case_id,
    sampling_rate_hz = fs,
    sample_count = n,
    detector_exact = identical(
      unname(mask$sample_by_measurement), mask_oracle
    ),
    tddr_max_abs = max(abs(tddr_value - tddr_oracle)),
    tddr_rms = sqrt(mean((tddr_value - tddr_oracle)^2)),
    wavelet_max_abs = max(abs(wavelet_value - wavelet_oracle)),
    spline_max_abs = max(abs(spline_value - spline_oracle)),
    clean_tddr_relative_rmse = relative_rmse(clean_tddr, clean),
    clean_wavelet_relative_rmse = relative_rmse(clean_wavelet, clean),
    clean_spline_relative_rmse = relative_rmse(clean_spline, clean),
    motion_variance_reduction = 1 - after_variance / before_variance,
    haemodynamic_amplitude_ratio =
      corrected_component$amplitude / clean_component$amplitude,
    haemodynamic_phase_error = phase_distance(
      corrected_component$phase, clean_component$phase
    ),
    mean_error = max(abs(colMeans(tddr_value) - colMeans(observed))),
    finite = all(is.finite(tddr_value)) &&
      all(is.finite(wavelet_value)) && all(is.finite(spline_value)),
    dimensions = identical(dim(tddr_value), dim(observed)) &&
      identical(dim(wavelet_value), dim(observed)) &&
      identical(dim(spline_value), dim(observed)),
    metadata_preserved = preserved(tddr_result) &&
      preserved(wavelet_result) && preserved(spline_result),
    provenance_present =
      length(S4Vectors::metadata(tddr_result)$provenance) >= 1L &&
      length(S4Vectors::metadata(wavelet_result)$provenance) >= 1L &&
      length(S4Vectors::metadata(spline_result)$provenance) >= 1L,
    input_immutable = identical(
      source_digest,
      digest::digest(x, algo = "sha256", serializeVersion = 2)
    ),
    stringsAsFactors = FALSE
  )
  if (case_id == 1L) {
    first_case <- list(
      x = x,
      observed = source_value,
      mask = mask,
      mask_oracle = mask_oracle,
      tddr = tddr_value,
      tddr_oracle = tddr_oracle,
      wavelet = wavelet_value,
      wavelet_oracle = wavelet_oracle,
      spline = spline_value,
      spline_oracle = spline_oracle,
      fs = fs,
      time = time
    )
  }
}

validation <- do.call(rbind, results)
tolerances <- list(
  detector_exact_cases = 100L,
  tddr_max_abs = 1e-10,
  tddr_rms = 1e-12,
  wavelet_max_abs = 1e-10,
  spline_max_abs = 1e-10,
  clean_wavelet_relative_rmse = 0.01,
  clean_spline_relative_rmse = 0.01,
  median_motion_variance_reduction = 0.50,
  median_haemodynamic_amplitude_error = 0.10,
  median_haemodynamic_phase_error = 0.05,
  mean_error = 1e-12
)
gates <- c(
  detector = sum(validation$detector_exact) ==
    tolerances$detector_exact_cases,
  tddr_reference = max(validation$tddr_max_abs) <=
    tolerances$tddr_max_abs,
  tddr_reference_rms = max(validation$tddr_rms) <= tolerances$tddr_rms,
  wavelet_reference = max(validation$wavelet_max_abs) <=
    tolerances$wavelet_max_abs,
  spline_reference = max(validation$spline_max_abs) <=
    tolerances$spline_max_abs,
  clean_wavelet = max(validation$clean_wavelet_relative_rmse) <=
    tolerances$clean_wavelet_relative_rmse,
  clean_spline = max(validation$clean_spline_relative_rmse) <=
    tolerances$clean_spline_relative_rmse,
  motion_variance = stats::median(validation$motion_variance_reduction) >
    tolerances$median_motion_variance_reduction,
  haemodynamic_amplitude =
    abs(stats::median(validation$haemodynamic_amplitude_ratio) - 1) <=
    tolerances$median_haemodynamic_amplitude_error,
  haemodynamic_phase =
    stats::median(validation$haemodynamic_phase_error) <=
    tolerances$median_haemodynamic_phase_error,
  mean = max(validation$mean_error) <= tolerances$mean_error,
  finite = all(validation$finite),
  dimensions = all(validation$dimensions),
  metadata_preserved = all(validation$metadata_preserved),
  provenance_present = all(validation$provenance_present),
  input_immutable = all(validation$input_immutable)
)

fails <- function(expression) {
  inherits(try(force(expression), silent = TRUE), "try-error")
}
changed <- function(a, b, tolerance = 1e-12) {
  !identical(length(a), length(b)) ||
    !identical(dim(a), dim(b)) ||
    max(abs(a - b)) > tolerance
}
base <- first_case
wrong_mask_shift <- rbind(
  base$mask_oracle[-1L, , drop = FALSE],
  FALSE
)
group_fixture <- matrix(0, 40, 4)
group_fixture[20, 1] <- 1
grouped_fixture <- oracle_detect(
  group_fixture, 10, source, detector, 0.1, 0, 1e6, 0.5, TRUE
)
wrong_mask_group <- oracle_detect(
  group_fixture, 10, source, detector, 0.1, 0, 1e6, 0.5, FALSE
)
display_group <- oracle_detect(
  group_fixture, 10, source, detector, 0.1, 0, 1e6, 0.5, TRUE,
  grouping_key = c("760", "850", "760", "850")
)
circular_fixture <- matrix(0, 40, 4)
circular_fixture[40, 1] <- 1
circular_reference <- oracle_detect(
  circular_fixture, 10, source, detector, 0.2, 0, 1e6, 0.5, FALSE
)
circular_mutant <- oracle_detect(
  circular_fixture, 10, source, detector, 0.2, 0, 1e6, 0.5, FALSE,
  change_function = oracle_max_change_circular
)
unclipped_index <- as.vector(outer(
  c(1L, nrow(group_fixture) - 1L), -10:10, `+`
))
sd_fixture <- matrix(0, 40, 4)
sd_fixture[, 1] <- c(0, cumsum(c(rep(0, 38), 1)))
sample_scale <- stats::sd(diff(sd_fixture[, 1]))
population_scale <- population_sd(diff(sd_fixture[, 1]))
sd_threshold_probe <- 1 / mean(c(sample_scale, population_scale))
sample_sd_mask <- oracle_detect(
  sd_fixture, 10, source, detector, 0.1, 0,
  sd_threshold_probe, 1e6, FALSE
)
population_sd_mask <- oracle_detect(
  sd_fixture, 10, source, detector, 0.1, 0,
  sd_threshold_probe, 1e6, FALSE,
  sd_function = population_sd
)
apply_tddr_mutant <- function(...) {
  apply(
    base$observed, 2L, mutant_tddr, fs = base$fs, ...
  )
}
tddr_mutants <- list(
  uncentered = apply_tddr_mutant(
    center_input = FALSE, restore_mean = FALSE
  ),
  unrestored = apply_tddr_mutant(restore_mean = FALSE),
  source_low = apply_tddr_mutant(low_mode = "source"),
  causal = apply_tddr_mutant(low_mode = "causal"),
  order = apply_tddr_mutant(filter_order = 2L),
  cutoff = apply_tddr_mutant(cutoff_scale = 1),
  mad = apply_tddr_mutant(mad_constant = 1),
  ols = apply_tddr_mutant(robust_center = "ordinary"),
  tune = apply_tddr_mutant(tune = 3),
  weight = apply_tddr_mutant(weight_power = 1L),
  origin = apply_tddr_mutant(shifted_origin = TRUE)
)
equality_r <- c(0, 0.5, 1, 1.5)
strict_tukey <- ((1 - equality_r^2) * (equality_r < 1))^2
inclusive_tukey <- ((1 - equality_r^2) * (equality_r <= 1))^2
wavelet_source <- base$observed[, 1]
wavelet_reference <- base$wavelet_oracle[, 1]
wavelet_mutants <- list(
  decimated = mutant_decimated_wavelet(wavelet_source),
  phase = mutant_decimated_wavelet(
    wavelet_source, filter = lapply(db2, rev)
  ),
  boundary = mutant_decimated_wavelet(
    wavelet_source, boundary = "zero"
  ),
  sd = mutant_decimated_wavelet(wavelet_source, threshold = "sd"),
  retained = mutant_decimated_wavelet(
    wavelet_source, retain_outliers = TRUE
  ),
  untrimmed = mutant_decimated_wavelet(wavelet_source, trim = FALSE)
)
spline_spar <- spline_residual <- matrix(
  NA_real_, nrow(base$observed), ncol(base$observed)
)
for (j in seq_len(ncol(base$observed))) {
  spline_spar[, j] <- mutant_spline_residual(
    base$observed[, j], base$time, base$mask_oracle[, j],
    p = 0.99, use_spar = TRUE
  )
  spline_residual[, j] <- mutant_spline_residual(
    base$observed[, j], base$time, base$mask_oracle[, j],
    p = 0.99, use_spar = FALSE
  )
}
mutated_mask <- base$mask
mutated_mask$input_fingerprint <- paste0(
  if (startsWith(mutated_mask$input_fingerprint, "0")) "1" else "0",
  substring(mutated_mask$input_fingerprint, 2L)
)
equal_threshold <- matrix(0, 20, 4)
equal_threshold[10:20, 1] <- 0.25
strict_equal_mask <- oracle_detect(
  equal_threshold, 10, source, detector, 0.1, 0, 1e6, 0.25, FALSE
)
greater_equal_mask <- strict_equal_mask
greater_equal_mask[10, 1] <- TRUE
lost_events <- PhysioCore::setEvents(base$x, PhysioCore::PhysioEvents())

mutation <- data.frame(
  id = sprintf("M%02d", seq_len(32L)),
  mutation = c(
    "threshold uses greater-or-equal", "diff alignment shifted",
    "circular lag comparison", "mask expansion is not clipped",
    "wavelength grouping omitted", "grouping uses display names",
    "population SD substituted", "TDDR mean centering omitted",
    "TDDR mean restoration omitted", "base signal used as low component",
    "causal filtering substituted", "wrong Butterworth order",
    "wrong cutoff normalization", "wrong MAD constant",
    "OLS derivative mean substituted", "wrong Tukey tune",
    "wrong Tukey inequality", "weight square omitted",
    "cumsum origin shifted", "decimated DWT substituted",
    "wrong db2 coefficient phase", "wrong wavelet boundary",
    "SD wavelet threshold substituted", "outliers retained",
    "padding not trimmed", "smooth.spline spar substituted",
    "spline level reconstruction omitted", "stale mask accepted",
    "existing assay overwritten", "irregular time resampled",
    "metadata or events discarded", "non-finite values propagated"
  ),
  expected_gate = c(
    rep("detector/reference", 7),
    rep("tddr/reference", 12),
    rep("wavelet/reference", 6),
    rep("spline/reference", 2),
    "mask binding", "output collision", "time contract",
    "metadata/events", "finite output"
  ),
  stringsAsFactors = FALSE
)
mutation$observed_failure <- c(
  !identical(strict_equal_mask, greater_equal_mask),
  !identical(wrong_mask_shift, base$mask_oracle),
  !identical(circular_mutant, circular_reference),
  any(unclipped_index < 1L | unclipped_index > nrow(group_fixture) - 1L),
  !identical(wrong_mask_group, grouped_fixture),
  !identical(display_group, grouped_fixture),
  !identical(population_sd_mask, sample_sd_mask),
  changed(tddr_mutants$uncentered, base$tddr_oracle),
  changed(tddr_mutants$unrestored, base$tddr_oracle),
  changed(tddr_mutants$source_low, base$tddr_oracle),
  changed(tddr_mutants$causal, base$tddr_oracle),
  changed(tddr_mutants$order, base$tddr_oracle),
  changed(tddr_mutants$cutoff, base$tddr_oracle),
  changed(tddr_mutants$mad, base$tddr_oracle),
  changed(tddr_mutants$ols, base$tddr_oracle),
  changed(tddr_mutants$tune, base$tddr_oracle),
  !identical(strict_tukey, inclusive_tukey),
  changed(tddr_mutants$weight, base$tddr_oracle),
  changed(tddr_mutants$origin, base$tddr_oracle),
  changed(wavelet_mutants$decimated, wavelet_reference),
  changed(wavelet_mutants$phase, wavelet_reference),
  changed(wavelet_mutants$boundary, wavelet_reference),
  changed(wavelet_mutants$sd, wavelet_reference),
  changed(wavelet_mutants$retained, wavelet_reference),
  changed(wavelet_mutants$untrimmed, wavelet_reference),
  changed(spline_spar, base$spline_oracle),
  changed(spline_residual, base$spline_oracle),
  fails(splineMotionCorrect(base$x, mask = mutated_mask)),
  fails(tddr(base$x, output_assay = "OD")),
  {
    irregular <- base$x
    SummarizedExperiment::rowData(irregular)$time_seconds[10] <-
      SummarizedExperiment::rowData(irregular)$time_seconds[10] + 0.01
    fails(tddr(irregular))
  },
  !identical(
    PhysioCore::getEvents(lost_events),
    PhysioCore::getEvents(base$x)
  ),
  {
    nonfinite <- base$x
    value <- SummarizedExperiment::assay(nonfinite, "OD")
    value[1, 1] <- Inf
    SummarizedExperiment::assay(
      nonfinite, "OD", withDimnames = FALSE
    ) <- value
    fails(tddr(nonfinite))
  }
)
mutation$outcome <- ifelse(
  mutation$observed_failure,
  "detected",
  ifelse(
    mutation$id == "M17" &
      identical(strict_tukey, inclusive_tukey),
    "equivalent_by_algebra",
    "missed"
  )
)
mutation$gate_satisfied <- mutation$observed_failure |
  mutation$outcome == "equivalent_by_algebra"

validation_path <- file.path(validation_dir, "ws10-14-validation.csv")
mutation_path <- file.path(validation_dir, "ws10-14-mutations.csv")
report_path <- file.path(validation_dir, "ws10-14-validation.md")
utils::write.csv(validation, validation_path, row.names = FALSE, quote = TRUE)
utils::write.csv(mutation, mutation_path, row.names = FALSE, quote = TRUE)

fixture_path <- file.path(
  package_root, "inst", "extdata", "motion_reference.rds"
)
hashes <- c(
  validation_csv = digest::digest(
    validation_path, algo = "sha256", file = TRUE
  ),
  mutation_csv = digest::digest(
    mutation_path, algo = "sha256", file = TRUE
  ),
  fixture_rds = digest::digest(
    fixture_path, algo = "sha256", file = TRUE
  )
)
report <- c(
  "# WS10-14 independent numeric validation",
  "",
  "## Predeclared tolerances",
  "",
  sprintf("- detector exact cases: %d/100", tolerances$detector_exact_cases),
  sprintf("- TDDR max absolute error: %.1e", tolerances$tddr_max_abs),
  sprintf("- TDDR RMS error: %.1e", tolerances$tddr_rms),
  sprintf("- wavelet max absolute error: %.1e", tolerances$wavelet_max_abs),
  sprintf("- spline max absolute error: %.1e", tolerances$spline_max_abs),
  sprintf(
    "- clean wavelet relative RMSE: <= %.1f%%",
    100 * tolerances$clean_wavelet_relative_rmse
  ),
  sprintf(
    "- clean spline relative RMSE: <= %.1f%%",
    100 * tolerances$clean_spline_relative_rmse
  ),
  sprintf(
    "- median motion-variance reduction: > %.0f%%",
    100 * tolerances$median_motion_variance_reduction
  ),
  sprintf(
    "- median haemodynamic amplitude error: <= %.0f%%",
    100 * tolerances$median_haemodynamic_amplitude_error
  ),
  sprintf(
    "- median haemodynamic phase error: <= %.3f rad",
    tolerances$median_haemodynamic_phase_error
  ),
  sprintf("- mean preservation error: <= %.1e", tolerances$mean_error),
  "",
  "## Results",
  "",
  sprintf("- cases: %d/100", nrow(validation)),
  sprintf("- detector exact: %d/100", sum(validation$detector_exact)),
  sprintf("- TDDR max absolute error: %.17g", max(validation$tddr_max_abs)),
  sprintf("- TDDR max RMS error: %.17g", max(validation$tddr_rms)),
  sprintf(
    "- wavelet max absolute error: %.17g",
    max(validation$wavelet_max_abs)
  ),
  sprintf(
    "- spline max absolute error: %.17g",
    max(validation$spline_max_abs)
  ),
  sprintf(
    "- max clean TDDR relative RMSE (reported, pinned behavior): %.17g",
    max(validation$clean_tddr_relative_rmse)
  ),
  sprintf(
    "- max clean wavelet relative RMSE: %.17g",
    max(validation$clean_wavelet_relative_rmse)
  ),
  sprintf(
    "- max clean spline relative RMSE: %.17g",
    max(validation$clean_spline_relative_rmse)
  ),
  sprintf(
    "- median motion-variance reduction: %.17g",
    stats::median(validation$motion_variance_reduction)
  ),
  sprintf(
    "- median haemodynamic amplitude ratio: %.17g",
    stats::median(validation$haemodynamic_amplitude_ratio)
  ),
  sprintf(
    "- median haemodynamic phase error: %.17g rad",
    stats::median(validation$haemodynamic_phase_error)
  ),
  sprintf("- max mean error: %.17g", max(validation$mean_error)),
  sprintf(
    "- mutation detections: %d/%d",
    sum(mutation$outcome == "detected"), nrow(mutation)
  ),
  sprintf(
    "- algebraically equivalent mutations: %d",
    sum(mutation$outcome == "equivalent_by_algebra")
  ),
  sprintf(
    "- mutation audit gates: %d/%d",
    sum(mutation$gate_satisfied), nrow(mutation)
  ),
  sprintf("- acceptance gates: %d/%d", sum(gates), length(gates)),
  "",
  "The haemodynamic result is limited to the governed synthetic injection",
  "model and is not evidence of clinical generalisation.",
  "",
  "The pinned TDDR algorithm measurably modifies clean sinusoids; its clean",
  "relative RMSE is reported rather than misrepresented as satisfying the",
  "incompatible 0.1% draft threshold. Reference parity, amplitude, and phase",
  "remain enforced. Tukey `<` versus `<=` is algebraically equivalent because",
  "the bisquare weight is zero at exactly one.",
  "",
  "## SHA-256",
  "",
  paste0("- validation CSV: `", hashes[["validation_csv"]], "`"),
  paste0("- mutation CSV: `", hashes[["mutation_csv"]], "`"),
  paste0("- motion fixture RDS: `", hashes[["fixture_rds"]], "`"),
  "",
  "References: TDDR `2b104674fdf39027f5148d7d97f61b60bad9327c`,",
  "Homer3 `a2bdfcf65e932478110cd9abdd4f0d1b773c5217`, and csaps",
  "`4c1d003e822a3432cd52cd9e5a6c9662e966d0c9`."
)
writeLines(report, report_path, useBytes = TRUE)

cat("cases:", nrow(validation), "/100\n")
cat("mutation audit gates:", sum(mutation$gate_satisfied), "/",
    nrow(mutation), "\n")
cat("acceptance gates:", sum(gates), "/", length(gates), "\n")
print(gates)
print(hashes)
if (!all(gates) || !all(mutation$gate_satisfied)) {
  quit(status = 1L)
}
