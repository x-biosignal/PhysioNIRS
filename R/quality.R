.nirs_quality_pair_score <- function(values, fun) {
  if (ncol(values) < 2L) {
    stop("Each optical-density pair requires at least two wavelengths",
         call. = FALSE)
  }
  pairs <- utils::combn(seq_len(ncol(values)), 2L, simplify = FALSE)
  scores <- vapply(pairs, function(index) {
    fun(values[, index[[1L]]], values[, index[[2L]]])
  }, numeric(1))
  if (!length(scores) || anyNA(scores) || any(!is.finite(scores))) 0 else {
    min(scores)
  }
}

.nirs_sci_pair <- function(x, y) {
  if (.nirs_near_constant(x) || .nirs_near_constant(y)) return(0)
  value <- suppressWarnings(stats::cor(x, y, method = "pearson"))
  if (length(value) != 1L || !is.finite(value)) 0 else as.numeric(value)
}

#' Calculate the governed fNIRS scalp-coupling index
#'
#' Optical-density measurements are grouped by exact source-detector identity,
#' zero-phase filtered in the cardiac band, and scored by the minimum
#' pairwise wavelength correlation. Scores are copied to every wavelength in
#' the group. The default threshold is a project default, not a clinical
#' cutoff.
#'
#' @param x A governed measurement-level optical-density
#'   `PhysioExperiment`.
#' @param assay_name Exact optical-density assay name.
#' @param l_freq,h_freq Cardiac-band limits in Hz.
#' @param window_seconds Optional positive complete-window duration.
#' @param step_seconds Optional positive window step. `NULL` uses the window.
#' @param threshold Finite pass threshold in `[-1, 1]`.
#' @param order Positive Butterworth filter order.
#'
#' @return A source-bound `nirs_quality` result.
#' @references Pollonini et al. (2014), DOI: 10.1117/1.JBO.19.8.086007.
#' @export
scalpCouplingIndex <- function(
    x,
    assay_name = "OD",
    l_freq = 0.7,
    h_freq = 1.5,
    window_seconds = NULL,
    step_seconds = NULL,
    threshold = 0.8,
    order = 3L) {
  threshold <- .nirs_quality_number(
    threshold, "threshold", lower = -1, upper = 1
  )
  context <- .nirs_quality_context(
    x, assay_name, window_seconds, step_seconds
  )
  if (!identical(context$identity_kind, "measurement") ||
      !identical(context$contract$kind, "optical_density")) {
    stop("SCI requires governed measurement-level optical density",
         call. = FALSE)
  }
  groups <- .nirs_quality_groups(context$identity)
  for (index in groups) {
    wavelength <- context$identity$wavelength_nm[index]
    if (length(unique(wavelength)) < 2L ||
        anyDuplicated(wavelength)) {
      stop(
        "Every SCI source-detector pair requires at least two distinct ",
        "wavelength measurements",
        call. = FALSE
      )
    }
  }
  filtered <- .nirs_quality_filter(
    context, l_freq, h_freq, order
  )
  score <- matrix(
    0, nrow = ncol(context$data), ncol = length(context$windows$id)
  )
  for (w in seq_along(context$windows$id)) {
    rows <- seq.int(
      context$windows$start[[w]], context$windows$end[[w]]
    )
    for (index in groups) {
      group_score <- .nirs_quality_pair_score(
        filtered$value[rows, index, drop = FALSE], .nirs_sci_pair
      )
      score[index, w] <- group_score
    }
  }
  .nirs_quality_result(
    context, "scalp_coupling_index", score, threshold,
    parameters = list(
      l_freq = filtered$l_freq,
      h_freq = filtered$h_freq,
      order = filtered$order,
      window_seconds = if (is.null(window_seconds)) NULL else {
        as.numeric(window_seconds)
      },
      step_seconds = if (is.null(step_seconds)) NULL else {
        as.numeric(step_seconds)
      },
      group_rule = "minimum finite pairwise wavelength correlation",
      zero_variance_score = 0
    ),
    implementation = list(
      package_version = "0.5.0",
      filter = "Butterworth forward/reverse padlen=0",
      reference = .nirs_quality_reference
    )
  )
}

.nirs_hamming <- function(n) {
  if (n < 2L) stop("Periodogram requires at least two samples", call. = FALSE)
  0.54 - 0.46 * cos(2 * pi * seq.int(0, n - 1L) / (n - 1L))
}

.nirs_periodogram <- function(x, fs) {
  x <- as.numeric(x)
  if (length(x) < 4L || anyNA(x) || any(!is.finite(x))) {
    stop("Periodogram input must contain at least four finite samples",
         call. = FALSE)
  }
  x <- x - mean(x)
  window <- .nirs_hamming(length(x))
  transformed <- stats::fft(x * window)
  last <- floor(length(x) / 2L) + 1L
  power <- Mod(transformed[seq_len(last)])^2 /
    (fs * sum(window^2))
  if (length(power) > 2L) {
    interior <- seq.int(2L, length(power) -
      if (length(x) %% 2L == 0L) 1L else 0L)
    power[interior] <- 2 * power[interior]
  }
  frequency <- seq.int(0L, last - 1L) * fs / length(x)
  list(frequency = as.numeric(frequency), power = as.numeric(power))
}

.nirs_full_cross_correlation <- function(x, y) {
  as.numeric(stats::convolve(x, y, type = "open")) / length(x)
}

.nirs_peak_power_pair <- function(x, y, fs) {
  if (.nirs_near_constant(x) || .nirs_near_constant(y)) return(0)
  x <- x / stats::sd(x)
  y <- y / stats::sd(y)
  correlation <- .nirs_full_cross_correlation(x, y)
  power <- .nirs_periodogram(correlation, fs)$power
  value <- max(power)
  if (!is.finite(value)) 0 else as.numeric(value)
}

.nirs_validate_bands <- function(bands, nyquist, arg) {
  if (!is.matrix(bands) || !is.numeric(bands) || is.complex(bands) ||
      ncol(bands) != 2L || nrow(bands) < 1L ||
      anyNA(bands) || any(!is.finite(bands)) ||
      any(bands[, 1L] <= 0) || any(bands[, 1L] >= bands[, 2L]) ||
      any(bands[, 2L] >= nyquist)) {
    stop("`", arg, "` must be a finite two-column matrix of positive ",
         "increasing bands strictly below Nyquist", call. = FALSE)
  }
  bands <- bands[order(bands[, 1L], bands[, 2L]), , drop = FALSE]
  if (nrow(bands) > 1L &&
      any(bands[-1L, 1L] <= bands[-nrow(bands), 2L])) {
    stop("`", arg, "` bands must be disjoint", call. = FALSE)
  }
  unname(bands)
}

.nirs_band_integral <- function(frequency, power, bands) {
  values <- vapply(seq_len(nrow(bands)), function(i) {
    lower <- bands[i, 1L]
    upper <- bands[i, 2L]
    interior <- frequency > lower & frequency < upper
    grid <- c(lower, frequency[interior], upper)
    ordinate <- stats::approx(
      frequency, power, xout = grid, method = "linear",
      rule = 1, ties = "ordered"
    )$y
    if (anyNA(ordinate) || any(!is.finite(ordinate))) {
      stop("Spectral band is outside the periodogram support",
           call. = FALSE)
    }
    sum(
      diff(grid) *
        (utils::head(ordinate, -1L) + utils::tail(ordinate, -1L)) / 2
    )
  }, numeric(1))
  sum(values)
}

.nirs_snr_score <- function(x, fs, cardiac, noise) {
  spectrum <- .nirs_periodogram(x, fs)
  signal_power <- .nirs_band_integral(
    spectrum$frequency, spectrum$power, cardiac
  )
  noise_power <- .nirs_band_integral(
    spectrum$frequency, spectrum$power, noise
  )
  if (!is.finite(signal_power) || signal_power <= 0 ||
      !is.finite(noise_power) || noise_power <= 0) {
    stop("SNR requires strictly positive signal and noise band power",
         call. = FALSE)
  }
  as.numeric(10 * log10(signal_power / noise_power))
}

#' Calculate governed fNIRS signal quality
#'
#' `method = "peak_power"` implements the pinned MNE-NIRS PHOEBE-style
#' cross-correlation peak-power metric. It is not the original PHOEBE
#' equation. `method = "snr"` returns a channel-specific cardiac-to-noise
#' spectral power ratio in decibels.
#'
#' @inheritParams scalpCouplingIndex
#' @param method Exactly `"peak_power"` or `"snr"`.
#' @param window_seconds Positive complete-window duration.
#' @param cardiac_range_hz Finite increasing cardiac band below Nyquist.
#' @param noise_range_hz For SNR, `NULL` or a disjoint two-column band matrix.
#' @param threshold Optional finite threshold. Defaults to `0.1` for peak
#'   power and `0` dB for SNR.
#'
#' @return A source-bound `nirs_quality` result.
#' @references Pollonini et al. (2016), DOI: 10.1364/BOE.7.005104.
#' @export
signalQualityIndex <- function(
    x,
    assay_name = "OD",
    method = c("peak_power", "snr"),
    window_seconds = 10,
    step_seconds = NULL,
    cardiac_range_hz = c(0.7, 1.5),
    noise_range_hz = NULL,
    threshold = NULL,
    order = 3L) {
  method <- if (missing(method)) {
    "peak_power"
  } else {
    .snirf_enum(method, c("peak_power", "snr"), "method")
  }
  if (!is.numeric(cardiac_range_hz) || is.complex(cardiac_range_hz) ||
      !is.null(dim(cardiac_range_hz)) ||
      length(cardiac_range_hz) != 2L || anyNA(cardiac_range_hz) ||
      any(!is.finite(cardiac_range_hz))) {
    stop("`cardiac_range_hz` must be two finite real values",
         call. = FALSE)
  }
  context <- .nirs_quality_context(
    x, assay_name, window_seconds, step_seconds
  )
  cardiac <- matrix(
    .nirs_validate_bands(
      matrix(cardiac_range_hz, nrow = 1L),
      context$fs / 2, "cardiac_range_hz"
    ),
    nrow = 1L
  )
  if (is.null(threshold)) {
    threshold <- if (method == "peak_power") 0.1 else 0
  }
  threshold <- .nirs_quality_number(threshold, "threshold")
  score <- matrix(
    0, nrow = ncol(context$data), ncol = length(context$windows$id)
  )

  if (method == "peak_power") {
    if (!identical(context$identity_kind, "measurement") ||
        !identical(context$contract$kind, "optical_density")) {
      stop("Peak power requires governed measurement-level optical density",
           call. = FALSE)
    }
    groups <- .nirs_quality_groups(context$identity)
    for (index in groups) {
      wavelength <- context$identity$wavelength_nm[index]
      if (length(unique(wavelength)) < 2L ||
          anyDuplicated(wavelength)) {
        stop("Peak power requires distinct wavelengths in every pair",
             call. = FALSE)
      }
    }
    filtered <- .nirs_quality_filter(
      context, cardiac[1L, 1L], cardiac[1L, 2L], order
    )
    for (w in seq_along(context$windows$id)) {
      rows <- seq.int(
        context$windows$start[[w]], context$windows$end[[w]]
      )
      for (index in groups) {
        value <- .nirs_quality_pair_score(
          filtered$value[rows, index, drop = FALSE],
          function(x, y) .nirs_peak_power_pair(x, y, context$fs)
        )
        score[index, w] <- value
      }
    }
    noise <- NULL
  } else {
    order <- .nirs_quality_number(
      order, "order", lower = 1, upper = 20, integer = TRUE
    )
    if (is.null(noise_range_hz)) {
      minimum_positive <- context$fs /
        max(vapply(seq_along(context$windows$id), function(w) {
          context$windows$end[[w]] - context$windows$start[[w]] + 1L
        }, integer(1)))
      candidate <- rbind(
        c(minimum_positive, cardiac[1L, 1L]),
        c(cardiac[1L, 2L], context$fs / 2 -
          .Machine$double.eps * context$fs)
      )
      candidate <- candidate[candidate[, 2L] > candidate[, 1L], ,
                             drop = FALSE]
      if (!nrow(candidate)) {
        stop("Default SNR noise band is empty", call. = FALSE)
      }
      noise <- candidate
    } else {
      noise <- .nirs_validate_bands(
        noise_range_hz, context$fs / 2, "noise_range_hz"
      )
    }
    overlap <- vapply(seq_len(nrow(noise)), function(i) {
      noise[i, 1L] < cardiac[1L, 2L] &&
        noise[i, 2L] > cardiac[1L, 1L]
    }, logical(1))
    if (any(overlap)) {
      stop("`noise_range_hz` must not overlap `cardiac_range_hz`",
           call. = FALSE)
    }
    for (w in seq_along(context$windows$id)) {
      rows <- seq.int(
        context$windows$start[[w]], context$windows$end[[w]]
      )
      for (j in seq_len(ncol(context$data))) {
        score[j, w] <- .nirs_snr_score(
          context$data[rows, j], context$fs, cardiac, noise
        )
      }
    }
  }

  .nirs_quality_result(
    context,
    if (method == "peak_power") "phoebe_peak_power" else "snr_db",
    score, threshold,
    parameters = list(
      method = method,
      cardiac_range_hz = as.numeric(cardiac),
      noise_range_hz = if (is.null(noise)) NULL else unname(noise),
      order = as.integer(order),
      window_seconds = as.numeric(window_seconds),
      step_seconds = if (is.null(step_seconds)) NULL else {
        as.numeric(step_seconds)
      },
      peak_group_rule = if (method == "peak_power") {
        "minimum pair score copied to every wavelength"
      } else {
        NULL
      }
    ),
    implementation = list(
      package_version = "0.5.0",
      periodogram = "one-sided Hamming PSD with interpolated band edges",
      filter = if (method == "peak_power") {
        "Butterworth forward/reverse padlen=0"
      } else {
        NULL
      },
      reference = .nirs_quality_reference
    )
  )
}

.nirs_quality_list <- function(quality, x) {
  if (inherits(quality, "nirs_quality")) {
    quality <- stats::setNames(list(quality), quality$metric)
  } else if (!is.list(quality) || is.object(quality) || !length(quality) ||
             is.null(names(quality)) || anyNA(names(quality)) ||
             any(!nzchar(names(quality))) || anyDuplicated(names(quality))) {
    stop("`quality` must be one result or a non-empty uniquely named list",
         call. = FALSE)
  }
  metrics <- character(length(quality))
  for (i in seq_along(quality)) {
    quality[[i]] <- .nirs_validate_quality(quality[[i]], x)
    metrics[[i]] <- quality[[i]]$metric
  }
  if (anyDuplicated(metrics)) {
    stop("Quality results must have distinct metric identities",
         call. = FALSE)
  }
  quality
}

.nirs_quality_fail_ids <- function(result, index) {
  windows <- colnames(result$pass)[!result$pass[index, ]]
  if (!length(windows)) "" else paste0(result$metric, ":", windows,
                                        collapse = ",")
}

.nirs_quality_update_drop_metadata <- function(out, x, keep, identity_kind) {
  metadata <- S4Vectors::metadata(out)
  if (identity_kind == "measurement") {
    measurement <- measurementList(x)[keep, , drop = FALSE]
    measurement$measurement_index <- seq_len(nrow(measurement))
    metadata$snirf$measurement_list <- measurement
    columns <- SummarizedExperiment::colData(out)
    if ("measurement_index" %in% names(columns)) {
      columns$measurement_index <- seq_len(nrow(columns))
    }
    SummarizedExperiment::colData(out) <- columns
  } else if (is.list(metadata$nirs$mbll)) {
    if (length(metadata$nirs$mbll$condition_numbers) == length(keep)) {
      metadata$nirs$mbll$condition_numbers <-
        metadata$nirs$mbll$condition_numbers[keep]
    }
    mapping <- metadata$nirs$mbll$input_to_pair_mapping
    if (is.list(mapping) && length(mapping) == length(keep)) {
      mapping <- mapping[keep]
      for (i in seq_along(mapping)) {
        mapping[[i]]$output_pair_index <- as.integer(i)
      }
      metadata$nirs$mbll$input_to_pair_mapping <- mapping
    }
  }
  S4Vectors::metadata(out) <- metadata
  out
}

#' Mark or drop channels using governed NIRS quality results
#'
#' Window decisions are collapsed conservatively per metric. Measurement-level
#' optical-density drops propagate to the complete source-detector pair.
#' Inputs and quality results are fingerprint-bound and are never modified.
#'
#' @param x The exact governed source `PhysioExperiment`.
#' @param quality One `nirs_quality` result or a non-empty named list of
#'   distinct results.
#' @param action Exactly `"mark"` or `"drop"`.
#' @param require_all Whether all metrics, rather than any metric, must pass.
#'
#' @return A cloned marked or consistently subset `PhysioExperiment`.
#' @export
pruneChannels <- function(
    x,
    quality,
    action = c("mark", "drop"),
    require_all = TRUE) {
  action <- if (missing(action)) {
    "mark"
  } else {
    .snirf_enum(action, c("mark", "drop"), "action")
  }
  require_all <- .nirs_motion_flag(require_all, "require_all")
  qualities <- .nirs_quality_list(quality, x)
  first <- qualities[[1L]]
  for (result in qualities[-1L]) {
    if (!identical(result$source_fingerprint, first$source_fingerprint) ||
        !identical(result$identity_fingerprint,
                   first$identity_fingerprint) ||
        !identical(result$channel_id, first$channel_id) ||
        !identical(result$identity_kind, first$identity_kind)) {
      stop("All quality results must bind the same ordered source identity",
           call. = FALSE)
    }
  }
  metric_pass <- vapply(qualities, function(result) {
    apply(result$pass, 1L, all)
  }, logical(length(first$channel_id)))
  if (is.null(dim(metric_pass))) {
    metric_pass <- matrix(
      metric_pass, ncol = 1L,
      dimnames = list(first$channel_id, names(qualities))
    )
  } else {
    rownames(metric_pass) <- first$channel_id
    colnames(metric_pass) <- names(qualities)
  }
  final_pass <- if (require_all) {
    apply(metric_pass, 1L, all)
  } else {
    apply(metric_pass, 1L, any)
  }
  pair_propagated <- rep(FALSE, length(final_pass))
  if (action == "drop" &&
      identical(first$identity_kind, "measurement")) {
    key <- paste(first$source_index, first$detector_index, sep = ":")
    for (index in split(seq_along(key), factor(key, levels = unique(key)))) {
      value <- all(final_pass[index])
      pair_propagated[index] <- final_pass[index] != value
      final_pass[index] <- value
    }
  }
  quality_fingerprints <- vapply(
    qualities,
    function(result) attr(result, "quality_fingerprint", exact = TRUE),
    character(1)
  )
  operation_fingerprint <- .nirs_sha256(list(
    source = first$source_fingerprint,
    quality = quality_fingerprints,
    action = action,
    require_all = require_all,
    final_pass = unname(final_pass)
  ))
  record <- list(
    schema_version = .nirs_quality_schema,
    action = action,
    require_all = require_all,
    quality_fingerprints = quality_fingerprints,
    metric_pass = unname(metric_pass),
    metric_names = colnames(metric_pass),
    retained_ids = first$channel_id[final_pass],
    dropped_ids = first$channel_id[!final_pass],
    pair_propagated = pair_propagated,
    source_fingerprint = first$source_fingerprint,
    operation_fingerprint = operation_fingerprint
  )

  if (action == "mark") {
    columns <- SummarizedExperiment::colData(x)
    new_names <- c(
      "nirs_quality_pass", "nirs_quality_metrics",
      "nirs_quality_failed_windows", "nirs_quality_fingerprint"
    )
    if (any(new_names %in% names(columns))) {
      stop("Governed NIRS quality columns already exist", call. = FALSE)
    }
    failures <- vapply(seq_along(final_pass), function(i) {
      values <- vapply(
        qualities, .nirs_quality_fail_ids, character(1), index = i
      )
      paste(values[nzchar(values)], collapse = ";")
    }, character(1))
    columns$nirs_quality_pass <- as.logical(final_pass)
    columns$nirs_quality_metrics <- rep(
      paste(colnames(metric_pass), collapse = ","), length(final_pass)
    )
    columns$nirs_quality_failed_windows <- failures
    columns$nirs_quality_fingerprint <- rep(
      operation_fingerprint, length(final_pass)
    )
    out <- x
    SummarizedExperiment::colData(out) <- columns
  } else {
    if (!any(final_pass)) {
      stop("Channel pruning would remove every channel", call. = FALSE)
    }
    contracts <- S4Vectors::metadata(x)$nirs$assays
    corrected <- vapply(contracts, function(contract) {
      isTRUE(contract$motion_corrected) ||
        isTRUE(contract$short_separation_corrected)
    }, logical(1))
    if (any(corrected)) {
      stop(
        "Dropping channels with stored correction provenance is unsupported; ",
        "mark first or recompute corrections after pruning",
        call. = FALSE
      )
    }
    out <- x[, final_pass, drop = FALSE]
    out <- .nirs_quality_update_drop_metadata(
      out, x, final_pass, first$identity_kind
    )
    if (ncol(out) < 1L ||
        any(vapply(
          SummarizedExperiment::assays(out),
          function(value) ncol(value) != ncol(out),
          logical(1)
        ))) {
      stop("Channel pruning produced inconsistent assay dimensions",
           call. = FALSE)
    }
  }
  metadata <- S4Vectors::metadata(out)
  if (is.null(metadata$nirs$quality)) metadata$nirs$quality <- list()
  metadata$nirs$quality[[length(metadata$nirs$quality) + 1L]] <- record
  S4Vectors::metadata(out) <- metadata
  out <- .nirs_append_step(
    out,
    "pruneChannels",
    params = list(
      implementation_version = "0.5.0",
      action = action,
      require_all = require_all,
      quality_fingerprints = quality_fingerprints,
      retained_ids = record$retained_ids,
      dropped_ids = record$dropped_ids,
      pair_propagated = pair_propagated,
      source_fingerprint = first$source_fingerprint,
      operation_fingerprint = operation_fingerprint
    ),
    input_assay = paste(
      unique(vapply(qualities, `[[`, character(1), "assay_name")),
      collapse = ","
    ),
    output_assay = paste(
      SummarizedExperiment::assayNames(out), collapse = ","
    )
  )
  methods::validObject(out)
  out
}
