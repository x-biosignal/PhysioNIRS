.nirs_resolve_short_map <- function(x, short, context) {
  if (is.null(short)) {
    short <- .nirs_make_short_map(context, 0.01)
  }
  .nirs_validate_short_map(short, x, context)
}

.nirs_short_compatible <- function(context, short_indices, target_index) {
  if (context$identity$identity_kind[[1L]] == "measurement") {
    wavelength <- context$identity$wavelength_nm
    short_indices[wavelength[short_indices] == wavelength[[target_index]]]
  } else {
    short_indices
  }
}

.nirs_nearest_short <- function(map, candidates, target_index) {
  target <- unlist(
    map[target_index, c("midpoint_x_m", "midpoint_y_m", "midpoint_z_m")],
    use.names = FALSE
  )
  short_midpoint <- as.matrix(map[
    candidates,
    c("midpoint_x_m", "midpoint_y_m", "midpoint_z_m"),
    drop = FALSE
  ])
  distance <- sqrt(rowSums(sweep(short_midpoint, 2L, target, "-")^2))
  minimum <- min(distance)
  tolerance <- max(1e-12, 1e-12 * max(1, minimum))
  tied <- which(abs(distance - minimum) <= tolerance)
  list(
    index = candidates[[tied[[1L]]]],
    distance_m = minimum,
    tied = as.integer(candidates[tied])
  )
}

.nirs_lag_rows <- function(n, lag_samples) {
  if (lag_samples > 0L) {
    list(
      target = seq.int(lag_samples + 1L, n),
      predictor = seq_len(n - lag_samples),
      unchanged = seq_len(lag_samples)
    )
  } else if (lag_samples < 0L) {
    amount <- -lag_samples
    list(
      target = seq_len(n - amount),
      predictor = seq.int(amount + 1L, n),
      unchanged = seq.int(n - amount + 1L, n)
    )
  } else {
    list(
      target = seq_len(n),
      predictor = seq_len(n),
      unchanged = integer()
    )
  }
}

.nirs_short_svd_fit <- function(y, x, condition_limit) {
  design <- cbind(intercept = 1, x)
  decomposition <- svd(
    design,
    nu = min(nrow(design), ncol(design)),
    nv = ncol(design)
  )
  tolerance <- max(dim(design)) * max(decomposition$d) *
    .Machine$double.eps
  rank <- sum(decomposition$d > tolerance)
  condition <- max(decomposition$d) / min(decomposition$d)
  if (rank < ncol(design) || !is.finite(condition) ||
      condition > condition_limit) {
    stop(
      "Short-separation design is rank-deficient or ill-conditioned ",
      "(rank=", rank, "/", ncol(design), ", condition=",
      format(condition, digits = 8L), ")",
      call. = FALSE
    )
  }
  if (nrow(design) - rank < 1L) {
    stop("Short-separation regression requires residual degrees of freedom",
         call. = FALSE)
  }
  coefficient <- decomposition$v %*%
    ((crossprod(decomposition$u, y)) / decomposition$d)
  coefficient <- as.numeric(coefficient)
  names(coefficient) <- colnames(design)
  if (anyNA(coefficient) || any(!is.finite(coefficient))) {
    stop("Short-separation regression produced non-finite coefficients",
         call. = FALSE)
  }
  list(
    coefficient = coefficient,
    rank = as.integer(rank),
    condition = as.numeric(condition),
    tolerance = as.numeric(tolerance)
  )
}

#' Regress governed short-separation nuisance signals
#'
#' Fits a centred short-channel regression separately for each long channel.
#' Short channels are copied unchanged. Positive lag means that the short
#' signal at `t - lag_seconds` predicts the long signal at `t`.
#'
#' @details For a target `y` and selected short signals `X`, the fitted
#'   intercept is retained and only `(X - colMeans(X)) beta` is subtracted.
#'   Regression uses a deterministic SVD solve and rejects rank-deficient or
#'   over-limit designs. Short columns and lag-edge rows remain unchanged.
#'
#' @param x A governed OD or pair-collapsed haemoglobin
#'   `PhysioExperiment`.
#' @param assay_name Exact source assay name.
#' @param short A compatible `nirs_short_channels` map returned by
#'   `identifyShortChannels()`, or `NULL` to use the strict 10-mm default.
#' @param method Exactly `"nearest"` or `"all"`.
#' @param lag_seconds A finite lag that corresponds to an exact sample count.
#' @param output_assay Exact name for the added corrected assay.
#' @param condition_limit Finite regression condition-number limit greater
#'   than one.
#'
#' @return A clone of `x` with one short-separation-corrected assay.
#' @references Saager and Berger (2005), DOI: 10.1364/JOSAA.22.001874.
#' @export
shortSeparationRegress <- function(
    x,
    assay_name = "OD",
    short = NULL,
    method = c("nearest", "all"),
    lag_seconds = 0,
    output_assay = paste0(assay_name, "_ssr"),
    condition_limit = 1e10) {
  method <- if (missing(method)) {
    "nearest"
  } else {
    .snirf_enum(method, c("nearest", "all"), "method")
  }
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  condition_limit <- .nirs_motion_scalar(
    condition_limit, "condition_limit",
    lower = 1 + .Machine$double.eps
  )
  lag_seconds <- .nirs_motion_scalar(lag_seconds, "lag_seconds")
  context <- .nirs_short_context(
    x, assay_name, min_samples = 3L, require_contract = TRUE
  )
  if (output_assay %in% SummarizedExperiment::assayNames(x)) {
    stop("`output_assay` already exists: ", output_assay, call. = FALSE)
  }
  map <- .nirs_resolve_short_map(x, short, context)
  short_indices <- which(map$is_short)
  target_indices <- which(!map$is_short)
  if (!length(short_indices) || !length(target_indices)) {
    stop(
      "Short-separation regression requires at least one short and one long ",
      "channel; distance range is [",
      format(min(map$distance_m), digits = 6L), ", ",
      format(max(map$distance_m), digits = 6L),
      "] m at threshold ",
      format(attr(map, "threshold_m"), digits = 6L), " m",
      call. = FALSE
    )
  }

  lag_exact <- lag_seconds * context$fs
  lag_rounded <- round(lag_exact)
  if (abs(lag_exact - lag_rounded) >
      sqrt(.Machine$double.eps) * max(1, abs(lag_exact)) ||
      abs(lag_rounded) > .Machine$integer.max) {
    stop("`lag_seconds` must correspond to an exact whole sample count",
         call. = FALSE)
  }
  lag_samples <- as.integer(lag_rounded)
  if (abs(lag_samples) >= nrow(context$data) - 2L) {
    stop("Absolute lag must be smaller than `n_time - 2` samples",
         call. = FALSE)
  }
  rows <- .nirs_lag_rows(nrow(context$data), lag_samples)
  corrected <- context$data
  diagnostics <- vector("list", length(target_indices))

  for (i in seq_along(target_indices)) {
    target <- target_indices[[i]]
    compatible <- .nirs_short_compatible(context, short_indices, target)
    if (!length(compatible)) {
      required <- if (context$identity$identity_kind[[1L]] == "measurement") {
        paste0("wavelength ", context$identity$wavelength_nm[[target]], " nm")
      } else {
        paste0("assay ", context$assay_name)
      }
      stop(
        "No compatible short channel for ",
        context$identity$channel_id[[target]], " (", required, ")",
        call. = FALSE
      )
    }
    nearest <- .nirs_nearest_short(map, compatible, target)
    selected <- if (method == "nearest") nearest$index else compatible
    predictors <- context$data[
      rows$predictor, selected, drop = FALSE
    ]
    if (any(vapply(
      seq_len(ncol(predictors)),
      function(j) {
        scale <- max(1, max(abs(predictors[, j])))
        stats::sd(predictors[, j]) <=
          sqrt(.Machine$double.eps) * scale
      },
      logical(1)
    ))) {
      stop(
        "Compatible short signal has zero or near-zero variance for target ",
        context$identity$channel_id[[target]],
        call. = FALSE
      )
    }
    colnames(predictors) <- context$identity$channel_id[selected]
    target_values <- context$data[rows$target, target]
    fit <- .nirs_short_svd_fit(
      target_values, predictors, condition_limit
    )
    centred <- sweep(predictors, 2L, colMeans(predictors), "-")
    value <- target_values -
      as.numeric(centred %*% fit$coefficient[-1L])
    if (anyNA(value) || any(!is.finite(value))) {
      stop("Short-separation correction produced non-finite values",
           call. = FALSE)
    }
    corrected[rows$target, target] <- value
    diagnostics[[i]] <- list(
      target_index = as.integer(target),
      target_id = context$identity$channel_id[[target]],
      compatible_short_indices = as.integer(compatible),
      selected_short_indices = as.integer(selected),
      selected_short_ids = context$identity$channel_id[selected],
      nearest_distance_m = nearest$distance_m,
      nearest_tied_indices = nearest$tied,
      coefficients = fit$coefficient,
      rank = fit$rank,
      condition_number = fit$condition,
      rank_tolerance = fit$tolerance
    )
  }
  dimnames(corrected) <- dimnames(context$data)
  out <- x
  SummarizedExperiment::assay(out, output_assay, withDimnames = FALSE) <-
    corrected
  source_contract <- context$contract
  output_contract <- source_contract
  output_contract$source_assay <- context$assay_name
  output_contract$method <- "short_separation_regression"
  output_contract$short_separation_corrected <- TRUE
  out <- .nirs_set_assay_contract(out, output_assay, output_contract)
  .nirs_append_step(
    out,
    "shortSeparationRegress",
    params = list(
      implementation_version = "0.4.0",
      reference = .nirs_short_reference,
      source_fingerprint = context$source_fingerprint,
      short_map_fingerprint = .nirs_sha256(map),
      probe_fingerprint = attr(map, "probe_fingerprint"),
      threshold_m = attr(map, "threshold_m"),
      method = method,
      lag_seconds = lag_seconds,
      lag_samples = lag_samples,
      sampling_rate_hz = context$fs,
      overlap_target_rows = as.integer(range(rows$target)),
      overlap_predictor_rows = as.integer(range(rows$predictor)),
      unchanged_edge_rows = as.integer(rows$unchanged),
      condition_limit = condition_limit,
      diagnostics = diagnostics
    ),
    input_assay = context$assay_name,
    output_assay = output_assay
  )
}
