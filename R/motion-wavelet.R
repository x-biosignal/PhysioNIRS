.nirs_db2_filters <- list(
  low = c(
    0.48296291314453416,
    0.83651630373780794,
    0.22414386804201339,
    -0.12940952255126037
  ),
  high = c(
    -0.12940952255126037,
    -0.22414386804201339,
    0.83651630373780794,
    -0.48296291314453416
  )
)

.nirs_db2_dwt_periodic <- function(x) {
  n <- length(x)
  if (n < 4L || n %% 2L != 0L) {
    stop("Periodic db2 DWT requires an even length of at least four",
         call. = FALSE)
  }
  half <- n %/% 2L
  low <- high <- numeric(half)
  for (k in seq_len(half)) {
    index <- ((2L * (k - 1L) + seq_along(.nirs_db2_filters$low) -
      2L) %% n) + 1L
    low[[k]] <- sum(.nirs_db2_filters$low * x[index])
    high[[k]] <- sum(.nirs_db2_filters$high * x[index])
  }
  list(approximation = low, detail = high)
}

.nirs_db2_idwt_periodic <- function(approximation, detail) {
  if (!is.numeric(approximation) || !is.numeric(detail) ||
      length(approximation) != length(detail) ||
      length(approximation) < 2L) {
    stop("Periodic db2 inverse requires equal coefficient vectors",
         call. = FALSE)
  }
  half <- length(approximation)
  n <- 2L * half
  out <- numeric(n)
  for (k in seq_len(half)) {
    index <- ((2L * (k - 1L) + seq_along(.nirs_db2_filters$low) -
      2L) %% n) + 1L
    out[index] <- out[index] +
      approximation[[k]] * .nirs_db2_filters$low +
      detail[[k]] * .nirs_db2_filters$high
  }
  out
}

.nirs_shift_invariant_db2 <- function(x, min_level) {
  n <- length(x)
  exponent <- log2(n)
  if (abs(exponent - round(exponent)) > 1e-12) {
    stop("Shift-invariant db2 transform requires power-of-two length",
         call. = FALSE)
  }
  exponent <- as.integer(round(exponent))
  depth <- exponent - min_level
  if (depth < 2L) {
    stop(
      "Signal is too short for `min_level`; padded length must be at least ",
      2^(min_level + 2L),
      call. = FALSE
    )
  }
  coefficients <- matrix(0, nrow = n, ncol = depth + 1L)
  coefficients[, 1L] <- x
  for (level in 0:(depth - 1L)) {
    n_block <- 2^level
    block_length <- n %/% n_block
    half <- block_length %/% 2L
    for (block in 0:(n_block - 1L)) {
      index <- block * block_length + seq_len(block_length)
      signal <- coefficients[index, 1L]
      shifted <- c(signal[[block_length]], signal[-block_length])
      ordinary <- .nirs_db2_dwt_periodic(signal)
      shifted_result <- .nirs_db2_dwt_periodic(shifted)
      coefficients[index[seq_len(half)], 1L] <-
        ordinary$approximation
      coefficients[index[half + seq_len(half)], 1L] <-
        shifted_result$approximation
      coefficients[index[seq_len(half)], level + 2L] <-
        ordinary$detail
      coefficients[index[half + seq_len(half)], level + 2L] <-
        shifted_result$detail
    }
  }
  coefficients
}

.nirs_inverse_shift_invariant_db2 <- function(coefficients) {
  n <- nrow(coefficients)
  depth <- ncol(coefficients) - 1L
  approximation <- coefficients[, 1L]
  for (level in seq.int(depth - 1L, 0L)) {
    n_block <- 2^level
    block_length <- n %/% n_block
    half <- block_length %/% 2L
    for (block in 0:(n_block - 1L)) {
      index <- block * block_length + seq_len(block_length)
      first <- seq_len(half)
      second <- half + seq_len(half)
      ordinary <- .nirs_db2_idwt_periodic(
        approximation[index[first]],
        coefficients[index[first], level + 2L]
      )
      shifted <- .nirs_db2_idwt_periodic(
        approximation[index[second]],
        coefficients[index[second], level + 2L]
      )
      unshifted <- c(shifted[-1L], shifted[[1L]])
      approximation[index] <- (ordinary + unshifted) / 2
    }
  }
  as.numeric(approximation)
}

.nirs_wavelet_channel <- function(x, iqr, min_level) {
  original_length <- length(x)
  exponent <- ceiling(log2(original_length))
  padded_length <- 2^exponent
  if (padded_length < 2^(min_level + 2L)) {
    stop(
      "Signal is too short for `min_level`; requires padded length ",
      2^(min_level + 2L), " but found ", padded_length,
      call. = FALSE
    )
  }
  padded <- numeric(padded_length)
  padded[seq_len(original_length)] <- x
  dc <- mean(padded)
  padded <- padded - dc
  coefficients <- .nirs_shift_invariant_db2(padded, min_level)
  depth <- ncol(coefficients) - 1L
  signal_length <- original_length
  rejected <- integer(depth)
  if (depth > 1L) {
    for (level in seq_len(depth - 1L)) {
      signal_length <- floor(signal_length / 2)
      n_block <- 2^level
      block_length <- padded_length %/% n_block
      for (block in 0:(n_block - 1L)) {
        index <- block * block_length + seq_len(block_length)
        valid_count <- min(signal_length, block_length)
        reference <- coefficients[
          index[seq_len(valid_count)], level + 1L
        ]
        quartile <- as.numeric(stats::quantile(
          reference,
          probs = c(0.25, 0.5, 0.75),
          type = 7,
          names = FALSE
        ))
        # Homer3's common coefficient normalisation cancels from both the
        # coefficients and these IQR fences, so threshold membership is exact.
        spread <- quartile[[3L]] - quartile[[1L]]
        lower <- quartile[[1L]] - iqr * spread
        upper <- quartile[[3L]] + iqr * spread
        local <- coefficients[index, level + 1L]
        outlier <- local < lower | local > upper
        rejected[[level]] <- rejected[[level]] + sum(outlier)
        local[outlier] <- 0
        coefficients[index, level + 1L] <- local
      }
    }
  }
  reconstructed <- .nirs_inverse_shift_invariant_db2(coefficients) + dc
  list(
    value = reconstructed[seq_len(original_length)],
    padded_length = padded_length,
    rejected_by_level = rejected
  )
}

#' Correct fNIRS motion with a shift-invariant db2 wavelet transform
#'
#' Outlying detail coefficients are removed using a per-block IQR rule and
#' the signal is reconstructed at its original length.
#'
#' @param x A governed `PhysioExperiment`.
#' @param assay_name Exact source assay name.
#' @param output_assay Exact new assay name.
#' @param iqr Positive IQR multiplier.
#' @param wavelet Exact wavelet name; currently only `"db2"`.
#' @param min_level Positive integer lowest analysis scale.
#'
#' @return A clone of `x` with one motion-corrected assay.
#' @export
waveletMotionCorrect <- function(
    x,
    assay_name = "OD",
    output_assay = paste0(assay_name, "_wavelet"),
    iqr = 1.5,
    wavelet = "db2",
    min_level = 4L) {
  context <- .nirs_motion_context(x, assay_name, min_samples = 2L)
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  iqr <- .nirs_motion_scalar(
    iqr, "iqr", lower = .Machine$double.xmin
  )
  wavelet <- .snirf_enum(wavelet, "db2", "wavelet")
  min_level <- .nirs_motion_scalar(
    min_level, "min_level", lower = 1,
    upper = .Machine$integer.max, integer = TRUE
  )
  corrected <- matrix(
    NA_real_,
    nrow = nrow(context$data),
    ncol = ncol(context$data),
    dimnames = dimnames(context$data)
  )
  diagnostics <- vector("list", ncol(context$data))
  for (j in seq_len(ncol(context$data))) {
    result <- .nirs_wavelet_channel(
      context$data[, j], iqr, min_level
    )
    corrected[, j] <- result$value
    diagnostics[[j]] <- list(
      measurement_id = context$identity$measurement_id[[j]],
      padded_length = result$padded_length,
      rejected_by_level = result$rejected_by_level
    )
  }
  .nirs_add_motion_assay(
    x,
    context,
    output_assay,
    corrected,
    "waveletMotionCorrect",
    params = list(
      iqr = iqr,
      wavelet = wavelet,
      min_level = min_level,
      transform = "Homer3 shift-invariant periodized db2",
      padding = "right zero padding to next power of two",
      quantile_type = 7L,
      diagnostics = diagnostics
    )
  )
}
