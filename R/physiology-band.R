.nirs_physiology_ranges <- list(
  mayer = c(0.07, 0.13),
  respiration = c(0.15, 0.40),
  cardiac = c(0.70, 2.00)
)

.nirs_physiology_state <- new.env(parent = emptyenv())
.nirs_physiology_state$default_warning_emitted <- FALSE

.nirs_band_range <- function(band, range_hz, fs) {
  band <- .snirf_enum(band, names(.nirs_physiology_ranges), "band")
  explicit <- !is.null(range_hz)
  range <- if (explicit) {
    .nirs_numeric_vector(range_hz, "range_hz", positive = TRUE)
  } else {
    .nirs_physiology_ranges[[band]]
  }
  if (length(range) != 2L || range[[1L]] >= range[[2L]]) {
    stop("`range_hz` must be a strictly increasing length-two vector",
         call. = FALSE)
  }
  if (range[[2L]] >= fs / 2) {
    stop("Both physiology-band edges must be strictly below Nyquist",
         call. = FALSE)
  }
  list(band = band, range_hz = as.numeric(range), explicit = explicit)
}

.nirs_butter_band <- function(range_hz, fs, order) {
  order <- .nirs_motion_scalar(
    order, "order", lower = 1, upper = .Machine$integer.max,
    integer = TRUE
  )
  filter <- signal::butter(
    n = order,
    W = 2 * range_hz / fs,
    type = "pass"
  )
  minimum <- 2L * max(length(filter$a), length(filter$b)) + 1L
  list(filter = filter, order = order, minimum_samples = minimum)
}

.nirs_zero_phase_band <- function(values, filter_spec) {
  if (nrow(values) < filter_spec$minimum_samples) {
    stop(
      "Physiology bandpass requires at least ",
      filter_spec$minimum_samples, " samples for these coefficients",
      call. = FALSE
    )
  }
  out <- matrix(
    0,
    nrow = nrow(values),
    ncol = ncol(values),
    dimnames = dimnames(values)
  )
  for (j in seq_len(ncol(values))) {
    source <- values[, j]
    if (all(source == source[[1L]])) {
      out[, j] <- 0
    } else {
      out[, j] <- as.numeric(signal::filtfilt(filter_spec$filter, source))
    }
  }
  if (anyNA(out) || any(!is.finite(out))) {
    stop("Physiology bandpass produced non-finite values", call. = FALSE)
  }
  out
}

.nirs_warn_default_physiology_range <- function() {
  if (!isTRUE(.nirs_physiology_state$default_warning_emitted)) {
    .nirs_physiology_state$default_warning_emitted <- TRUE
    warning(
      "Default adult physiology ranges can be inappropriate for paediatric, ",
      "exercise, autonomic, or impaired populations; supply `range_hz` when ",
      "a governed population-specific range is required",
      call. = FALSE
    )
  }
}

#' Extract an explicit fNIRS physiology frequency band
#'
#' Applies a zero-phase Butterworth band-pass independently to every source
#' assay column. Band names identify common nuisance ranges; they are not
#' physiological diagnoses.
#'
#' @details Default passbands in Hz are Mayer `0.07--0.13`, respiration
#'   `0.15--0.40`, and cardiac `0.70--2.00`. Filtering uses an order-`order`
#'   Butterworth band-pass through `signal::filtfilt()`. Explicit ranges must
#'   remain strictly below the governed Nyquist frequency.
#'
#' @param x A governed `PhysioExperiment`.
#' @param assay_name Exact source assay name.
#' @param band Exactly `"mayer"`, `"respiration"`, or `"cardiac"`.
#' @param range_hz Optional explicit positive passband below Nyquist.
#' @param output_assay Exact name for the added filtered assay.
#' @param order Positive integer Butterworth order.
#'
#' @return A clone of `x` with one zero-phase filtered assay.
#' @export
physiologyBandpass <- function(
    x,
    assay_name = "OD",
    band = c("mayer", "respiration", "cardiac"),
    range_hz = NULL,
    output_assay = paste0(assay_name, "_", band[[1L]]),
    order = 3L) {
  band <- if (missing(band)) {
    "mayer"
  } else {
    .snirf_enum(band, names(.nirs_physiology_ranges), "band")
  }
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  context <- .nirs_short_context(
    x, assay_name, min_samples = 3L, require_contract = TRUE
  )
  if (output_assay %in% SummarizedExperiment::assayNames(x)) {
    stop("`output_assay` already exists: ", output_assay, call. = FALSE)
  }
  resolved <- .nirs_band_range(band, range_hz, context$fs)
  filter_spec <- .nirs_butter_band(
    resolved$range_hz, context$fs, order
  )
  filtered <- .nirs_zero_phase_band(context$data, filter_spec)
  if (!resolved$explicit) .nirs_warn_default_physiology_range()

  out <- x
  SummarizedExperiment::assay(out, output_assay, withDimnames = FALSE) <-
    filtered
  output_contract <- context$contract
  output_contract$source_assay <- context$assay_name
  output_contract$method <- "physiology_bandpass"
  output_contract$band <- resolved$band
  output_contract$range_hz <- resolved$range_hz
  out <- .nirs_set_assay_contract(out, output_assay, output_contract)
  .nirs_append_step(
    out,
    "physiologyBandpass",
    params = list(
      implementation_version = "0.4.0",
      source_fingerprint = context$source_fingerprint,
      band = resolved$band,
      range_hz = resolved$range_hz,
      explicit_range = resolved$explicit,
      order = filter_spec$order,
      sampling_rate_hz = context$fs,
      filter_type = "Butterworth band-pass, signal::filtfilt",
      filter_b = as.numeric(filter_spec$filter$b),
      filter_a = as.numeric(filter_spec$filter$a),
      minimum_samples = filter_spec$minimum_samples,
      edge_handling = "signal::filtfilt zero extension"
    ),
    input_assay = context$assay_name,
    output_assay = output_assay
  )
}
