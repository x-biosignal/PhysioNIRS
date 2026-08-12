.nirs_design_bands <- function(physiology_bands, band_ranges, fs) {
  if (!is.character(physiology_bands) || !is.null(dim(physiology_bands)) ||
      anyNA(physiology_bands) || any(!nzchar(physiology_bands)) ||
      anyDuplicated(physiology_bands) ||
      any(!physiology_bands %in% names(.nirs_physiology_ranges))) {
    stop(
      "`physiology_bands` must contain unique exact physiology-band names",
      call. = FALSE
    )
  }
  if (is.null(band_ranges)) {
    band_ranges <- stats::setNames(
      rep(list(NULL), length(physiology_bands)),
      physiology_bands
    )
  } else {
    if (!is.list(band_ranges) || is.null(names(band_ranges)) ||
        anyNA(names(band_ranges)) || any(!nzchar(names(band_ranges))) ||
        anyDuplicated(names(band_ranges)) ||
        !identical(sort(names(band_ranges)), sort(physiology_bands))) {
      stop(
        "`band_ranges` must name every requested band exactly with no extras",
        call. = FALSE
      )
    }
  }
  lapply(physiology_bands, function(band) {
    .nirs_band_range(band, band_ranges[[band]], fs)
  })
}

.nirs_short_design_base <- function(context, map, aggregation) {
  short <- which(map$is_short)
  if (!length(short) || all(map$is_short)) {
    stop(
      "A nuisance design requires at least one short and one long channel",
      call. = FALSE
    )
  }
  if (aggregation == "all") {
    value <- context$data[, short, drop = FALSE]
    names <- paste0("short_", context$identity$channel_id[short])
    roles <- data.frame(
      column = names,
      role = rep("short_channel", length(short)),
      source_indices = vapply(
        short, function(i) as.character(i), character(1)
      ),
      compatibility = if (
        context$identity$identity_kind[[1L]] == "measurement"
      ) {
        paste0(
          "wavelength_", .nirs_short_number_id(
            context$identity$wavelength_nm[short]
          )
        )
      } else {
        rep(context$assay_name, length(short))
      },
      stringsAsFactors = FALSE
    )
  } else if (context$identity$identity_kind[[1L]] == "measurement") {
    wavelength <- context$identity$wavelength_nm[short]
    levels <- unique(wavelength)
    groups <- lapply(levels, function(value) short[wavelength == value])
    value <- vapply(
      groups,
      function(index) rowMeans(context$data[, index, drop = FALSE]),
      numeric(nrow(context$data))
    )
    if (!is.matrix(value)) {
      value <- matrix(value, ncol = length(groups))
    }
    names <- paste0(
      "short_mean_wl", .nirs_short_number_id(levels)
    )
    roles <- data.frame(
      column = names,
      role = rep("short_mean", length(groups)),
      source_indices = vapply(
        groups, function(index) paste(index, collapse = ","), character(1)
      ),
      compatibility = paste0(
        "wavelength_", .nirs_short_number_id(levels)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    value <- matrix(
      rowMeans(context$data[, short, drop = FALSE]),
      ncol = 1L
    )
    names <- paste0("short_mean_", context$assay_name)
    roles <- data.frame(
      column = names,
      role = "short_mean",
      source_indices = paste(short, collapse = ","),
      compatibility = context$assay_name,
      stringsAsFactors = FALSE
    )
  }
  colnames(value) <- names
  rownames(value) <- rownames(context$data)
  list(value = value, roles = roles)
}

#' Construct a governed short-separation nuisance design
#'
#' Returns short-channel signals, optionally aggregated by compatibility class,
#' with requested physiology-band copies. No intercept, task, or drift
#' regressor is added.
#'
#' @param x A governed OD or pair-collapsed haemoglobin
#'   `PhysioExperiment`.
#' @param assay_name Exact source assay.
#' @param short A compatible `nirs_short_channels` map returned by
#'   `identifyShortChannels()`, or `NULL` for the strict 10-mm default.
#' @param aggregation Exactly `"all"` or `"mean"`.
#' @param physiology_bands Unique exact band names to append.
#' @param band_ranges Optional named list providing one explicit range per
#'   requested band.
#' @param standardize Whether to centre and sample-SD standardize every column.
#'
#' @return A time-aligned `nirs_nuisance_design` numeric matrix.
#' @export
shortSeparationDesign <- function(
    x,
    assay_name = "OD",
    short = NULL,
    aggregation = c("all", "mean"),
    physiology_bands = character(),
    band_ranges = NULL,
    standardize = TRUE) {
  aggregation <- if (missing(aggregation)) {
    "all"
  } else {
    .snirf_enum(aggregation, c("all", "mean"), "aggregation")
  }
  standardize <- .snirf_flag(standardize, "standardize")
  context <- .nirs_short_context(
    x, assay_name, min_samples = 3L, require_contract = TRUE
  )
  map <- .nirs_resolve_short_map(x, short, context)
  resolved_bands <- .nirs_design_bands(
    physiology_bands, band_ranges, context$fs
  )
  base <- .nirs_short_design_base(context, map, aggregation)
  design <- base$value
  roles <- base$roles
  roles$band <- rep(NA_character_, nrow(roles))
  roles$range_low_hz <- rep(NA_real_, nrow(roles))
  roles$range_high_hz <- rep(NA_real_, nrow(roles))

  for (i in seq_along(resolved_bands)) {
    resolved <- resolved_bands[[i]]
    filter_spec <- .nirs_butter_band(
      resolved$range_hz, context$fs, order = 3L
    )
    filtered <- .nirs_zero_phase_band(base$value, filter_spec)
    names <- paste0(resolved$band, "_", colnames(base$value))
    colnames(filtered) <- names
    design <- cbind(design, filtered)
    roles <- rbind(
      roles,
      data.frame(
        column = names,
        role = rep("physiology_band", ncol(filtered)),
        source_indices = base$roles$source_indices,
        compatibility = base$roles$compatibility,
        band = rep(resolved$band, ncol(filtered)),
        range_low_hz = rep(resolved$range_hz[[1L]], ncol(filtered)),
        range_high_hz = rep(resolved$range_hz[[2L]], ncol(filtered)),
        stringsAsFactors = FALSE
      )
    )
    if (!resolved$explicit) .nirs_warn_default_physiology_range()
  }
  if (anyDuplicated(colnames(design))) {
    stop("Nuisance-design column names are not unique", call. = FALSE)
  }
  centre <- rep(0, ncol(design))
  scale <- rep(1, ncol(design))
  if (standardize) {
    centre <- colMeans(design)
    scale <- apply(design, 2L, stats::sd)
    near_zero <- scale <= sqrt(.Machine$double.eps) *
      pmax(1, apply(abs(design), 2L, max))
    if (anyNA(scale) || any(!is.finite(scale)) || any(near_zero)) {
      stop(
        "Nuisance design contains a zero or near-zero variance column: ",
        paste(colnames(design)[which(near_zero)], collapse = ", "),
        call. = FALSE
      )
    }
    design <- sweep(design, 2L, centre, "-")
    design <- sweep(design, 2L, scale, "/")
  }
  if (anyNA(design) || any(!is.finite(design))) {
    stop("Nuisance design produced non-finite values", call. = FALSE)
  }
  names(centre) <- colnames(design)
  names(scale) <- colnames(design)
  rownames(design) <- rownames(context$data)
  attr(design, "time_seconds") <- context$time
  attr(design, "sampling_rate_hz") <- context$fs
  attr(design, "column_roles") <- roles
  attr(design, "bands") <- resolved_bands
  attr(design, "standardization") <- list(
    applied = standardize,
    center = centre,
    scale = scale,
    sample_sd = TRUE
  )
  attr(design, "short_map_fingerprint") <- .nirs_sha256(map)
  attr(design, "source_fingerprint") <- context$source_fingerprint
  attr(design, "implementation_version") <- "0.4.0"
  class(design) <- c("nirs_nuisance_design", "matrix")
  design
}
