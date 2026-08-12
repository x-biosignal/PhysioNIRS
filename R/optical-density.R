.nirs_scalar_name <- function(x, arg, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be one non-empty character value", call. = FALSE)
  }
  x
}

.nirs_numeric_vector <- function(x, arg, positive = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
      anyNA(x) || any(!is.finite(x))) {
    stop("`", arg, "` must contain finite numeric values", call. = FALSE)
  }
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- original_names
  if (positive && any(x <= 0)) {
    stop("`", arg, "` must contain only positive values", call. = FALSE)
  }
  x
}

.nirs_cw_measurements <- function(x) {
  measurement <- measurementList(x)
  bad <- which(is.na(measurement$data_type) | measurement$data_type != 1L)
  if (length(bad)) {
    stop(
      "Continuous-wave amplitude data_type 1 is required; measurement(s) ",
      paste(bad, collapse = ", "), " have unsupported data_type",
      call. = FALSE
    )
  }
  wavelength <- as.numeric(measurement$wavelength_nm)
  if (anyNA(wavelength) || any(!is.finite(wavelength)) ||
      any(wavelength <= 0)) {
    stop("Continuous-wave measurements require positive wavelengths",
         call. = FALSE)
  }
  measurement
}

.nirs_assay_matrix <- function(x, assay_name) {
  assay_name <- .nirs_scalar_name(assay_name, "assay_name")
  available <- SummarizedExperiment::assayNames(x)
  if (!assay_name %in% available) {
    stop(
      "`assay_name` does not exist: ", assay_name,
      "; available assays: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  value <- SummarizedExperiment::assay(x, assay_name)
  if (!is.numeric(value) || length(dim(value)) != 2L) {
    stop("Assay `", assay_name, "` must be a numeric matrix", call. = FALSE)
  }
  value <- as.matrix(value)
  if (!identical(dim(value), dim(x))) {
    stop("Assay `", assay_name, "` dimensions disagree with `x`",
         call. = FALSE)
  }
  value
}

.nirs_resolve_intensity_assay <- function(x, assay_name) {
  if (!is.null(assay_name)) {
    return(.nirs_scalar_name(assay_name, "assay_name"))
  }
  available <- SummarizedExperiment::assayNames(x)
  contracts <- S4Vectors::metadata(x)$nirs$assays
  kinds <- vapply(
    available,
    function(name) {
      contract <- contracts[[name]]
      if (is.null(contract)) {
        "continuous_wave_intensity"
      } else {
        kind <- contract$kind
        if (!is.list(contract) || !is.character(kind) ||
            length(kind) != 1L || is.na(kind) || !nzchar(kind)) {
          stop(
            "Malformed assay contract for `", name,
            "`: `kind` must be one non-empty character value",
            call. = FALSE
          )
        }
        kind
      }
    },
    character(1)
  )
  eligible <- available[kinds == "continuous_wave_intensity"]
  if (length(eligible) != 1L) {
    stop(
      "`assay_name = NULL` requires exactly one eligible continuous-wave ",
      "intensity assay; found ", length(eligible), ": ",
      if (length(eligible)) paste(eligible, collapse = ", ") else "<none>",
      call. = FALSE
    )
  }
  eligible
}

.nirs_time_seconds <- function(x) {
  row_data <- SummarizedExperiment::rowData(x)
  if (!"time_seconds" %in% names(row_data)) {
    stop("`x` must carry `rowData(x)$time_seconds`", call. = FALSE)
  }
  time <- as.numeric(row_data$time_seconds)
  if (length(time) != nrow(x) || anyNA(time) || any(!is.finite(time)) ||
      (length(time) > 1L && any(diff(time) <= 0))) {
    stop("`rowData(x)$time_seconds` must be finite and strictly increasing",
         call. = FALSE)
  }
  time
}

.nirs_set_assay_contract <- function(x, name, contract) {
  metadata <- S4Vectors::metadata(x)
  if (is.null(metadata$nirs)) metadata$nirs <- list()
  if (is.null(metadata$nirs$assays)) metadata$nirs$assays <- list()
  metadata$nirs$assays[[name]] <- contract
  S4Vectors::metadata(x) <- metadata
  x
}

.nirs_append_step <- function(x, activity, params, input_assay, output_assay) {
  PhysioCore::appendProvenance(
    x,
    activity = activity,
    params = params,
    input_assay = input_assay,
    output_assay = output_assay,
    software_version = "0.4.0",
    package = "PhysioNIRS"
  )
}

#' Convert continuous-wave fNIRS intensity to optical density
#'
#' Optical density is calculated with the natural-log convention
#' `-log(I / I0)`. The reference intensity and any replacement of non-positive
#' cells are explicit and recorded in provenance.
#'
#' @param x A governed continuous-wave `PhysioExperiment`.
#' @param assay_name Exact input assay name. `NULL` is allowed only when one
#'   eligible intensity assay exists.
#' @param baseline Exactly `"mean"` or `"median"`.
#' @param baseline_interval Optional closed interval in seconds.
#' @param reference Optional positive scalar or vector with one value per
#'   measurement.
#' @param nonpositive Exactly `"error"` or `"replace"`.
#' @param replacement Explicit positive scalar used only with
#'   `nonpositive = "replace"`.
#' @param output_assay Exact name for the new optical-density assay.
#'
#' @return `x` with a unitless optical-density assay appended.
#' @export
intensityToOD <- function(
    x,
    assay_name = NULL,
    baseline = c("mean", "median"),
    baseline_interval = NULL,
    reference = NULL,
    nonpositive = c("error", "replace"),
    replacement = NULL,
    output_assay = "OD") {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  baseline <- if (missing(baseline)) {
    "mean"
  } else {
    .snirf_enum(baseline, c("mean", "median"), "baseline")
  }
  nonpositive <- if (missing(nonpositive)) {
    "error"
  } else {
    .snirf_enum(nonpositive, c("error", "replace"), "nonpositive")
  }
  assay_name <- .nirs_resolve_intensity_assay(x, assay_name)
  output_assay <- .nirs_scalar_name(output_assay, "output_assay")
  if (output_assay %in% SummarizedExperiment::assayNames(x)) {
    stop("`output_assay` already exists: ", output_assay, call. = FALSE)
  }
  measurement <- .nirs_cw_measurements(x)
  intensity <- .nirs_assay_matrix(x, assay_name)
  if (anyNA(intensity) || any(!is.finite(intensity))) {
    bad <- which(!is.finite(intensity), arr.ind = TRUE)[1L, ]
    stop(
      "Intensity must be finite; assay `", assay_name, "` row ",
      bad[[1L]], ", measurement ", bad[[2L]], " is non-finite",
      call. = FALSE
    )
  }

  nonpositive_index <- which(intensity <= 0, arr.ind = TRUE)
  replacement_count <- nrow(nonpositive_index)
  if (replacement_count && nonpositive == "error") {
    bad <- nonpositive_index[1L, ]
    stop(
      "Intensity must be positive; assay `", assay_name, "` row ",
      bad[[1L]], ", measurement ", bad[[2L]], " is non-positive",
      call. = FALSE
    )
  }
  if (nonpositive == "replace") {
    if (!is.numeric(replacement) || length(replacement) != 1L ||
        is.na(replacement) || !is.finite(replacement) || replacement <= 0) {
      stop(
        "`replacement` must be one explicit finite positive value when ",
        "`nonpositive = \"replace\"`",
        call. = FALSE
      )
    }
    if (replacement_count) {
      intensity[nonpositive_index] <- as.numeric(replacement)
      warning(
        "Replaced ", replacement_count, " non-positive intensity cell(s) ",
        "across ", length(unique(nonpositive_index[, 2L])),
        " measurement(s)",
        call. = FALSE
      )
    }
  } else if (!is.null(replacement)) {
    stop("`replacement` is only valid with `nonpositive = \"replace\"`",
         call. = FALSE)
  }

  if (!is.null(reference) && !is.null(baseline_interval)) {
    stop("`reference` and `baseline_interval` cannot both be supplied",
         call. = FALSE)
  }
  if (!is.null(reference)) {
    reference <- .nirs_numeric_vector(reference, "reference", positive = TRUE)
    if (!length(reference) %in% c(1L, ncol(intensity))) {
      stop(
        "`reference` length must be one or the number of measurements (",
        ncol(intensity), ")",
        call. = FALSE
      )
    }
    if (length(reference) == 1L) {
      reference <- rep(reference, ncol(intensity))
    }
    selected_rows <- integer()
    reference_source <- "explicit"
  } else {
    selected_rows <- seq_len(nrow(intensity))
    if (!is.null(baseline_interval)) {
      if (!is.numeric(baseline_interval) ||
          !is.null(dim(baseline_interval)) ||
          length(baseline_interval) != 2L ||
          anyNA(baseline_interval) ||
          any(!is.finite(baseline_interval)) ||
          baseline_interval[[1L]] > baseline_interval[[2L]]) {
        stop(
          "`baseline_interval` must be finite c(start_seconds, end_seconds) ",
          "with start <= end",
          call. = FALSE
        )
      }
      time <- .nirs_time_seconds(x)
      selected_rows <- which(
        time >= baseline_interval[[1L]] &
          time <= baseline_interval[[2L]]
      )
      if (!length(selected_rows)) {
        stop("`baseline_interval` selects no time rows", call. = FALSE)
      }
    }
    selected <- intensity[selected_rows, , drop = FALSE]
    reference <- if (baseline == "mean") {
      colMeans(selected)
    } else {
      apply(selected, 2L, stats::median)
    }
    reference_source <- baseline
  }
  if (anyNA(reference) || any(!is.finite(reference)) ||
      any(reference <= 0)) {
    bad <- which(!is.finite(reference) | reference <= 0)[[1L]]
    stop(
      "Derived reference must be finite and positive; measurement ",
      bad, " is invalid",
      call. = FALSE
    )
  }

  optical_density <- -sweep(
    log(intensity), 2L, log(reference), "-"
  )
  if (anyNA(optical_density) || any(!is.finite(optical_density))) {
    bad <- which(!is.finite(optical_density), arr.ind = TRUE)[1L, ]
    stop(
      "Optical-density calculation is non-finite at row ", bad[[1L]],
      ", measurement ", bad[[2L]],
      call. = FALSE
    )
  }
  dimnames(optical_density) <- dimnames(intensity)
  out <- x
  assays <- SummarizedExperiment::assays(out)
  assays[[output_assay]] <- optical_density
  SummarizedExperiment::assays(out) <- assays
  out <- .nirs_set_assay_contract(out, output_assay, list(
    kind = "optical_density",
    unit = "1",
    log_convention = "natural",
    source_assay = assay_name,
    measurement_count = nrow(measurement)
  ))
  .nirs_append_step(
    out,
    "intensityToOD",
    params = list(
      assay_name = assay_name,
      baseline = baseline,
      baseline_interval = if (is.null(baseline_interval)) {
        NULL
      } else {
        as.numeric(baseline_interval)
      },
      reference_source = reference_source,
      selected_rows = as.integer(selected_rows),
      reference = as.numeric(reference),
      nonpositive = nonpositive,
      replacement = if (is.null(replacement)) NULL else as.numeric(replacement),
      replacement_count = replacement_count,
      log_convention = "natural"
    ),
    input_assay = assay_name,
    output_assay = output_assay
  )
}
