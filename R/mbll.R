.nirs_pair_contract <- function(measurement) {
  key <- paste(measurement$source_index, measurement$detector_index, sep = ":")
  pair_key <- unique(key)
  pair_index <- match(key, pair_key)
  groups <- lapply(seq_along(pair_key), function(i) which(pair_index == i))
  for (i in seq_along(groups)) {
    index <- groups[[i]]
    wavelengths <- as.numeric(measurement$wavelength_nm[index])
    if (anyDuplicated(wavelengths)) {
      stop(
        "Source-detector pair ", pair_key[[i]],
        " has duplicate wavelength measurements",
        call. = FALSE
      )
    }
    if (length(wavelengths) < 2L) {
      stop(
        "Source-detector pair ", pair_key[[i]],
        " requires at least two distinct wavelengths",
        call. = FALSE
      )
    }
  }
  list(key = key, pair_key = pair_key, pair_index = pair_index, groups = groups)
}

.nirs_resolve_by_measurement <- function(
    value, measurement, pair, arg, allow_wavelength = FALSE) {
  value <- .nirs_numeric_vector(value, arg, positive = TRUE)
  n_measurement <- nrow(measurement)
  n_pair <- length(pair$groups)
  if (length(value) == 1L) {
    return(rep(value, n_measurement))
  }
  wavelengths <- sort(unique(as.numeric(measurement$wavelength_nm)))
  if (allow_wavelength && !is.null(names(value))) {
    requested <- as.character(measurement$wavelength_nm)
    if (length(value) != length(wavelengths) ||
        anyNA(names(value)) || any(!nzchar(names(value))) ||
        anyDuplicated(names(value)) ||
        !setequal(names(value), as.character(wavelengths))) {
      stop(
        "Named `", arg,
        "` must name each unique measurement wavelength exactly",
        call. = FALSE
      )
    }
    return(unname(value[requested]))
  }
  if (length(value) == n_measurement) {
    return(value)
  }
  if (allow_wavelength && length(value) == length(wavelengths)) {
    return(value[match(measurement$wavelength_nm, wavelengths)])
  }
  if (!allow_wavelength && length(value) == n_pair) {
    return(value[pair$pair_index])
  }
  allowed <- c(1L, n_measurement)
  if (allow_wavelength) {
    allowed <- c(allowed, length(wavelengths))
  } else {
    allowed <- c(allowed, n_pair)
  }
  stop(
    "`", arg, "` length must be one of: ",
    paste(unique(allowed), collapse = ", "),
    call. = FALSE
  )
}

.nirs_pair_distances <- function(distance, pair) {
  result <- numeric(length(pair$groups))
  for (i in seq_along(pair$groups)) {
    values <- distance[pair$groups[[i]]]
    tolerance <- max(1e-12, 1e-8 * max(values))
    if (max(values) - min(values) > tolerance) {
      stop(
        "Source-detector pair ", pair$pair_key[[i]],
        " has wavelength-specific distance disagreement",
        call. = FALSE
      )
    }
    result[[i]] <- mean(values)
  }
  result
}

.nirs_mbll_coldata <- function(measurement, pair, distance, ppf) {
  source <- vapply(
    pair$groups,
    function(index) as.integer(measurement$source_index[index[[1L]]]),
    integer(1)
  )
  detector <- vapply(
    pair$groups,
    function(index) as.integer(measurement$detector_index[index[[1L]]]),
    integer(1)
  )
  wavelengths <- lapply(
    pair$groups,
    function(index) as.numeric(measurement$wavelength_nm[index])
  )
  factors <- lapply(pair$groups, function(index) as.numeric(ppf[index]))
  channel_id <- paste0("S", source, "_D", detector)
  if (anyDuplicated(channel_id)) {
    stop("Derived source-detector `channel_id` values are not unique",
         call. = FALSE)
  }
  S4Vectors::DataFrame(
    channel_id = channel_id,
    source_index = source,
    detector_index = detector,
    source_detector_distance_m = distance,
    wavelengths_nm = S4Vectors::List(wavelengths),
    pathlength_factors = S4Vectors::List(factors),
    concentration_unit = rep("uM", length(pair$groups)),
    row.names = channel_id
  )
}

#' Convert optical density to haemoglobin concentration change
#'
#' Solves the modified Beer-Lambert system independently for each exact
#' source-detector pair. The returned `PhysioExperiment` has one column per
#' pair and assays `HbO`, `HbR`, and `HbT` in micromolar.
#'
#' @param x A governed `PhysioExperiment` carrying natural-log optical density.
#' @param assay_name Exact optical-density assay name.
#' @param pathlength_factor Positive scalar/vector, or `NULL` to derive DPF from
#'   `age_years`.
#' @param age_years Age used only when `pathlength_factor = NULL`.
#' @param dpf_model Exact DPF model name passed to [dpf()].
#' @param distance_m Optional positive scalar, per-measurement vector, or
#'   per-pair vector overriding probe geometry.
#' @param output_unit Exactly `"uM"`.
#'
#' @return A pair-collapsed `PhysioExperiment` with concentration assays.
#' @export
mbll <- function(
    x,
    assay_name = "OD",
    pathlength_factor = 6,
    age_years = NULL,
    dpf_model = "scholkmann2013",
    distance_m = NULL,
    output_unit = "uM") {
  if (!inherits(x, "PhysioExperiment")) {
    stop("`x` must be a PhysioExperiment", call. = FALSE)
  }
  assay_name <- .nirs_scalar_name(assay_name, "assay_name")
  dpf_model <- .snirf_enum(
    dpf_model, c("scholkmann2013", "duncan1996"), "dpf_model"
  )
  output_unit <- .snirf_enum(output_unit, "uM", "output_unit")
  if (!is.null(pathlength_factor) && !is.null(age_years)) {
    stop(
      "`pathlength_factor` and `age_years` cannot both be supplied; set ",
      "`pathlength_factor = NULL` to derive DPF",
      call. = FALSE
    )
  }
  measurement <- .nirs_cw_measurements(x)
  optical_density <- .nirs_assay_matrix(x, assay_name)
  if (anyNA(optical_density) || any(!is.finite(optical_density))) {
    bad <- which(!is.finite(optical_density), arr.ind = TRUE)[1L, ]
    stop(
      "Optical density must be finite; row ", bad[[1L]],
      ", measurement ", bad[[2L]], " is non-finite",
      call. = FALSE
    )
  }
  assay_contract <- S4Vectors::metadata(x)$nirs$assays[[assay_name]]
  if (is.null(assay_contract) ||
      !identical(as.character(assay_contract$kind), "optical_density") ||
      !identical(
        as.character(assay_contract$log_convention), "natural"
      )) {
    stop(
      "`assay_name` must carry a governed natural-log optical-density ",
      "contract from `intensityToOD()`: ", assay_name,
      call. = FALSE
    )
  }

  pair <- .nirs_pair_contract(measurement)
  if (is.null(pathlength_factor)) {
    if (is.null(age_years)) {
      stop(
        "`age_years` is required when `pathlength_factor = NULL`",
        call. = FALSE
      )
    }
    ppf <- dpf(
      as.numeric(measurement$wavelength_nm),
      age_years,
      model = dpf_model
    )
    dpf_source <- dpf_model
    dpf_extrapolated <- attr(ppf, "extrapolated")
    ppf <- as.numeric(ppf)
  } else {
    ppf <- .nirs_resolve_by_measurement(
      pathlength_factor, measurement, pair, "pathlength_factor",
      allow_wavelength = TRUE
    )
    dpf_source <- "explicit"
    dpf_extrapolated <- rep(FALSE, length(ppf))
  }

  if (is.null(distance_m)) {
    distance <- as.numeric(
      sourceDetectorDistances(x, unit = "m")$distance
    )
    distance <- .nirs_numeric_vector(
      distance, "sourceDetectorDistances(x)", positive = TRUE
    )
    distance_source <- "probe_geometry"
  } else {
    distance <- .nirs_resolve_by_measurement(
      distance_m, measurement, pair, "distance_m"
    )
    distance_source <- "explicit_override"
    warning(
      "`distance_m` overrides governed probe geometry",
      call. = FALSE
    )
  }
  pair_distance <- .nirs_pair_distances(distance, pair)

  n_time <- nrow(optical_density)
  n_pair <- length(pair$groups)
  hbo <- matrix(NA_real_, nrow = n_time, ncol = n_pair)
  hbr <- matrix(NA_real_, nrow = n_time, ncol = n_pair)
  condition_number <- numeric(n_pair)
  mapping <- vector("list", n_pair)
  for (i in seq_along(pair$groups)) {
    index <- pair$groups[[i]]
    order <- order(as.numeric(measurement$wavelength_nm[index]))
    algebra_index <- index[order]
    wavelength <- as.numeric(measurement$wavelength_nm[algebra_index])
    extinction <- extinctionCoefficients(wavelength, unit = "m-1 M-1")
    design <- extinction * (pair_distance[[i]] * ppf[algebra_index])
    if (anyNA(design) || any(!is.finite(design))) {
      stop(
        "Non-finite MBLL design for source-detector pair ",
        pair$pair_key[[i]],
        call. = FALSE
      )
    }
    singular <- svd(design, nu = 0L, nv = 0L)$d
    tolerance <- max(dim(design)) * max(singular) * .Machine$double.eps
    rank <- sum(singular > tolerance)
    condition_number[[i]] <- max(singular) / min(singular)
    if (rank < 2L || !is.finite(condition_number[[i]]) ||
        condition_number[[i]] > 1e8) {
      stop(
        "MBLL design is rank-deficient or ill-conditioned for ",
        "source-detector pair ", pair$pair_key[[i]],
        " (rank=", rank, ", condition=",
        format(condition_number[[i]], digits = 6L), ")",
        call. = FALSE
      )
    }
    concentration <- t(qr.solve(
      design,
      t(optical_density[, algebra_index, drop = FALSE]),
      tol = tolerance
    )) * 1e6
    if (anyNA(concentration) || any(!is.finite(concentration))) {
      stop(
        "MBLL solve produced non-finite concentration for ",
        "source-detector pair ", pair$pair_key[[i]],
        call. = FALSE
      )
    }
    hbo[, i] <- concentration[, 1L]
    hbr[, i] <- concentration[, 2L]
    mapping[[i]] <- list(
      output_pair_index = i,
      channel_id = paste0(
        "S", measurement$source_index[index[[1L]]],
        "_D", measurement$detector_index[index[[1L]]]
      ),
      input_measurement_indices = as.integer(index),
      algebra_measurement_indices = as.integer(algebra_index)
    )
  }
  column_names <- vapply(mapping, `[[`, character(1), "channel_id")
  dimnames(hbo) <- list(rownames(optical_density), column_names)
  dimnames(hbr) <- list(rownames(optical_density), column_names)
  hbt <- hbo + hbr
  if (anyNA(hbt) || any(!is.finite(hbt))) {
    stop("HbT calculation produced non-finite concentration",
         call. = FALSE)
  }
  pair_data <- .nirs_mbll_coldata(
    measurement, pair, pair_distance, ppf
  )
  metadata <- S4Vectors::metadata(x)
  if (is.null(metadata$nirs)) metadata$nirs <- list()
  source_assay_contract <- metadata$nirs$assays[[assay_name]]
  extinction_meta <- .nirs_extinction_metadata()
  metadata$nirs$mbll <- c(
    extinction_meta,
    list(
      source_assay = assay_name,
      source_assay_contract = source_assay_contract,
      log_convention = "natural",
      dpf_source = dpf_source,
      dpf_model = if (dpf_source == "explicit") NULL else dpf_model,
      dpf_extrapolated = as.logical(dpf_extrapolated),
      distance_source = distance_source,
      output_unit = output_unit,
      condition_numbers = condition_number,
      input_to_pair_mapping = mapping
    )
  )
  metadata$nirs$assays <- list(
    HbO = list(
      kind = "haemoglobin_concentration",
      chromophore = "HbO",
      unit = "uM",
      source_assay = assay_name,
      method = "mbll"
    ),
    HbR = list(
      kind = "haemoglobin_concentration",
      chromophore = "HbR",
      unit = "uM",
      source_assay = assay_name,
      method = "mbll"
    ),
    HbT = list(
      kind = "haemoglobin_concentration",
      chromophore = "HbT",
      unit = "uM",
      source_assay = c("HbO", "HbR"),
      method = "sum"
    )
  )
  out <- PhysioCore::PhysioExperiment(
    assays = list(HbO = hbo, HbR = hbr, HbT = hbt),
    rowData = SummarizedExperiment::rowData(x),
    colData = pair_data,
    metadata = metadata,
    samplingRate = PhysioCore::samplingRate(x)
  )
  methods::validObject(out)
  .nirs_append_step(
    out,
    "mbll",
    params = list(
      assay_name = assay_name,
      pathlength_factor = if (dpf_source == "explicit") {
        as.numeric(ppf)
      } else {
        NULL
      },
      age_years = if (is.null(age_years)) NULL else as.numeric(age_years),
      dpf_model = if (dpf_source == "explicit") NULL else dpf_model,
      distance_source = distance_source,
      distance_m = as.numeric(pair_distance),
      output_unit = output_unit,
      condition_numbers = condition_number,
      coefficient_sha256 = extinction_meta$table_sha256,
      log_convention = "natural"
    ),
    input_assay = assay_name,
    output_assay = "HbO,HbR,HbT"
  )
}
