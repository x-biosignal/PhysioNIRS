#' Calculate source-detector distances for each SNIRF measurement
#'
#' @param x A governed `PhysioExperiment`.
#' @param dimension Geometry selection: exactly `"auto"`, `"3d"`, or `"2d"`.
#' @param unit Output unit: exactly `"m"` or `"native"`.
#'
#' @return A `DataFrame` with one row per measurement and distance identity.
#' @export
sourceDetectorDistances <- function(
    x,
    dimension = c("auto", "3d", "2d"),
    unit = c("m", "native")) {
  dimension <- if (missing(dimension)) {
    "auto"
  } else {
    .snirf_enum(dimension, c("auto", "3d", "2d"), "dimension")
  }
  unit <- if (missing(unit)) {
    "m"
  } else {
    .snirf_enum(unit, c("m", "native"), "unit")
  }
  measurement <- measurementList(x)
  snirf <- S4Vectors::metadata(x)$snirf
  probe <- .snirf_validate_probe(snirf$probe, measurement)
  use <- dimension
  if (use == "auto") {
    use <- if (all(c("sourcePos3D", "detectorPos3D") %in% names(probe))) {
      "3d"
    } else {
      "2d"
    }
  }
  source_field <- if (use == "3d") "sourcePos3D" else "sourcePos2D"
  detector_field <- if (use == "3d") "detectorPos3D" else "detectorPos2D"
  if (!all(c(source_field, detector_field) %in% names(probe))) {
    stop("Requested probe geometry is not available: ", use, call. = FALSE)
  }
  source <- probe[[source_field]][measurement$source_index, , drop = FALSE]
  detector <- probe[[detector_field]][measurement$detector_index, ,
                                      drop = FALSE]
  distance <- sqrt(rowSums((source - detector)^2))
  native_unit <- probe$LengthUnit
  output_unit <- native_unit
  if (unit == "m") {
    distance <- distance * .snirf_length_factor(native_unit)
    output_unit <- "m"
  } else if (!is.character(native_unit) || length(native_unit) != 1L ||
             is.na(native_unit) || !nzchar(native_unit)) {
    stop("A declared SNIRF LengthUnit is required", call. = FALSE)
  }
  out <- measurement[, c(
    "measurement_index", "source_index", "detector_index",
    "wavelength_index", "wavelength_nm", "channel_label"
  )]
  out$distance <- as.numeric(distance)
  out$distance_unit <- rep(output_unit, nrow(out))
  out$geometry_dimension <- rep(use, nrow(out))
  out
}
