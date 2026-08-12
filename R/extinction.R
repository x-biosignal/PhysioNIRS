.nirs_extinction_cache <- new.env(parent = emptyenv())

.nirs_extinction_table <- function() {
  if (exists("table", envir = .nirs_extinction_cache, inherits = FALSE)) {
    return(get("table", envir = .nirs_extinction_cache, inherits = FALSE))
  }
  path <- system.file(
    "extdata", "haemoglobin-extinction.csv", package = "PhysioNIRS"
  )
  if (!nzchar(path) || !file.exists(path)) {
    stop("Packaged haemoglobin extinction table is unavailable",
         call. = FALSE)
  }
  table <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    colClasses = c("numeric", "numeric", "numeric", "character", "character"),
    check.names = FALSE
  )
  required <- c(
    "wavelength_nm", "hbo2_cm1_M1", "hbr_cm1_M1",
    "source", "source_version"
  )
  if (!identical(names(table), required) || !nrow(table) ||
      anyNA(table) || any(!is.finite(table$wavelength_nm)) ||
      any(table$wavelength_nm <= 0) ||
      any(diff(table$wavelength_nm) <= 0) ||
      any(!is.finite(table$hbo2_cm1_M1)) ||
      any(!is.finite(table$hbr_cm1_M1)) ||
      any(table$hbo2_cm1_M1 < 0) || any(table$hbr_cm1_M1 < 0) ||
      any(!nzchar(table$source)) || any(!nzchar(table$source_version)) ||
      length(unique(table$source)) != 1L ||
      length(unique(table$source_version)) != 1L) {
    stop("Packaged haemoglobin extinction table violates its schema",
         call. = FALSE)
  }
  assign("table", table, envir = .nirs_extinction_cache)
  table
}

.nirs_extinction_metadata <- function(table = .nirs_extinction_table()) {
  hash_path <- system.file(
    "extdata", "haemoglobin-extinction.sha256", package = "PhysioNIRS"
  )
  hash <- if (nzchar(hash_path) && file.exists(hash_path)) {
    strsplit(readLines(hash_path, warn = FALSE)[[1L]], "[[:space:]]+")[[1L]][[1L]]
  } else {
    NA_character_
  }
  if (!is.character(hash) || length(hash) != 1L ||
      !grepl("^[0-9a-f]{64}$", hash)) {
    stop("Packaged extinction SHA-256 manifest is malformed", call. = FALSE)
  }
  list(
    source = unique(table$source),
    source_version = unique(table$source_version),
    table_sha256 = hash,
    source_unit = "cm-1 M-1",
    source_log_convention = "base10",
    natural_unit_conversion = "100 * ln(10)",
    wavelength_domain_nm = range(table$wavelength_nm)
  )
}

.nirs_linear_lookup <- function(x, grid, value) {
  result <- numeric(length(x))
  exact <- match(x, grid)
  is_exact <- !is.na(exact)
  result[is_exact] <- value[exact[is_exact]]
  pending <- which(!is_exact)
  for (i in pending) {
    if (x[[i]] < grid[[1L]]) {
      left <- 1L
      right <- 2L
    } else if (x[[i]] > grid[[length(grid)]]) {
      right <- length(grid)
      left <- right - 1L
    } else {
      right <- which(grid > x[[i]])[[1L]]
      left <- right - 1L
    }
    fraction <- (x[[i]] - grid[[left]]) / (grid[[right]] - grid[[left]])
    result[[i]] <- value[[left]] +
      fraction * (value[[right]] - value[[left]])
  }
  result
}

#' Get governed haemoglobin extinction coefficients
#'
#' Coefficients are linearly interpolated from the package-owned table. Stored
#' coefficients use the base-10 molar convention in `cm^-1 M^-1`; natural-log
#' coefficients in `m^-1 M^-1` include the explicit `100 * log(10)` conversion.
#'
#' @param wavelength_nm Positive wavelengths in nanometres.
#' @param unit Exactly `"m-1 M-1"` or `"cm-1 M-1"`.
#' @param extrapolate Whether to permit terminal-slope linear extrapolation.
#'
#' @return A numeric matrix with exact columns `HbO` and `HbR`.
#' @export
extinctionCoefficients <- function(
    wavelength_nm,
    unit = c("m-1 M-1", "cm-1 M-1"),
    extrapolate = FALSE) {
  unit <- if (missing(unit)) {
    "m-1 M-1"
  } else {
    .snirf_enum(unit, c("m-1 M-1", "cm-1 M-1"), "unit")
  }
  extrapolate <- .snirf_flag(extrapolate, "extrapolate")
  wavelength_nm <- .nirs_numeric_vector(
    wavelength_nm, "wavelength_nm", positive = TRUE
  )
  table <- .nirs_extinction_table()
  domain <- range(table$wavelength_nm)
  outside <- wavelength_nm < domain[[1L]] | wavelength_nm > domain[[2L]]
  if (any(outside) && !extrapolate) {
    stop(
      "`wavelength_nm` is outside the extinction table domain ",
      domain[[1L]], "-", domain[[2L]],
      " nm; set `extrapolate = TRUE` explicitly",
      call. = FALSE
    )
  }
  if (any(outside)) {
    warning(
      "Extinction coefficients extrapolated for ", sum(outside),
      " wavelength(s)",
      call. = FALSE
    )
  }
  result <- cbind(
    HbO = .nirs_linear_lookup(
      wavelength_nm, table$wavelength_nm, table$hbo2_cm1_M1
    ),
    HbR = .nirs_linear_lookup(
      wavelength_nm, table$wavelength_nm, table$hbr_cm1_M1
    )
  )
  if (unit == "m-1 M-1") {
    result <- result * (100 * log(10))
  }
  if (anyNA(result) || any(!is.finite(result)) || any(result < 0)) {
    stop(
      "Extinction interpolation produced a non-finite or negative ",
      "coefficient",
      call. = FALSE
    )
  }
  rownames(result) <- NULL
  attr(result, "wavelength_nm") <- as.numeric(wavelength_nm)
  attr(result, "unit") <- unit
  attr(result, "extrapolated") <- as.logical(outside)
  result
}
