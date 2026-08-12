#' Estimate a differential pathlength factor
#'
#' Implements published population equations. These values are population
#' estimates, not subject-specific pathlength measurements.
#'
#' @param wavelength_nm Positive wavelength values in nanometres.
#' @param age_years Non-negative age values in years.
#' @param model Exactly `"scholkmann2013"` or `"duncan1996"`.
#' @param extrapolate Whether to allow age/wavelength extrapolation where the
#'   selected model defines an unambiguous equation.
#'
#' @return A positive numeric vector with model and extrapolation attributes.
#' @export
dpf <- function(
    wavelength_nm,
    age_years,
    model = c("scholkmann2013", "duncan1996"),
    extrapolate = FALSE) {
  model <- if (missing(model)) {
    "scholkmann2013"
  } else {
    .snirf_enum(model, c("scholkmann2013", "duncan1996"), "model")
  }
  extrapolate <- .snirf_flag(extrapolate, "extrapolate")
  wavelength_nm <- .nirs_numeric_vector(
    wavelength_nm, "wavelength_nm", positive = TRUE
  )
  age_years <- .nirs_numeric_vector(age_years, "age_years")
  if (any(age_years < 0)) {
    stop("`age_years` must contain only non-negative values", call. = FALSE)
  }
  lengths <- c(length(wavelength_nm), length(age_years))
  if (lengths[[1L]] != lengths[[2L]] && all(lengths != 1L)) {
    stop(
      "`wavelength_nm` and `age_years` lengths must be one or equal",
      call. = FALSE
    )
  }
  n <- max(lengths)
  wavelength_nm <- rep(wavelength_nm, length.out = n)
  age_years <- rep(age_years, length.out = n)

  if (model == "scholkmann2013") {
    outside <- wavelength_nm < 650 | wavelength_nm > 950 |
      age_years > 70
    if (any(outside) && !extrapolate) {
      stop(
        "Scholkmann 2013 governed domain is wavelength 650-950 nm and ",
        "age 0-70 years; set `extrapolate = TRUE` explicitly",
        call. = FALSE
      )
    }
    value <- 223.3 + 0.05624 * age_years^0.8493 -
      5.723e-7 * wavelength_nm^3 +
      0.001245 * wavelength_nm^2 -
      0.9025 * wavelength_nm
  } else {
    supported <- c(690, 744, 807, 832)
    if (any(!wavelength_nm %in% supported)) {
      stop(
        "Duncan 1996 supports only exact wavelengths 690, 744, 807, and ",
        "832 nm; use `model = \"scholkmann2013\"` or explicit ",
        "`pathlength_factor`",
        call. = FALSE
      )
    }
    outside <- age_years > 50
    if (any(outside) && !extrapolate) {
      stop(
        "Duncan 1996 governed age domain is 0-50 years; set ",
        "`extrapolate = TRUE` explicitly",
        call. = FALSE
      )
    }
    value <- numeric(n)
    index <- wavelength_nm == 690
    value[index] <- 5.38 + 0.049 * age_years[index]^0.877
    index <- wavelength_nm == 744
    value[index] <- 5.11 + 0.106 * age_years[index]^0.723
    index <- wavelength_nm == 807
    value[index] <- 4.99 + 0.067 * age_years[index]^0.814
    index <- wavelength_nm == 832
    value[index] <- 4.67 + 0.062 * age_years[index]^0.819
  }
  if (any(outside)) {
    warning(
      "DPF equation extrapolated for ", sum(outside), " value(s)",
      call. = FALSE
    )
  }
  if (anyNA(value) || any(!is.finite(value)) || any(value <= 0)) {
    stop("DPF equation produced a non-finite or non-positive value",
         call. = FALSE)
  }
  value <- as.numeric(value)
  attr(value, "model") <- model
  attr(value, "extrapolated") <- as.logical(outside)
  value
}
