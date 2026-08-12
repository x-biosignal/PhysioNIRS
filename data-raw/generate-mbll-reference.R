# Generate deterministic, package-owned MBLL fixtures without calling
# PhysioNIRS production functions.

root <- normalizePath("physio-ecosystem/PhysioNIRS")
table_path <- file.path(
  root, "inst", "extdata", "haemoglobin-extinction.csv"
)
output <- file.path(root, "inst", "extdata", "mbll_reference.rds")
generator_path <- file.path(root, "data-raw", "generate-mbll-reference.R")
table <- utils::read.csv(table_path, stringsAsFactors = FALSE)

lookup <- function(wavelength_nm) {
  cbind(
    HbO = stats::approx(
      table$wavelength_nm, table$hbo2_cm1_M1, wavelength_nm,
      method = "linear", ties = "ordered"
    )$y,
    HbR = stats::approx(
      table$wavelength_nm, table$hbr_cm1_M1, wavelength_nm,
      method = "linear", ties = "ordered"
    )$y
  ) * (100 * log(10))
}

scholkmann <- function(wavelength_nm, age_years) {
  223.3 + 0.05624 * age_years^0.8493 -
    5.723e-7 * wavelength_nm^3 +
    0.001245 * wavelength_nm^2 -
    0.9025 * wavelength_nm
}

duncan <- function(wavelength_nm, age_years) {
  constants <- list(
    `690` = c(5.38, 0.049, 0.877),
    `744` = c(5.11, 0.106, 0.723),
    `807` = c(4.99, 0.067, 0.814),
    `832` = c(4.67, 0.062, 0.819)
  )
  vapply(seq_along(wavelength_nm), function(i) {
    value <- constants[[as.character(wavelength_nm[[i]])]]
    value[[1L]] + value[[2L]] * age_years^value[[3L]]
  }, numeric(1))
}

make_case <- function(
    id, wavelengths, distance_m, ppf, dpf_model, age_years = NULL) {
  time <- seq.int(0, 300) / 10
  concentration_uM <- cbind(
    HbO = 3.2 * sin(2 * pi * 0.05 * time) +
      0.6 * (time >= 10 & time <= 20),
    HbR = -1.1 * sin(2 * pi * 0.05 * time + 0.15) -
      0.25 * (time >= 10 & time <= 20)
  )
  coefficient <- lookup(wavelengths)
  optical_density <- (concentration_uM * 1e-6) %*%
    t(coefficient * (distance_m * ppf))
  reference_intensity <- seq(900, 1200, length.out = length(wavelengths))
  intensity <- sweep(exp(-optical_density), 2L, reference_intensity, "*")
  list(
    id = id,
    time_seconds = time,
    wavelengths_nm = wavelengths,
    distance_m = distance_m,
    pathlength_factor = ppf,
    dpf_model = dpf_model,
    age_years = age_years,
    reference_intensity = reference_intensity,
    intensity = intensity,
    optical_density = optical_density,
    HbO_uM = concentration_uM[, "HbO"],
    HbR_uM = concentration_uM[, "HbR"],
    HbT_uM = rowSums(concentration_uM)
  )
}

cases <- list(
  make_case(
    "two_explicit", c(760, 850), 0.02, c(5.8, 6.2), "explicit"
  ),
  make_case(
    "three_scholkmann", c(690, 760, 850), 0.03,
    scholkmann(c(690, 760, 850), 30),
    "scholkmann2013", 30
  ),
  make_case(
    "four_explicit", c(690, 760, 830, 850), 0.045,
    c(5.5, 5.8, 6.1, 6.4), "explicit"
  ),
  make_case(
    "three_duncan", c(744, 807, 832), 0.035,
    duncan(c(744, 807, 832), 20),
    "duncan1996", 20
  )
)

fixture <- list(
  schema_version = "1.0.0",
  generated_by = basename(generator_path),
  generator_sha256 = digest::digest(
    generator_path, algo = "sha256", file = TRUE
  ),
  extinction_table_sha256 = digest::digest(
    table_path, algo = "sha256", file = TRUE
  ),
  coefficient_convention = list(
    stored = "base10 cm-1 M-1",
    solver = "natural m-1 M-1",
    conversion = "100 * ln(10)",
    concentration = "uM"
  ),
  rng_seed = NULL,
  cases = cases
)

saveRDS(fixture, output, version = 3, compress = "xz")
message("fixture_sha256=",
        digest::digest(output, algo = "sha256", file = TRUE))
