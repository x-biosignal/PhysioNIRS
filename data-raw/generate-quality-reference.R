#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- if (length(file_arg)) {
  normalizePath(file_arg[[1L]], mustWork = TRUE)
} else {
  normalizePath(
    sys.frame(1)$ofile %||% "data-raw/generate-quality-reference.R",
    mustWork = TRUE
  )
}
package_root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
output_dir <- file.path(package_root, "inst", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

curated_fs <- 20
curated_time <- seq.int(0, 400L - 1L) / curated_fs
curated <- lapply(c(2L, 3L, 4L), function(n_wavelength) {
  coupled <- vapply(seq_len(n_wavelength), function(j) {
    sin(2 * pi * curated_time + 0.02 * j) +
      0.01 * sin(2 * pi * (2.1 + j / 20) * curated_time)
  }, numeric(length(curated_time)))
  decoupled <- vapply(seq_len(n_wavelength), function(j) {
    sin(2 * pi * (0.15 + 0.11 * j) * curated_time + j)
  }, numeric(length(curated_time)))
  list(
    wavelength_nm = c(730, 760, 810, 850)[seq_len(n_wavelength)],
    coupled = coupled,
    decoupled = decoupled
  )
})
names(curated) <- c("two_wavelength", "three_wavelength", "four_wavelength")

snr_time <- seq.int(0, 400L - 1L) / curated_fs
snr_signal <- 2 * sin(2 * pi * 1 * snr_time) +
  0.2 * sin(2 * pi * 3 * snr_time)

live_time <- seq.int(0, 30L - 1L) / 10
live_values <- cbind(
  S1_D1 = 1 + ifelse(live_time >= 1, live_time, 0),
  S2_D2 = 1 + ifelse(live_time >= 1, 0.5 * live_time, 0),
  S3_D3 = 1
)

fixture <- list(
  schema_version = 1L,
  work_package = "WS10-16",
  license = "CC0-1.0 package-owned synthetic data",
  redistribution = paste(
    "Only deterministic package-owned synthetic values and upstream source",
    "identifiers are distributed; no upstream datasets or source code are",
    "included."
  ),
  references = list(
    pollonini_sci_doi = "10.1117/1.JBO.19.8.086007",
    pollonini_phoebe_doi = "10.1364/BOE.7.005104",
    kober_neurofeedback_doi = "10.3389/fnhum.2017.00081",
    mne = list(
      version = "1.3.1",
      source_url = paste0(
        "https://github.com/mne-tools/mne-python/blob/v1.3.1/",
        "mne/preprocessing/nirs/_scalp_coupling_index.py"
      ),
      source_sha256 =
        "403842797280f8e64d181763bfa1a48ba4b0e504108c9c2bb5852425e64e8315"
    ),
    mne_nirs = list(
      commit = "0a5081735144b902a3953e81d010420e1210c556",
      segmented_source_sha256 =
        "6c1b6a01adf01acedff097bf485c8d40116b77b839aaaaa5be66fa7792535e35",
      peak_power_source_sha256 =
        "da6ff09d847fb671d6dc922375d17ec23d753568837214c341a5c7a5eec9caf2"
    )
  ),
  generation = list(
    command = "Rscript data-raw/generate-quality-reference.R",
    validation_command = "Rscript inst/validation/validate-ws10-16.R",
    r_version = R.version.string,
    platform = R.version$platform,
    package_versions = c(
      digest = as.character(utils::packageVersion("digest")),
      signal = as.character(utils::packageVersion("signal")),
      PhysioStream = as.character(utils::packageVersion("PhysioStream"))
    ),
    validation_seed_range = 16001:16120
  ),
  curated = list(
    sampling_rate_hz = curated_fs,
    time_seconds = curated_time,
    wavelength_cases = curated,
    expected_class = c(coupled = TRUE, decoupled = FALSE)
  ),
  analytic_snr = list(
    sampling_rate_hz = curated_fs,
    time_seconds = snr_time,
    values = snr_signal,
    signal_band_hz = c(0.8, 1.2),
    noise_band_hz = matrix(c(2.8, 3.2), nrow = 1L),
    expected_snr_db = 20
  ),
  live = list(
    sampling_rate_hz = 10,
    time_seconds = live_time,
    values = live_values,
    channel_id = colnames(live_values),
    source_index = 1:3,
    detector_index = 1:3,
    regions = list(left = c("S1_D1", "S2_D2"), right = "S3_D3"),
    contrast = c(left = 1, right = -1),
    baseline_seconds = 1,
    update_seconds = 0.2,
    smoothing_seconds = 0.4,
    expected_baseline = c(left = 1, right = 1)
  )
)

rds_path <- file.path(output_dir, "quality_reference.rds")
saveRDS(fixture, rds_path, version = 3, compress = "xz")
rds_sha <- digest::digest(rds_path, algo = "sha256", file = TRUE)
generator_sha <- digest::digest(script, algo = "sha256", file = TRUE)
manifest <- c(
  paste(rds_sha, "quality_reference.rds"),
  paste(generator_sha, "generate-quality-reference.R"),
  paste(
    fixture$references$mne$source_sha256,
    "MNE 1.3.1 scalp-coupling source SHA256"
  ),
  paste(
    fixture$references$mne_nirs$commit,
    "MNE-NIRS reference commit"
  ),
  paste(
    fixture$references$mne_nirs$segmented_source_sha256,
    "MNE-NIRS segmented SCI source SHA256"
  ),
  paste(
    fixture$references$mne_nirs$peak_power_source_sha256,
    "MNE-NIRS peak-power source SHA256"
  )
)
writeLines(
  manifest,
  file.path(output_dir, "quality_reference.sha256"),
  useBytes = TRUE
)
cat("quality_reference.rds SHA256:", rds_sha, "\n")
