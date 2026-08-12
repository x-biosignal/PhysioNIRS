#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- if (length(file_arg)) {
  normalizePath(file_arg[[1L]], mustWork = TRUE)
} else {
  normalizePath(
    sys.frame(1)$ofile %||%
      "data-raw/generate-short-separation-reference.R",
    mustWork = TRUE
  )
}
package_root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
output_dir <- file.path(package_root, "inst", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fs <- 10
time <- seq.int(0, 240 * fs - 1L) / fs
task <- as.numeric(time >= 80 & time < 140)
superficial_1 <- sin(2 * pi * 0.10 * time) +
  0.25 * sin(2 * pi * 0.28 * time + 0.2)
superficial_2 <- cos(2 * pi * 0.09 * time + 0.4) +
  0.20 * sin(2 * pi * 0.32 * time)
neural <- 0.4 * task + 0.08 * sin(2 * pi * 0.03 * time)
noise <- outer(
  time,
  seq_len(6L),
  function(t, j) 0.01 * sin(2 * pi * (0.41 + 0.013 * j) * t + j)
)
observed <- cbind(
  superficial_1,
  0.82 * superficial_1,
  superficial_2,
  0.75 * superficial_2,
  neural + 1.6 * superficial_1,
  0.7 * neural + 1.3 * 0.82 * superficial_1
) + noise
colnames(observed) <- c(
  "S1_D1_760", "S1_D1_850", "S2_D2_760",
  "S2_D2_850", "S3_D3_760", "S3_D3_850"
)

source_position_m <- rbind(c(0, 0), c(0.08, 0), c(0.018, 0))
detector_position_m <- rbind(c(0, 0.008), c(0.08, 0.009), c(0.018, 0.03))
midpoint_m <- (source_position_m + detector_position_m) / 2
distance_m <- sqrt(rowSums((source_position_m - detector_position_m)^2))
nearest_pair <- vapply(3L, function(target) {
  candidates <- 1:2
  candidates[[which.min(sqrt(rowSums(
    sweep(midpoint_m[candidates, , drop = FALSE], 2L,
          midpoint_m[target, ], "-")^2
  )))]]
}, integer(1))

oracle_corrected <- observed
for (target in 5:6) {
  predictor <- if (target == 5L) 1L else 2L
  beta <- stats::cov(observed[, target], observed[, predictor]) /
    stats::var(observed[, predictor])
  oracle_corrected[, target] <- observed[, target] -
    (observed[, predictor] - mean(observed[, predictor])) * beta
}

fixture <- list(
  schema_version = 1L,
  work_package = "WS10-15",
  license = "CC0-1.0 package-owned synthetic data",
  redistribution = paste(
    "Only deterministic package-owned synthetic values are distributed;",
    "no upstream source code or datasets are included."
  ),
  references = list(
    mne_nirs = list(
      url = "https://github.com/mne-tools/mne-nirs",
      commit = "0a5081735144b902a3953e81d010420e1210c556",
      short_source_sha256 =
        "c9c0f621301e6c0aafca29bb9e987e2d455cf06c59ae1f8e3418e6c39f88991f",
      correction_source_sha256 =
        "0d8c34bd559e49839ff41a8b3699a556bb4efb0b17bf6ecc2c0867b875c406c4"
    ),
    saager_berger_doi = "10.1364/JOSAA.22.001874",
    gagnon_doi = "10.1016/j.neuroimage.2011.08.095",
    scholkmann_doi = "10.1016/j.neuroimage.2013.05.004"
  ),
  generation = list(
    command = "Rscript data-raw/generate-short-separation-reference.R",
    validation_command = "Rscript inst/validation/validate-ws10-15.R",
    r_version = R.version.string,
    platform = R.version$platform,
    package_versions = c(
      digest = as.character(utils::packageVersion("digest")),
      signal = as.character(utils::packageVersion("signal"))
    ),
    fixture_random_seed = NA_integer_,
    validation_random_seeds = 15000L + seq_len(100L)
  ),
  parameters = list(
    sampling_rate_hz = fs,
    threshold_m = 0.01,
    wavelengths_nm = c(760, 850),
    source_index = rep(1:3, each = 2L),
    detector_index = rep(1:3, each = 2L),
    distance_m = rep(distance_m, each = 2L),
    source_position_m = source_position_m,
    detector_position_m = detector_position_m,
    nearest_short_pair = nearest_pair
  ),
  time_seconds = time,
  task = task,
  neural = neural,
  superficial = cbind(superficial_1, superficial_2),
  observed = observed,
  oracle_corrected = oracle_corrected
)

rds_path <- file.path(output_dir, "short_separation_reference.rds")
saveRDS(fixture, rds_path, version = 3, compress = "xz")
rds_sha <- digest::digest(rds_path, algo = "sha256", file = TRUE)
generator_sha <- digest::digest(script, algo = "sha256", file = TRUE)
manifest <- c(
  paste(rds_sha, "short_separation_reference.rds"),
  paste(generator_sha, "generate-short-separation-reference.R"),
  paste(
    "0a5081735144b902a3953e81d010420e1210c556",
    "MNE-NIRS reference commit"
  ),
  paste(
    "c9c0f621301e6c0aafca29bb9e987e2d455cf06c59ae1f8e3418e6c39f88991f",
    "MNE-NIRS _short.py SHA256"
  ),
  paste(
    "0d8c34bd559e49839ff41a8b3699a556bb4efb0b17bf6ecc2c0867b875c406c4",
    "MNE-NIRS _short_channel_correction.py SHA256"
  )
)
writeLines(
  manifest,
  file.path(output_dir, "short_separation_reference.sha256"),
  useBytes = TRUE
)
cat("short_separation_reference.rds SHA256:", rds_sha, "\n")
