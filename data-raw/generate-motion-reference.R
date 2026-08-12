#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- if (length(file_arg)) {
  normalizePath(file_arg[[1L]], mustWork = TRUE)
} else {
  normalizePath(
    sys.frame(1)$ofile %||% "data-raw/generate-motion-reference.R",
    mustWork = TRUE
  )
}
package_root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
output_dir <- file.path(package_root, "inst", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fs <- 10
time <- seq.int(0, 300 * fs - 1L) / fs
n <- length(time)
clean <- cbind(
  S1_D1_760 = 0.08 * sin(2 * pi * 0.05 * time) +
    0.025 * sin(2 * pi * 0.10 * time),
  S1_D1_850 = 0.07 * sin(2 * pi * 0.05 * time + 0.04) +
    0.020 * sin(2 * pi * 0.10 * time),
  S2_D2_760 = 0.09 * sin(2 * pi * 0.05 * time + 0.08) +
    0.020 * sin(2 * pi * 0.10 * time + 0.03),
  S2_D2_850 = 0.08 * sin(2 * pi * 0.05 * time + 0.12) +
    0.018 * sin(2 * pi * 0.10 * time + 0.03)
)
deterministic_noise <- outer(
  time,
  seq_len(ncol(clean)),
  function(t, j) 0.002 * sin(2 * pi * (0.31 + 0.013 * j) * t + j)
)
observed <- clean + deterministic_noise
artifact_mask <- matrix(
  FALSE, nrow = n, ncol = ncol(clean), dimnames = dimnames(clean)
)

observed[801, 1:2] <- observed[801, 1:2] + c(1.4, 1.2)
artifact_mask[796:806, 1:2] <- TRUE
observed[1401:n, 1] <- observed[1401:n, 1] + 0.9
observed[1401:n, 2] <- observed[1401:n, 2] + 0.75
artifact_mask[1396:1406, 1:2] <- TRUE
burst <- 0.8 * sin(seq(0, 8 * pi, length.out = 31))
observed[2100:2130, 3] <- observed[2100:2130, 3] + burst
observed[2100:2130, 4] <- observed[2100:2130, 4] + 0.85 * burst
artifact_mask[2095:2135, 3:4] <- TRUE

fixture <- list(
  schema_version = 1L,
  work_package = "WS10-14",
  license = "CC0-1.0 package-owned synthetic data",
  references = list(
    tddr = list(
      url = "https://github.com/frankfishburn/TDDR",
      commit = "2b104674fdf39027f5148d7d97f61b60bad9327c"
    ),
    homer3 = list(
      url = "https://github.com/BUNPC/Homer3",
      commit = "a2bdfcf65e932478110cd9abdd4f0d1b773c5217"
    ),
    csaps = list(
      url = "https://github.com/espdev/csaps",
      commit = "4c1d003e822a3432cd52cd9e5a6c9662e966d0c9"
    )
  ),
  generation = list(
    command = "Rscript data-raw/generate-motion-reference.R",
    validation_command = "Rscript inst/validation/validate-ws10-14.R",
    external_reference_command = paste(
      "none; pinned sources were inspected and independently ported;",
      "no external output is copied into this fixture"
    ),
    r_version = R.version.string,
    platform = R.version$platform,
    package_versions = c(
      digest = as.character(utils::packageVersion("digest")),
      signal = as.character(utils::packageVersion("signal"))
    ),
    fixture_random_seed = NA_integer_,
    validation_random_seeds = 14000L + seq_len(100L)
  ),
  redistribution = paste(
    "Only package-owned deterministic synthetic values are distributed;",
    "no upstream source code, datasets, or generated reference outputs",
    "are included."
  ),
  parameters = list(
    sampling_rate_hz = fs,
    duration_seconds = 300,
    wavelengths_nm = c(760, 850),
    source_index = c(1L, 1L, 2L, 2L),
    detector_index = c(1L, 1L, 2L, 2L),
    source_detector_distance_m = c(0.03, 0.03, 0.04, 0.04),
    t_motion = 0.5,
    t_mask = 1,
    sd_threshold = 20,
    amplitude_threshold = 0.5,
    tddr_cutoff_hz = 0.5,
    wavelet_iqr = 1.5,
    spline_p = 0.99
  ),
  time_seconds = time,
  clean = clean,
  observed = observed,
  artifact_mask = artifact_mask
)

rds_path <- file.path(output_dir, "motion_reference.rds")
saveRDS(fixture, rds_path, version = 3, compress = "xz")
rds_sha <- digest::digest(rds_path, algo = "sha256", file = TRUE)
generator_sha <- digest::digest(script, algo = "sha256", file = TRUE)
manifest <- c(
  paste(rds_sha, "motion_reference.rds"),
  paste(generator_sha, "generate-motion-reference.R"),
  paste(
    "2b104674fdf39027f5148d7d97f61b60bad9327c",
    "TDDR reference commit"
  ),
  paste(
    "a2bdfcf65e932478110cd9abdd4f0d1b773c5217",
    "Homer3 reference commit"
  ),
  paste(
    "4c1d003e822a3432cd52cd9e5a6c9662e966d0c9",
    "csaps reference commit"
  )
)
writeLines(
  manifest,
  file.path(output_dir, "motion_reference.sha256"),
  useBytes = TRUE
)
cat("motion_reference.rds SHA256:", rds_sha, "\n")
