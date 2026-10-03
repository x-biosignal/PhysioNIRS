test_that("package-owned MBLL fixture is deterministic and governed", {
  path <- system.file(
    "extdata", "mbll_reference.rds", package = "PhysioNIRS"
  )
  expect_true(file.exists(path))
  fixture <- readRDS(path)
  expect_identical(fixture$schema_version, "1.0.0")
  expect_match(fixture$generator_sha256, "^[0-9a-f]{64}$")
  expect_match(fixture$extinction_table_sha256, "^[0-9a-f]{64}$")
  expect_identical(
    fixture$coefficient_convention$conversion,
    "100 * ln(10)"
  )
  expect_length(fixture$cases, 4L)
  expect_identical(
    vapply(fixture$cases, `[[`, character(1), "id"),
    c(
      "two_explicit", "three_scholkmann",
      "four_explicit", "three_duncan"
    )
  )
  expect_identical(
    vapply(
      fixture$cases,
      function(case) length(case$wavelengths_nm),
      integer(1)
    ),
    c(2L, 3L, 4L, 3L)
  )
  for (case in fixture$cases) {
    # Self-consistency: stored intensity must be derivable from stored
    # optical density (I = I0 * exp(-OD)). This recomputes exp(), a
    # transcendental whose final bit is libm/platform-dependent, so an
    # exact (tolerance = 0) compare is not portable across build toolchains.
    # Measured cross-platform noise from a 1-ULP exp() difference on these
    # fixtures is <= 4.3e-16 relative (~2x machine eps, ~4.5e-13 absolute on
    # intensities of order 1e3); 1e-12 sits ~2300x above that floor while
    # still failing any meaningful (>= 1e-6 relative) corruption.
    expect_equal(
      case$intensity,
      sweep(
        exp(-case$optical_density),
        2L,
        case$reference_intensity,
        "*"
      ),
      tolerance = 1e-12
    )
    expect_equal(
      case$HbT_uM,
      case$HbO_uM + case$HbR_uM,
      tolerance = 1e-15
    )
  }
})

test_that("motion reference fixture is package-owned and deterministic", {
  path <- system.file(
    "extdata", "motion_reference.rds", package = "PhysioNIRS"
  )
  manifest <- system.file(
    "extdata", "motion_reference.sha256", package = "PhysioNIRS"
  )
  expect_true(file.exists(path))
  expect_true(file.exists(manifest))
  expected <- strsplit(
    readLines(manifest, warn = FALSE)[[1L]], " "
  )[[1L]][1L]
  expect_identical(
    digest::digest(path, algo = "sha256", file = TRUE),
    expected
  )
  fixture <- readRDS(path)
  expect_identical(fixture$schema_version, 1L)
  expect_identical(fixture$work_package, "WS10-14")
  expect_identical(
    fixture$references$tddr$commit,
    "2b104674fdf39027f5148d7d97f61b60bad9327c"
  )
  expect_identical(
    fixture$references$homer3$commit,
    "a2bdfcf65e932478110cd9abdd4f0d1b773c5217"
  )
  expect_identical(
    fixture$references$csaps$commit,
    "4c1d003e822a3432cd52cd9e5a6c9662e966d0c9"
  )
  expect_match(fixture$references$tddr$url, "^https://")
  expect_match(fixture$generation$command, "generate-motion-reference")
  expect_match(fixture$generation$validation_command, "validate-ws10-14")
  expect_length(fixture$generation$validation_random_seeds, 100L)
  expect_match(fixture$redistribution, "no upstream", ignore.case = TRUE)
  expect_identical(dim(fixture$clean), c(3000L, 4L))
  expect_identical(dim(fixture$observed), c(3000L, 4L))
  expect_identical(dim(fixture$artifact_mask), c(3000L, 4L))
  expect_true(all(is.finite(fixture$observed)))
  expect_match(fixture$license, "package-owned")
})

test_that("short-separation reference and evidence are governed", {
  path <- system.file(
    "extdata", "short_separation_reference.rds", package = "PhysioNIRS"
  )
  manifest <- system.file(
    "extdata", "short_separation_reference.sha256", package = "PhysioNIRS"
  )
  expect_true(file.exists(path))
  expect_true(file.exists(manifest))
  lines <- readLines(manifest, warn = FALSE)
  expect_gte(length(lines), 5L)
  expected <- strsplit(lines[[1L]], " ", fixed = TRUE)[[1L]][[1L]]
  expect_identical(
    digest::digest(path, algo = "sha256", file = TRUE),
    expected
  )

  fixture <- readRDS(path)
  expect_identical(fixture$schema_version, 1L)
  expect_identical(fixture$work_package, "WS10-15")
  expect_identical(
    fixture$references$mne_nirs$commit,
    "0a5081735144b902a3953e81d010420e1210c556"
  )
  expect_match(
    fixture$references$mne_nirs$short_source_sha256,
    "^[0-9a-f]{64}$"
  )
  expect_match(
    fixture$references$mne_nirs$correction_source_sha256,
    "^[0-9a-f]{64}$"
  )
  expect_match(
    fixture$generation$command,
    "generate-short-separation-reference"
  )
  expect_match(fixture$generation$validation_command, "validate-ws10-15")
  expect_length(fixture$generation$validation_random_seeds, 100L)
  expect_match(fixture$license, "package-owned")
  expect_match(fixture$redistribution, "no upstream", ignore.case = TRUE)
  expect_identical(dim(fixture$observed), c(2400L, 6L))
  expect_identical(dim(fixture$oracle_corrected), c(2400L, 6L))
  expect_equal(
    colMeans(fixture$oracle_corrected),
    colMeans(fixture$observed),
    tolerance = 1e-14
  )
  expect_true(all(is.finite(fixture$oracle_corrected)))

  validation <- utils::read.csv(system.file(
    "validation", "ws10-15-validation.csv", package = "PhysioNIRS"
  ))
  mutations <- utils::read.csv(system.file(
    "validation", "ws10-15-mutations.csv", package = "PhysioNIRS"
  ))
  report <- system.file(
    "validation", "ws10-15-validation.md", package = "PhysioNIRS"
  )
  expect_identical(nrow(validation), 100L)
  expect_gte(nrow(mutations), 25L)
  expect_true(all(mutations$detected))
  expect_true(file.exists(report))
})

test_that("quality and neurofeedback reference evidence is governed", {
  path <- system.file(
    "extdata", "quality_reference.rds", package = "PhysioNIRS"
  )
  manifest <- system.file(
    "extdata", "quality_reference.sha256", package = "PhysioNIRS"
  )
  expect_true(file.exists(path))
  expect_true(file.exists(manifest))
  lines <- readLines(manifest, warn = FALSE)
  expect_gte(length(lines), 6L)
  expected <- strsplit(lines[[1L]], " ", fixed = TRUE)[[1L]][[1L]]
  expect_identical(
    digest::digest(path, algo = "sha256", file = TRUE),
    expected
  )

  fixture <- readRDS(path)
  expect_identical(fixture$schema_version, 1L)
  expect_identical(fixture$work_package, "WS10-16")
  expect_identical(
    fixture$references$mne_nirs$commit,
    "0a5081735144b902a3953e81d010420e1210c556"
  )
  expect_match(fixture$references$mne$source_sha256, "^[0-9a-f]{64}$")
  expect_match(
    fixture$references$mne_nirs$peak_power_source_sha256,
    "^[0-9a-f]{64}$"
  )
  expect_match(fixture$generation$command, "generate-quality-reference")
  expect_match(fixture$generation$validation_command, "validate-ws10-16")
  expect_length(fixture$generation$validation_seed_range, 120L)
  expect_identical(
    names(fixture$curated$wavelength_cases),
    c("two_wavelength", "three_wavelength", "four_wavelength")
  )
  expect_identical(fixture$analytic_snr$expected_snr_db, 20)
  expect_identical(
    colnames(fixture$live$values), fixture$live$channel_id
  )
  expect_match(fixture$license, "package-owned")
  expect_match(fixture$redistribution, "no upstream", ignore.case = TRUE)
})
