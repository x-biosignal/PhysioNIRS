test_that("package metadata and public surface are governed", {
  description <- utils::packageDescription("PhysioNIRS")
  expect_identical(description$Package, "PhysioNIRS")
  expect_identical(description$Version, "0.5.1")
  authors <- eval(parse(text = description$`Authors@R`))
  expect_length(authors, 1L)
  expect_identical(authors$given, "Yusuke")
  expect_identical(authors$family, "Matsui")
  expect_identical(authors$role, c("aut", "cre"))
  expect_identical(authors$email, "mail.to.matsui@gmail.com")
  expect_match(description$Imports, "PhysioStream \\(>= 0.9.0\\)")
  expect_match(
    description$Collate,
    paste0(
      "'PhysioNIRS-package.R'.*'snirf-contract.R'.*'io-snirf.R'.*",
      "'probe.R'.*'optical-density.R'.*'dpf.R'.*'extinction.R'.*",
      "'mbll.R'.*'motion-contract.R'.*'motion-detect.R'.*",
      "'motion-tddr.R'.*'motion-wavelet.R'.*'motion-spline.R'.*",
      "'short-separation-contract.R'.*'short-separation.R'.*",
      "'physiology-band.R'.*'nuisance-design.R'.*",
      "'quality-contract.R'.*'quality.R'.*'neurofeedback.R'.*'zzz.R'"
    )
  )
  expect_setequal(
    getNamespaceExports("PhysioNIRS"),
    c(
      "readSNIRF", "writeSNIRF", "measurementList",
      "sourceDetectorDistances", "intensityToOD", "dpf",
      "extinctionCoefficients", "mbll", "motionArtifactDetect", "tddr",
      "waveletMotionCorrect", "splineMotionCorrect",
      "identifyShortChannels", "shortSeparationRegress",
      "physiologyBandpass", "shortSeparationDesign",
      "scalpCouplingIndex", "signalQualityIndex", "pruneChannels",
      "nirsNeurofeedback", "nirsNeurofeedbackStart",
      "nirsNeurofeedbackStep", "nirsNeurofeedbackState",
      "nirsNeurofeedbackStop", "nirsNeurofeedbackScope"
    )
  )
})

test_that("the canonical public inventory includes PhysioNIRS", {
  repository <- normalizePath(
    test_path("..", "..", ".."), mustWork = FALSE
  )
  install_path <- file.path(repository, "install_ecosystem.R")
  check_path <- file.path(repository, "publishing", "check-ecosystem.yml")
  sync_path <- file.path(repository, "publishing", "sync_public.sh")
  packages_path <- file.path(repository, "publishing", "packages.json")
  paths <- c(install_path, check_path, sync_path, packages_path)
  skip_if_not(all(file.exists(paths)), "monorepo publishing files are not installed")

  expect_match(paste(readLines(install_path), collapse = "\n"), "PhysioNIRS")
  expect_match(paste(readLines(packages_path), collapse = "\n"),
               '"package": "PhysioNIRS"', fixed = TRUE)
  expect_match(paste(readLines(check_path), collapse = "\n"),
               "build_ecosystem.R", fixed = TRUE)
  sync_text <- paste(readLines(sync_path), collapse = "\n")
  expect_match(sync_text, "packages.json", fixed = TRUE)
  expect_match(sync_text, "ALL_PACKAGES", fixed = TRUE)
})

test_that("print output does not disclose subject identifiers", {
  x <- make_snirf_experiment()
  output <- capture.output(print(x))
  expect_false(any(grepl("fixture-subject", output, fixed = TRUE)))
})
