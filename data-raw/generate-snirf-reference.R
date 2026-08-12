# Generate the package-owned, non-identifying SNIRF reference fixture.
#
# Run from the repository root:
# Rscript physio-ecosystem/PhysioNIRS/data-raw/generate-snirf-reference.R

package_dir <- normalizePath("physio-ecosystem/PhysioNIRS")
devtools::load_all(package_dir, quiet = TRUE)
source(file.path(package_dir, "tests/testthat/helper-snirf.R"), local = TRUE)

x <- make_snirf_experiment()
S4Vectors::metadata(x)$snirf$metadata_tags$SubjectID <- "anonymous"
output <- file.path(package_dir, "inst/extdata/snirf_reference.snirf")
writeSNIRF(x, output, overwrite = TRUE)

semantic <- list(
  data = unname(SummarizedExperiment::assay(x)),
  time = as.numeric(SummarizedExperiment::rowData(x)$time_seconds),
  measurement = as.list(measurementList(x)),
  probe = S4Vectors::metadata(x)$snirf$probe,
  stim = S4Vectors::metadata(x)$snirf$stim,
  tags = S4Vectors::metadata(x)$snirf$metadata_tags
)
semantic_file <- tempfile(fileext = ".rds")
on.exit(unlink(semantic_file), add = TRUE)
saveRDS(semantic, semantic_file, version = 2, compress = FALSE)
semantic_digest <- unname(tools::md5sum(semantic_file))
sha256 <- unname(system2("sha256sum", output, stdout = TRUE))
sha256 <- sub(" .*", "", sha256)

dcf <- data.frame(
  Generator = "Rscript physio-ecosystem/PhysioNIRS/data-raw/generate-snirf-reference.R",
  R = paste(R.version$major, R.version$minor, sep = "."),
  rhdf5 = as.character(utils::packageVersion("rhdf5")),
  HDF5 = Rhdf5lib::getHdf5Version(),
  SNIRF = "1.0 (v1.1 contract)",
  Dimensions = "17 time x 6 measurements",
  Wavelengths = "760, 850 nm",
  Geometry = "2-D and 3-D, cm",
  Stimuli = "left (4 columns), right (3 columns)",
  License = "CC0/public domain synthetic fixture",
  Provenance = "Deterministic non-identifying analytic values",
  SemanticDigest = semantic_digest,
  stringsAsFactors = FALSE
)
write.dcf(dcf, file.path(package_dir, "inst/extdata/snirf_reference.dcf"))
writeLines(
  paste(sha256, " snirf_reference.snirf"),
  file.path(package_dir, "inst/extdata/snirf_reference.sha256")
)
