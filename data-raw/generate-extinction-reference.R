# Regenerate the governed extinction table from the BSD-3-Clause MNE-Python
# source asset. This script is not part of the installed package.

upstream_commit <- "b1fb8933a2ba62028b1dac3ae8bc24da0220d00e"
upstream_url <- paste0(
  "https://raw.githubusercontent.com/mne-tools/mne-python/",
  upstream_commit,
  "/mne/data/extinction_coef.mat"
)
upstream_sha256 <-
  "6b7fe0d64af92e542741b1a05e9c4a02233d088329e6e9f6ef9c5896f8871843"
retrieval_date <- "2026-07-30"
source_label <- "MNE-Python extinction_coef.mat (Prahl/OMLC)"

root <- normalizePath("physio-ecosystem/PhysioNIRS")
output <- file.path(root, "inst", "extdata", "haemoglobin-extinction.csv")
hash_output <- file.path(
  root, "inst", "extdata", "haemoglobin-extinction.sha256"
)
mat_path <- tempfile(fileext = ".mat")
helper <- tempfile(fileext = ".py")
on.exit(unlink(c(mat_path, helper)), add = TRUE)

utils::download.file(upstream_url, mat_path, mode = "wb", quiet = TRUE)
actual_upstream <- digest::digest(mat_path, algo = "sha256", file = TRUE)
stopifnot(identical(actual_upstream, upstream_sha256))

python <- Sys.getenv("PHYSIONIRS_PYTHON", "python3")
writeLines(c(
  "import csv, sys",
  "from scipy.io import loadmat",
  "mat_path, output, source, version = sys.argv[1:]",
  "values = loadmat(mat_path)['extinct_coef']",
  "with open(output, 'w', newline='', encoding='ascii') as stream:",
  "    writer = csv.writer(stream, lineterminator='\\n')",
  "    writer.writerow(['wavelength_nm', 'hbo2_cm1_M1', 'hbr_cm1_M1', 'source', 'source_version'])",
  "    for wavelength, hbo, hbr in values:",
  "        writer.writerow([format(wavelength, '.15g'), format(hbo, '.15g'), format(hbr, '.15g'), source, version])"
), helper)
status <- system2(
  python,
  shQuote(c(helper, mat_path, output, source_label, upstream_commit))
)
stopifnot(identical(status, 0L))

table <- utils::read.csv(output, stringsAsFactors = FALSE)
stopifnot(
  identical(
    names(table),
    c(
      "wavelength_nm", "hbo2_cm1_M1", "hbr_cm1_M1",
      "source", "source_version"
    )
  ),
  nrow(table) == 376L,
  isTRUE(all.equal(range(table$wavelength_nm), c(250, 1000))),
  all(diff(table$wavelength_nm) == 2),
  all(table$source == source_label),
  all(table$source_version == upstream_commit)
)
output_sha256 <- digest::digest(output, algo = "sha256", file = TRUE)
writeLines(paste(output_sha256, basename(output)), hash_output, useBytes = TRUE)

message("retrieval_date=", retrieval_date)
message("upstream_url=", upstream_url)
message("upstream_sha256=", upstream_sha256)
message("rows=", nrow(table), " range=", paste(range(table$wavelength_nm),
                                               collapse = "-"))
message("output_sha256=", output_sha256)
