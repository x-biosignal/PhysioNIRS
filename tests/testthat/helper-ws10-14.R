make_motion_nirs <- function(
    n_time = 128L,
    fs = 10,
    values = NULL,
    assay_name = "OD") {
  n_measurement <- 4L
  if (is.null(values)) {
    time <- seq.int(0, n_time - 1L) / fs
    values <- cbind(
      sin(2 * pi * 0.08 * time),
      0.8 * sin(2 * pi * 0.08 * time + 0.1),
      0.5 * sin(2 * pi * 0.1 * time),
      0.4 * sin(2 * pi * 0.1 * time + 0.2)
    )
  }
  values <- as.matrix(values)
  stopifnot(ncol(values) == n_measurement)
  x <- make_cw_nirs(
    wavelengths = c(760, 850),
    distances = c(0.03, 0.04),
    n_time = nrow(values),
    assay_name = assay_name,
    values = values
  )
  time <- seq.int(0, nrow(values) - 1L) / fs
  SummarizedExperiment::rowData(x)$time_seconds <- time
  methods::slot(x, "samplingRate") <- as.numeric(fs)
  metadata <- S4Vectors::metadata(x)
  metadata$snirf$time_sampling <- "uniform"
  if (assay_name == "OD") {
    metadata$nirs <- list(assays = list(OD = list(
      kind = "optical_density",
      unit = "1",
      log_convention = "natural",
      source_assay = "synthetic"
    )))
  }
  S4Vectors::metadata(x) <- metadata
  x
}

motion_assay <- function(x, name) {
  SummarizedExperiment::assay(x, name, withDimnames = FALSE)
}

make_manual_motion_mask <- function(x, bad, assay_name = "OD") {
  context <- PhysioNIRS:::.nirs_motion_context(x, assay_name)
  bad <- as.matrix(bad)
  stopifnot(identical(dim(bad), dim(context$data)))
  storage.mode(bad) <- "logical"
  out <- list(
    sample_by_measurement = bad,
    global = apply(bad, 1L, any),
    intervals = PhysioNIRS:::.nirs_motion_intervals(
      bad, context$time, context$identity
    ),
    parameters = list(source = "test fixture"),
    assay_name = assay_name,
    sampling_rate_hz = context$fs,
    measurement_id = as.character(context$identity$measurement_id),
    input_fingerprint = context$fingerprint
  )
  class(out) <- "nirs_motion_mask"
  out$fingerprint <- PhysioNIRS:::.nirs_motion_mask_fingerprint(out)
  out
}
