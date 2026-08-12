make_short_nirs <- function(
    n_time = 1800L,
    fs = 10,
    distances = c(0.008, 0.009, 0.03),
    values = NULL) {
  time <- seq.int(0, n_time - 1L) / fs
  if (is.null(values)) {
    short_1 <- sin(2 * pi * 0.10 * time) +
      0.1 * sin(2 * pi * 0.27 * time)
    short_2 <- cos(2 * pi * 0.09 * time) +
      0.1 * sin(2 * pi * 0.31 * time)
    neural <- sin(2 * pi * 0.035 * time)
    values <- cbind(
      short_1,
      0.8 * short_1 + 0.05 * sin(2 * pi * 0.4 * time),
      short_2,
      0.7 * short_2 + 0.05 * cos(2 * pi * 0.4 * time),
      neural + 1.8 * short_1,
      0.7 * neural + 1.2 *
        (0.8 * short_1 + 0.05 * sin(2 * pi * 0.4 * time))
    )
  }
  x <- make_cw_nirs(
    wavelengths = c(760, 850),
    distances = distances,
    n_time = nrow(values),
    assay_name = "OD",
    values = values
  )
  SummarizedExperiment::rowData(x)$time_seconds <-
    seq.int(0, nrow(values) - 1L) / fs
  methods::slot(x, "samplingRate") <- as.numeric(fs)
  x
}

make_short_hb <- function(n_time = 600L, fs = 10) {
  time <- seq.int(0, n_time - 1L) / fs
  od <- cbind(
    0.01 * sin(2 * pi * 0.10 * time),
    0.009 * sin(2 * pi * 0.10 * time + 0.1),
    0.02 * sin(2 * pi * 0.10 * time) +
      0.005 * sin(2 * pi * 0.03 * time),
    0.018 * sin(2 * pi * 0.10 * time + 0.1) +
      0.004 * sin(2 * pi * 0.03 * time)
  )
  x <- make_cw_nirs(
    wavelengths = c(760, 850),
    distances = c(0.008, 0.03),
    n_time = n_time,
    assay_name = "OD",
    values = od
  )
  SummarizedExperiment::rowData(x)$time_seconds <- time
  methods::slot(x, "samplingRate") <- as.numeric(fs)
  suppressWarnings(mbll(x, pathlength_factor = 6))
}

short_assay <- function(x, name) {
  SummarizedExperiment::assay(x, name, withDimnames = FALSE)
}
