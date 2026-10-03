#' Build a tiny synthetic continuous-wave fNIRS object for examples
#'
#' Constructs a governed continuous-wave `PhysioExperiment` with two
#' source-detector pairs (one long, one short) measured at two wavelengths,
#' deterministic cardiac, respiration, drift, and task components, a single
#' injected motion spike, and a `"task"` stimulus. It exists only so help-page
#' examples have a self-contained, offline object to operate on; it reads no
#' files and draws no random numbers.
#'
#' @param n_time Number of samples.
#' @param fs Sampling rate in Hz.
#' @return A continuous-wave intensity `PhysioExperiment`.
#' @keywords internal
#' @noRd
.nirs_demo_object <- function(n_time = 160L, fs = 10) {
  t <- seq_len(n_time) / fs
  cardiac <- 0.06 * sin(2 * pi * 1.0 * t)
  resp <- 0.03 * sin(2 * pi * 0.25 * t)
  drift <- 0.02 * (t - mean(t)) / max(t)
  task <- 0.08 * exp(-((t - 8)^2) / (2 * 2^2))
  base_od <- c(1.2, 1.0, 1.1, 0.9)
  wl_gain <- c(1.0, 0.8, 1.0, 0.8)
  od <- matrix(0, n_time, 4L)
  for (j in seq_len(4L)) {
    od[, j] <- base_od[j] + cardiac + resp + drift + wl_gain[j] * task
  }
  od[round(n_time / 2), ] <- od[round(n_time / 2), ] + 0.7
  intensity <- sweep(exp(-od), 2L, c(5000, 4000, 5200, 4100), "*")
  wavelengths <- c(760, 850)
  source_index <- c(1L, 1L, 2L, 2L)
  detector_index <- c(1L, 1L, 2L, 2L)
  wavelength_index <- c(1L, 2L, 1L, 2L)
  channel_label <- c("S1_D1_760", "S1_D1_850", "S2_D2_760", "S2_D2_850")
  colnames(intensity) <- channel_label
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(4L),
    source_index = source_index,
    detector_index = detector_index,
    wavelength_index = wavelength_index,
    wavelength_nm = wavelengths[wavelength_index],
    wavelength_actual_nm = rep(NA_real_, 4L),
    data_type = rep(1L, 4L),
    data_type_index = rep(1L, 4L),
    data_type_label = rep(NA_character_, 4L),
    data_unit = rep(NA_character_, 4L),
    source_power = rep(NA_real_, 4L),
    detector_gain = rep(NA_real_, 4L),
    source_label = paste0("S", source_index),
    detector_label = paste0("D", detector_index),
    channel_label = channel_label
  )
  probe <- list(
    wavelengths = wavelengths,
    sourceLabels = c("S1", "S2"),
    detectorLabels = c("D1", "D2"),
    sourcePos2D = rbind(c(0, 0), c(0.1, 0)),
    detectorPos2D = rbind(c(0, 0.03), c(0.1, 0.008)),
    LengthUnit = "m"
  )
  x <- PhysioExperiment::PhysioExperiment(
    assays = list(raw = intensity),
    rowData = S4Vectors::DataFrame(time_seconds = t),
    colData = cbind(measurement, label = channel_label),
    metadata = list(snirf = list(
      measurement_list = measurement, probe = probe,
      metadata_tags = list(LengthUnit = "m", TimeUnit = "s"),
      stim = list(), time_sampling = "uniform"
    )),
    samplingRate = fs
  )
  PhysioExperiment::setEvents(
    x,
    PhysioExperiment::PhysioEvents(
      onset = 6, duration = 4, type = "task", value = "task"
    )
  )
}
