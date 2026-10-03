make_quality_nirs <- function(
    n_time = 600L,
    fs = 10,
    wavelengths = c(760, 850),
    coupled = c(TRUE, FALSE),
    offset = 0,
    scale = 1) {
  time <- seq.int(0, n_time - 1L) / fs
  values <- matrix(
    0, nrow = n_time, ncol = length(wavelengths) * length(coupled)
  )
  for (pair in seq_along(coupled)) {
    base <- sin(2 * pi * 1.0 * time + pair / 10) +
      0.05 * sin(2 * pi * 0.2 * time)
    for (w in seq_along(wavelengths)) {
      column <- (pair - 1L) * length(wavelengths) + w
      values[, column] <- if (coupled[[pair]]) {
        base + 0.01 * sin(2 * pi * (2 + w / 10) * time)
      } else {
        sin(2 * pi * (0.15 + 0.12 * w) * time + w)
      }
    }
  }
  values <- offset + scale * values
  x <- make_cw_nirs(
    wavelengths = wavelengths,
    distances = rep(0.03, length(coupled)),
    n_time = n_time,
    assay_name = "OD",
    values = values
  )
  SummarizedExperiment::rowData(x)$time_seconds <- time
  methods::slot(x, "samplingRate") <- as.numeric(fs)
  x
}

quality_assay <- function(x, name = "OD") {
  SummarizedExperiment::assay(x, name, withDimnames = FALSE)
}

make_feedback_source <- function(
    fs = 10,
    channel_id = c("S1_D1", "S2_D2", "S3_D3"),
    open = TRUE,
    metadata = NULL,
    units = rep("uM", length(channel_id))) {
  if (is.null(metadata)) {
    metadata <- list(nirs = list(
      assay_name = "HbO",
      assay_kind = "haemoglobin_concentration",
      identity_kind = "pair",
      channel_id = channel_id,
      source_index = seq_along(channel_id),
      detector_index = seq_along(channel_id)
    ))
  }
  info <- PhysioStream::streamInfo(
    "nirs-feedback-fixture",
    type = "NIRS",
    channel_names = channel_id,
    nominal_srate = fs,
    dtype = "float64",
    source_id = "nirs-feedback-fixture",
    clock_domain = "nirs-fixture-clock",
    channel_units = units,
    metadata = metadata
  )
  source <- PhysioStream::loopbackSource(info, 8192L)
  if (open) PhysioStream::streamOpen(source) else source
}

make_feedback_controller <- function(
    source = make_feedback_source(),
    contrast = c(left = 1, right = -1),
    baseline_seconds = 1,
    update_seconds = 0.2,
    smoothing_seconds = 0.4) {
  nirsNeurofeedback(
    source,
    regions = list(
      left = c("S1_D1", "S2_D2"),
      right = "S3_D3"
    ),
    contrast = contrast,
    baseline_seconds = baseline_seconds,
    update_seconds = update_seconds,
    smoothing_seconds = smoothing_seconds
  )
}

feedback_values <- function(time) {
  cbind(
    S1_D1 = 1 + ifelse(time >= 1, time, 0),
    S2_D2 = 1 + ifelse(time >= 1, 0.5 * time, 0),
    S3_D3 = 1
  )
}

run_feedback_fixture <- function(n_time = 30L, fs = 10, max_chunks = 4L) {
  source <- make_feedback_source(fs = fs)
  controller <- make_feedback_controller(source)
  nirsNeurofeedbackStart(controller)
  time <- seq.int(0, n_time - 1L) / fs
  PhysioStream::loopbackFeed(source, feedback_values(time), time)
  repeat {
    result <- nirsNeurofeedbackStep(controller, max_chunks)
    if (!result$updated) break
  }
  list(
    source = source,
    controller = controller,
    state = nirsNeurofeedbackState(controller),
    time = time
  )
}
