.nirs_feedback_schema <- "1.0.0"
.nirs_feedback_update_limit <- 512L
.nirs_feedback_receipt_limit <- 1024L
.nirs_feedback_registry <- new.env(parent = emptyenv())
.nirs_feedback_registry$counter <- 0

.nirs_feedback_runtime <- function(controller) {
  if (!inherits(controller, "NIRSNeurofeedback") ||
      !is.environment(controller) ||
      !identical(names(controller), "id") ||
      !bindingIsLocked("id", controller) ||
      !is.character(controller$id) || length(controller$id) != 1L ||
      is.na(controller$id) || !nzchar(controller$id) ||
      !exists(
        controller$id, envir = .nirs_feedback_registry, inherits = FALSE
      )) {
    stop("`controller` must be a valid NIRSNeurofeedback runtime",
         call. = FALSE)
  }
  get(controller$id, envir = .nirs_feedback_registry, inherits = FALSE)
}

.nirs_feedback_assert <- function(controller) {
  runtime <- .nirs_feedback_runtime(controller)
  if (!is.environment(runtime) ||
      !inherits(runtime$scope, "BiofeedbackScope") ||
      !inherits(runtime$pipeline, "StreamPipeline") ||
      !is.logical(runtime$busy) || length(runtime$busy) != 1L ||
      is.na(runtime$busy)) {
    stop("`controller` must be a valid NIRSNeurofeedback runtime",
         call. = FALSE)
  }
  invisible(runtime)
}

.nirs_feedback_strings <- function(x, arg) {
  if (!is.character(x) || !is.null(dim(x)) || !length(x) || anyNA(x) ||
      any(!nzchar(x)) || anyDuplicated(x)) {
    stop("`", arg, "` must contain unique non-empty strings",
         call. = FALSE)
  }
  as.character(x)
}

.nirs_feedback_contract <- function(stream, assay_name) {
  if (!methods::is(stream, "StreamSource")) {
    stop("`stream` must be a PhysioStream StreamSource", call. = FALSE)
  }
  info <- PhysioStream::streamInfo(stream)
  if (!is.numeric(info@nominal_srate) ||
      length(info@nominal_srate) != 1L ||
      !is.finite(info@nominal_srate) || info@nominal_srate <= 0 ||
      identical(info@dtype, "string")) {
    stop("NIRS neurofeedback requires a regular-rate numeric stream",
         call. = FALSE)
  }
  assay_name <- .nirs_scalar_name(assay_name, "assay_name")
  contract <- info@metadata$nirs
  required <- c(
    "assay_name", "assay_kind", "identity_kind", "channel_id",
    "source_index", "detector_index"
  )
  if (!is.list(contract) || !all(required %in% names(contract)) ||
      !identical(contract$assay_name, assay_name) ||
      !identical(contract$assay_kind, "haemoglobin_concentration") ||
      !identical(contract$identity_kind, "pair")) {
    stop(
      "Stream metadata must declare governed pair-level haemoglobin ",
      "concentration for the requested assay",
      call. = FALSE
    )
  }
  channel_id <- .nirs_feedback_strings(contract$channel_id, "channel_id")
  if (!identical(channel_id, as.character(info@channel_names))) {
    stop("Stream NIRS identity must exactly match ordered channel names",
         call. = FALSE)
  }
  units <- info@channel_units
  if (!is.character(units) || !identical(units, rep("uM", length(channel_id)))) {
    stop("Every live HbO stream channel must declare unit `uM`",
         call. = FALSE)
  }
  source <- contract$source_index
  detector <- contract$detector_index
  if (!is.numeric(source) || is.complex(source) ||
      !is.null(dim(source)) ||
      !is.numeric(detector) || is.complex(detector) ||
      !is.null(dim(detector)) ||
      length(source) != length(channel_id) ||
      length(detector) != length(channel_id) ||
      anyNA(source) || anyNA(detector) ||
      any(!is.finite(source)) || any(!is.finite(detector)) ||
      any(source < 1 | source != floor(source)) ||
      any(detector < 1 | detector != floor(detector)) ||
      anyDuplicated(paste(source, detector, sep = ":"))) {
    stop("Live HbO source-detector identity is malformed", call. = FALSE)
  }
  list(
    info = info,
    assay_name = assay_name,
    channel_id = channel_id,
    source_index = as.integer(source),
    detector_index = as.integer(detector),
    sample_rate_hz = as.numeric(info@nominal_srate),
    source_fingerprint = .nirs_sha256(list(
      stream_name = info@name,
      stream_type = info@type,
      source_id = info@source_id,
      clock_domain = info@clock_domain,
      channel_id = channel_id,
      source_index = as.integer(source),
      detector_index = as.integer(detector),
      assay_name = assay_name,
      assay_kind = contract$assay_kind,
      identity_kind = contract$identity_kind,
      sample_rate_hz = as.numeric(info@nominal_srate)
    ))
  )
}

.nirs_feedback_regions <- function(regions, channel_id) {
  if (!is.list(regions) || is.object(regions) || !length(regions) ||
      is.null(names(regions)) || anyNA(names(regions)) ||
      any(!nzchar(names(regions))) || anyDuplicated(names(regions))) {
    stop("`regions` must be a non-empty uniquely named plain list",
         call. = FALSE)
  }
  positions <- vector("list", length(regions))
  names(positions) <- names(regions)
  used <- character()
  for (i in seq_along(regions)) {
    ids <- .nirs_feedback_strings(
      regions[[i]], paste0("regions$", names(regions)[[i]])
    )
    if (any(!ids %in% channel_id)) {
      stop("Every region member must exactly identify an HbO channel",
           call. = FALSE)
    }
    if (any(ids %in% used)) {
      stop("An HbO channel may belong to at most one region",
           call. = FALSE)
    }
    used <- c(used, ids)
    positions[[i]] <- as.integer(match(ids, channel_id))
  }
  positions
}

.nirs_feedback_contrast <- function(contrast, region_names) {
  if (is.numeric(contrast) && is.null(dim(contrast))) {
    if (is.complex(contrast) || length(contrast) != length(region_names) ||
        anyNA(contrast) || any(!is.finite(contrast)) ||
        is.null(names(contrast)) || anyNA(names(contrast)) ||
        anyDuplicated(names(contrast)) ||
        !setequal(names(contrast), region_names)) {
      stop(
        "Numeric `contrast` must be finite and name every region exactly",
        call. = FALSE
      )
    }
    contrast <- matrix(
      as.numeric(contrast[region_names]),
      nrow = length(region_names), ncol = 1L,
      dimnames = list(region_names, "feedback")
    )
  } else {
    if (!is.matrix(contrast) || !is.numeric(contrast) ||
        is.complex(contrast) || nrow(contrast) != length(region_names) ||
        ncol(contrast) < 1L || anyNA(contrast) ||
        any(!is.finite(contrast)) ||
        is.null(rownames(contrast)) || is.null(colnames(contrast)) ||
        anyNA(rownames(contrast)) || anyNA(colnames(contrast)) ||
        any(!nzchar(rownames(contrast))) ||
        any(!nzchar(colnames(contrast))) ||
        anyDuplicated(rownames(contrast)) ||
        anyDuplicated(colnames(contrast)) ||
        !setequal(rownames(contrast), region_names)) {
      stop(
        "Matrix `contrast` must have exact region rows and unique target ",
        "columns",
        call. = FALSE
      )
    }
    target_names <- colnames(contrast)
    contrast <- unname(contrast[region_names, , drop = FALSE])
    dimnames(contrast) <- list(region_names, target_names)
  }
  contrast
}

.nirs_feedback_region_values <- function(samples, positions) {
  value <- vapply(
    positions,
    function(index) rowMeans(samples[, index, drop = FALSE]),
    numeric(nrow(samples))
  )
  value <- as.matrix(value)
  if (nrow(value) != nrow(samples)) value <- t(value)
  colnames(value) <- names(positions)
  if (anyNA(value) || any(!is.finite(value))) {
    stop("Live HbO regional aggregation produced non-finite values",
         call. = FALSE)
  }
  value
}

.nirs_feedback_trim_buffer <- function(time, value, cutoff) {
  before <- which(time < cutoff)
  keep_before <- if (length(before)) {
    utils::tail(before, 1L)
  } else {
    integer()
  }
  keep <- unique(c(keep_before, which(time >= cutoff)))
  list(time = time[keep], value = value[keep, , drop = FALSE])
}

.nirs_feedback_weighted_mean <- function(time, value, end, duration) {
  start <- end - duration
  tolerance <- max(1e-9, duration * 1e-9)
  if (!length(time) || min(time) > start + tolerance ||
      max(time) < end - tolerance) {
    return(NULL)
  }
  interior <- time > start & time < end
  grid <- c(start, time[interior], end)
  output <- numeric(ncol(value))
  for (j in seq_len(ncol(value))) {
    ordinate <- stats::approx(
      time, value[, j], xout = grid, method = "linear",
      rule = 1, ties = "ordered"
    )$y
    if (anyNA(ordinate) || any(!is.finite(ordinate))) return(NULL)
    output[[j]] <- sum(
      diff(grid) *
        (utils::head(ordinate, -1L) + utils::tail(ordinate, -1L)) / 2
    ) / duration
  }
  names(output) <- colnames(value)
  output
}

.nirs_feedback_append_bounded <- function(x, value, limit) {
  x[[length(x) + 1L]] <- value
  if (length(x) > limit) x <- utils::tail(x, limit)
  x
}

.nirs_feedback_callback <- function(chunk, state, context) {
  samples <- chunk$samples
  timestamps <- chunk$timestamps
  if (!is.matrix(samples) || !is.numeric(samples) ||
      is.complex(samples) || anyNA(samples) ||
      any(!is.finite(samples)) || !is.numeric(timestamps) ||
      length(timestamps) != nrow(samples) || anyNA(timestamps) ||
      any(!is.finite(timestamps)) ||
      (length(timestamps) > 1L && any(diff(timestamps) <= 0)) ||
      !identical(colnames(samples), state$configuration$channel_id)) {
    stop("Neurofeedback received a malformed HbO chunk", call. = FALSE)
  }
  expected_step <- 1 / state$configuration$sample_rate_hz
  tolerance <- max(1e-9, 0.05 * expected_step)
  delta <- diff(c(state$last_timestamp, timestamps))
  if (length(delta) && any(abs(delta - expected_step) > tolerance)) {
    stop("Live HbO clock discontinuity exceeds the governed tolerance",
         call. = FALSE)
  }
  region <- .nirs_feedback_region_values(
    samples, state$configuration$region_positions
  )
  if (is.null(state$start_timestamp)) {
    state$start_timestamp <- timestamps[[1L]]
    state$baseline_end_timestamp <-
      state$start_timestamp + state$configuration$baseline_seconds
  }

  baseline_rows <- !state$baseline_complete &
    timestamps < state$baseline_end_timestamp
  if (any(baseline_rows)) {
    state$baseline_sum <- state$baseline_sum +
      colSums(region[baseline_rows, , drop = FALSE])
    state$baseline_count <- state$baseline_count + sum(baseline_rows)
  }
  state$buffer_time <- c(state$buffer_time, timestamps)
  state$buffer_values <- rbind(state$buffer_values, region)
  state$last_timestamp <- utils::tail(timestamps, 1L)

  if (!state$baseline_complete &&
      state$last_timestamp >= state$baseline_end_timestamp) {
    if (state$baseline_count < 2L) {
      stop("The governed baseline interval contains too few samples",
           call. = FALSE)
    }
    state$baseline <- state$baseline_sum / state$baseline_count
    names(state$baseline) <- state$configuration$region_names
    state$baseline_complete <- TRUE
    state$next_update_timestamp <- max(
      state$baseline_end_timestamp,
      state$start_timestamp + state$configuration$smoothing_seconds
    )
  }

  emitted <- 0L
  if (state$baseline_complete &&
      state$last_timestamp + tolerance >= state$next_update_timestamp) {
    due <- state$next_update_timestamp
    missed <- max(
      0L,
      as.integer(floor(
        (state$last_timestamp - due) /
          state$configuration$update_seconds + 1e-10
      ))
    )
    if (missed > 0L) {
      state$skipped_updates <- .nirs_feedback_append_bounded(
        state$skipped_updates,
        list(
          reason = "obsolete_schedule",
          first_timestamp = due,
          count = as.numeric(missed)
        ),
        .nirs_feedback_update_limit
      )
      due <- due + missed * state$configuration$update_seconds
    }
    smoothed <- .nirs_feedback_weighted_mean(
      state$buffer_time, state$buffer_values, state$last_timestamp,
      state$configuration$smoothing_seconds
    )
    if (is.null(smoothed)) {
      state$skipped_updates <- .nirs_feedback_append_bounded(
        state$skipped_updates,
        list(
          reason = "incomplete_smoothing_window",
          first_timestamp = due,
          count = 1
        ),
        .nirs_feedback_update_limit
      )
    } else {
      if (state$emitted_total >= 2^53) {
        stop("Neurofeedback update sequence exceeded exact numeric identity",
             call. = FALSE)
      }
      centred <- smoothed - state$baseline
      values <- as.numeric(
        crossprod(state$configuration$contrast, centred)
      )
      names(values) <- state$configuration$target_names
      if (anyNA(values) || any(!is.finite(values))) {
        stop("Neurofeedback target calculation produced non-finite values",
             call. = FALSE)
      }
      update <- list(
        sequence = as.numeric(state$emitted_total),
        timestamp = as.numeric(state$last_timestamp),
        values = values,
        region_values = centred,
        scheduled_timestamp = as.numeric(due)
      )
      state$updates <- .nirs_feedback_append_bounded(
        state$updates, update, .nirs_feedback_update_limit
      )
      state$emitted_total <- state$emitted_total + 1
      emitted <- 1L
    }
    state$next_update_timestamp <- due +
      state$configuration$update_seconds
  }
  cutoff <- state$last_timestamp -
    max(state$configuration$smoothing_seconds,
        2 * state$configuration$update_seconds,
        4 * expected_step)
  trimmed <- .nirs_feedback_trim_buffer(
    state$buffer_time, state$buffer_values, cutoff
  )
  state$buffer_time <- trimmed$time
  state$buffer_values <- trimmed$value
  list(
    output = chunk,
    state = state,
    events = list(),
    diagnostics = list(
      baseline_complete = state$baseline_complete,
      emitted_updates = emitted,
      emitted_total = state$emitted_total,
      skipped_total = length(state$skipped_updates)
    )
  )
}

.nirs_feedback_operation_state <- function(controller) {
  runtime <- .nirs_feedback_assert(controller)
  state <- PhysioStream::pipelineState(runtime$pipeline)
  names <- vapply(state$operations, `[[`, character(1), "name")
  index <- which(names == runtime$operation_name)
  if (length(index) != 1L) {
    stop("Neurofeedback pipeline operation identity changed",
         call. = FALSE)
  }
  state$operations[[index]]$state
}

.nirs_feedback_new_updates <- function(controller, state) {
  runtime <- .nirs_feedback_assert(controller)
  expected <- runtime$last_delivered + 1
  updates <- Filter(
    function(update) update$sequence >= expected,
    state$updates
  )
  if (!length(updates)) return(list())
  sequence <- vapply(updates, `[[`, numeric(1), "sequence")
  if (!identical(sequence, seq(expected, length.out = length(sequence)))) {
    stop(
      "Bounded neurofeedback state no longer contains a contiguous ",
      "undelivered update sequence",
      call. = FALSE
    )
  }
  updates
}

#' Construct a governed live fNIRS neurofeedback controller
#'
#' The controller builds a public PhysioStream pipeline and
#' `BiofeedbackScope`. It computes frozen-baseline regional HbO contrasts in a
#' committed pipeline operation and publishes each eligible multi-target value
#' through one atomic `biofeedbackUpdate()` call. Construction is
#' side-effect-free.
#'
#' The stream must declare `metadata$nirs` fields `assay_name`,
#' `assay_kind = "haemoglobin_concentration"`, `identity_kind = "pair"`,
#' `channel_id`, `source_index`, and `detector_index`; ordered channel names
#' and units must be the exact HbO identities and `"uM"`.
#'
#' @param stream A created or open regular-rate numeric `StreamSource`.
#' @param regions A named list of non-overlapping exact HbO channel IDs.
#' @param contrast A region-named numeric vector or region-by-target matrix.
#' @param assay_name Exact live assay name, currently `"HbO"`.
#' @param baseline_seconds Positive frozen-baseline duration.
#' @param update_seconds Positive scheduled update interval.
#' @param smoothing_seconds Positive trailing time-weighted mean duration.
#'
#' @return A side-effect-free `NIRSNeurofeedback` runtime.
#' @references Kober et al. (2017), DOI: 10.3389/fnhum.2017.00081.
#' @importFrom PhysioStream biofeedbackScope biofeedbackStart biofeedbackStep biofeedbackUpdate biofeedbackState biofeedbackStop onChunk pipelineState streamInfo streamState streamPipeline
#' @export
nirsNeurofeedback <- function(
    stream,
    regions,
    contrast,
    assay_name = "HbO",
    baseline_seconds = 20,
    update_seconds = 0.1,
    smoothing_seconds = 2) {
  contract <- .nirs_feedback_contract(stream, assay_name)
  if (!identical(assay_name, "HbO")) {
    stop("Neurofeedback currently requires exact `assay_name = \"HbO\"`",
         call. = FALSE)
  }
  region_positions <- .nirs_feedback_regions(
    regions, contract$channel_id
  )
  contrast <- .nirs_feedback_contrast(
    contrast, names(region_positions)
  )
  baseline_seconds <- .nirs_quality_number(
    baseline_seconds, "baseline_seconds",
    lower = 0, lower_open = TRUE, upper = 1800
  )
  update_seconds <- .nirs_quality_number(
    update_seconds, "update_seconds",
    lower = 1 / 120, upper = 1
  )
  smoothing_seconds <- .nirs_quality_number(
    smoothing_seconds, "smoothing_seconds",
    lower = 0, lower_open = TRUE, upper = 600
  )
  lifecycle <- PhysioStream::streamState(stream)
  source_lifecycle <- if (identical(lifecycle, "created")) {
    "own"
  } else if (identical(lifecycle, "open")) {
    "borrow"
  } else {
    stop("`stream` must be in the created or open state", call. = FALSE)
  }
  chunk_size <- max(
    1L,
    as.integer(floor(contract$sample_rate_hz * update_seconds))
  )
  pipeline <- PhysioStream::streamPipeline(
    source = stream,
    chunk_size = chunk_size,
    queue_capacity = 16L,
    backpressure = "error",
    latency_budget_ms = 1000 * update_seconds
  )
  operation_name <- "nirs_neurofeedback_hbo"
  initial_state <- list(
    configuration = list(
      schema_version = .nirs_feedback_schema,
      assay_name = contract$assay_name,
      channel_id = contract$channel_id,
      region_names = names(region_positions),
      region_positions = region_positions,
      target_names = colnames(contrast),
      contrast = unname(contrast),
      sample_rate_hz = contract$sample_rate_hz,
      baseline_seconds = baseline_seconds,
      update_seconds = update_seconds,
      smoothing_seconds = smoothing_seconds,
      source_fingerprint = contract$source_fingerprint,
      clock_domain = contract$info@clock_domain
    ),
    start_timestamp = NULL,
    baseline_end_timestamp = NULL,
    last_timestamp = NULL,
    baseline_sum = stats::setNames(
      numeric(length(region_positions)), names(region_positions)
    ),
    baseline_count = 0,
    baseline_complete = FALSE,
    baseline = NULL,
    next_update_timestamp = NULL,
    buffer_time = numeric(),
    buffer_values = matrix(
      numeric(), nrow = 0L, ncol = length(region_positions),
      dimnames = list(NULL, names(region_positions))
    ),
    emitted_total = 0,
    updates = list(),
    skipped_updates = list()
  )
  PhysioStream::onChunk(
    pipeline, .nirs_feedback_callback, state = initial_state,
    name = operation_name, kind = "feature"
  )
  derived <- lapply(colnames(contrast), function(name) {
    list(type = "external", unit = "uM")
  })
  names(derived) <- colnames(contrast)
  scope <- PhysioStream::biofeedbackScope(
    source = stream,
    pipeline = pipeline,
    channels = contract$channel_id,
    derived = derived,
    window_seconds = min(
      3600, max(10, baseline_seconds + smoothing_seconds)
    ),
    update_hz = 1 / update_seconds,
    max_points = 2000L,
    source_lifecycle = source_lifecycle,
    launch = FALSE
  )
  runtime <- new.env(parent = emptyenv())
  runtime$scope <- scope
  runtime$pipeline <- pipeline
  runtime$operation_name <- operation_name
  runtime$configuration <- list(
    schema_version = .nirs_feedback_schema,
    assay_name = contract$assay_name,
    channel_id = contract$channel_id,
    source_index = contract$source_index,
    detector_index = contract$detector_index,
    regions = lapply(regions, as.character),
    target_names = colnames(contrast),
    contrast = unname(contrast),
    baseline_seconds = baseline_seconds,
    update_seconds = update_seconds,
    smoothing_seconds = smoothing_seconds,
    source_fingerprint = contract$source_fingerprint,
    clock_provenance = list(
      clock_domain = contract$info@clock_domain,
      source_id = contract$info@source_id,
      nominal_srate = contract$sample_rate_hz
    )
  )
  runtime$lifecycle <- "created"
  runtime$last_delivered <- -1
  runtime$receipts <- list()
  runtime$busy <- FALSE
  runtime$last_error_class <- NULL
  .nirs_feedback_registry$counter <-
    .nirs_feedback_registry$counter + 1
  id <- digest::digest(
    list(runtime$configuration, .nirs_feedback_registry$counter),
    algo = "sha256", serialize = TRUE, serializeVersion = 2
  )
  assign(id, runtime, envir = .nirs_feedback_registry)
  controller <- new.env(parent = emptyenv())
  controller$id <- id
  class(controller) <- c("NIRSNeurofeedback", "NIRSNeurofeedbackRuntime")
  lockEnvironment(controller, bindings = TRUE)
  reg.finalizer(controller, function(value) {
    id <- value$id
    if (exists(id, envir = .nirs_feedback_registry, inherits = FALSE)) {
      rm(list = id, envir = .nirs_feedback_registry)
    }
  }, onexit = TRUE)
  controller
}

#' Control and inspect live NIRS neurofeedback
#'
#' `nirsNeurofeedbackStep()` processes at most `max_chunks` synchronously,
#' flushing committed target values after each chunk so the external update
#' timestamp remains the latest scope signal time.
#'
#' @param controller A `NIRSNeurofeedback` runtime.
#' @param max_chunks Positive exact per-call work bound.
#' @return Lifecycle functions return `controller` invisibly. Step and state
#'   functions return portable plain lists. `nirsNeurofeedbackScope()` returns
#'   the governed `BiofeedbackScope` for a viewer.
#' @name nirsNeurofeedback-lifecycle
NULL

#' @rdname nirsNeurofeedback-lifecycle
#' @export
nirsNeurofeedbackStart <- function(controller) {
  runtime <- .nirs_feedback_assert(controller)
  if (runtime$busy || !identical(runtime$lifecycle, "created")) {
    stop("nirsNeurofeedbackStart requires an idle created controller",
         call. = FALSE)
  }
  runtime$busy <- TRUE
  on.exit(runtime$busy <- FALSE, add = TRUE)
  tryCatch(
    PhysioStream::biofeedbackStart(runtime$scope),
    error = function(e) {
      runtime$lifecycle <- "error_stopped"
      runtime$last_error_class <- class(e)[[1L]]
      stop(e)
    }
  )
  runtime$lifecycle <- "running"
  invisible(controller)
}

#' @rdname nirsNeurofeedback-lifecycle
#' @export
nirsNeurofeedbackStep <- function(controller, max_chunks = 16L) {
  runtime <- .nirs_feedback_assert(controller)
  if (runtime$busy || !identical(runtime$lifecycle, "running")) {
    stop("nirsNeurofeedbackStep requires an idle running controller",
         call. = FALSE)
  }
  max_chunks <- .nirs_quality_number(
    max_chunks, "max_chunks", lower = 1, upper = 100000, integer = TRUE
  )
  runtime$busy <- TRUE
  on.exit(runtime$busy <- FALSE, add = TRUE)
  processed <- 0L
  delivered <- 0L
  latest_frame <- NULL
  tryCatch({
    for (i in seq_len(max_chunks)) {
      stepped <- PhysioStream::biofeedbackStep(runtime$scope, 1L)
      latest_frame <- stepped$frame_id
      if (!isTRUE(stepped$updated)) break
      processed <- processed + 1L
      pipeline_state <- PhysioStream::pipelineState(runtime$pipeline)
      if (isTRUE(utils::tail(
        pipeline_state$latency$budget_exceeded, 1L
      ))) {
        stop(
          "Neurofeedback pipeline exceeded its governed latency budget",
          call. = FALSE
        )
      }
      operation_state <- .nirs_feedback_operation_state(controller)
      updates <- .nirs_feedback_new_updates(controller, operation_state)
      for (update in updates) {
        receipt <- PhysioStream::biofeedbackUpdate(
          runtime$scope,
          update$values,
          timestamp = update$timestamp,
          names = names(update$values),
          units = rep("uM", length(update$values)),
          sequence = update$sequence
        )
        latency <- utils::tail(
          pipeline_state$latency$end_to_end_ms, 1L
        )
        runtime$receipts <- .nirs_feedback_append_bounded(
          runtime$receipts,
          list(
            receipt = receipt,
            values = update$values,
            scheduled_timestamp = update$scheduled_timestamp,
            region_values = update$region_values,
            pipeline_latency_ms = as.numeric(latency)
          ),
          .nirs_feedback_receipt_limit
        )
        runtime$last_delivered <- update$sequence
        delivered <- delivered + 1L
      }
    }
  }, error = function(e) {
    runtime$lifecycle <- "error_stopped"
    runtime$last_error_class <- class(e)[[1L]]
    stop(e)
  })
  list(
    updated = processed > 0L,
    n_processed = as.numeric(processed),
    n_delivered = as.numeric(delivered),
    last_sequence = if (runtime$last_delivered < 0) {
      NULL
    } else {
      as.numeric(runtime$last_delivered)
    },
    frame_id = latest_frame,
    schema_version = .nirs_feedback_schema
  )
}

#' @rdname nirsNeurofeedback-lifecycle
#' @export
nirsNeurofeedbackState <- function(controller) {
  runtime <- .nirs_feedback_assert(controller)
  operation <- .nirs_feedback_operation_state(controller)
  receipts <- unserialize(serialize(
    runtime$receipts, NULL, version = 3L
  ))
  list(
    schema_version = .nirs_feedback_schema,
    lifecycle = runtime$lifecycle,
    configuration = unserialize(serialize(
      runtime$configuration, NULL, version = 3L
    )),
    baseline_complete = operation$baseline_complete,
    baseline = operation$baseline,
    baseline_sample_count = as.numeric(operation$baseline_count),
    update_timestamps = vapply(
      receipts,
      function(value) value$receipt$timestamp,
      numeric(1)
    ),
    target_values = lapply(
      receipts,
      function(value) value$values
    ),
    receipts = receipts,
    skipped_updates = operation$skipped_updates,
    emitted_total = as.numeric(operation$emitted_total),
    last_successful_sequence = if (runtime$last_delivered < 0) {
      NULL
    } else {
      as.numeric(runtime$last_delivered)
    },
    last_error_class = runtime$last_error_class
  )
}

#' @rdname nirsNeurofeedback-lifecycle
#' @export
nirsNeurofeedbackStop <- function(controller) {
  runtime <- .nirs_feedback_assert(controller)
  if (runtime$busy || !identical(runtime$lifecycle, "running")) {
    stop("nirsNeurofeedbackStop requires an idle running controller",
         call. = FALSE)
  }
  runtime$busy <- TRUE
  on.exit(runtime$busy <- FALSE, add = TRUE)
  tryCatch(
    PhysioStream::biofeedbackStop(runtime$scope),
    error = function(e) {
      runtime$lifecycle <- "error_stopped"
      runtime$last_error_class <- class(e)[[1L]]
      stop(e)
    }
  )
  runtime$lifecycle <- "stopped"
  invisible(controller)
}

#' @rdname nirsNeurofeedback-lifecycle
#' @export
nirsNeurofeedbackScope <- function(controller) {
  runtime <- .nirs_feedback_assert(controller)
  runtime$scope
}

#' Display a governed NIRS neurofeedback runtime
#'
#' @param x A `NIRSNeurofeedback` runtime.
#' @param ... Unused.
#' @export
print.NIRSNeurofeedback <- function(x, ...) {
  state <- nirsNeurofeedbackState(x)
  cat(
    "NIRS neurofeedback <", state$lifecycle, ">: ",
    length(state$configuration$regions), " regions -> ",
    length(state$configuration$target_names), " targets; delivered ",
    length(state$receipts), "\n", sep = ""
  )
  invisible(x)
}
