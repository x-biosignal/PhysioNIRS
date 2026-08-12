test_that("neurofeedback construction is side-effect-free and governed", {
  source <- make_feedback_source()
  controller <- make_feedback_controller(source)

  expect_s3_class(controller, "NIRSNeurofeedback")
  expect_identical(PhysioStream::streamState(source), "open")
  state <- nirsNeurofeedbackState(controller)
  expect_identical(state$lifecycle, "created")
  expect_false(state$baseline_complete)
  expect_length(state$receipts, 0L)
  expect_s3_class(
    nirsNeurofeedbackScope(controller), "BiofeedbackScope"
  )
  expect_error(controller$last_delivered <- 100, "locked environment")
  expect_error(controller$scope <- NULL, "locked environment")
  expect_match(capture.output(print(controller)), "delivered 0")
})

test_that("live HbO loop freezes baseline and emits each target once", {
  fixture <- run_feedback_fixture()
  state <- fixture$state

  expect_true(state$baseline_complete)
  expect_equal(state$baseline, c(left = 1, right = 1), tolerance = 1e-12)
  expect_identical(state$baseline_sample_count, 10)
  expect_length(state$receipts, 10L)
  expect_identical(
    vapply(state$receipts, function(x) x$receipt$sequence, numeric(1)),
    as.numeric(0:9)
  )
  expect_true(all(diff(state$update_timestamps) > 0))
  expect_true(min(state$update_timestamps) >= 1)
  expect_identical(state$last_successful_sequence, 9)
  expect_true(all(vapply(
    state$receipts,
    function(x) identical(x$receipt$names, "feedback"),
    logical(1)
  )))
  expect_lt(length(serialize(state, NULL, version = 3L)), 1024^2)
  expect_false(any(vapply(
    state,
    function(value) is.matrix(value) && nrow(value) == length(fixture$time),
    logical(1)
  )))

  expected_last <- mean(c(2.7, 2.8, 2.9, 3.0)) * 0.75
  expect_equal(
    tail(state$target_values, 1L)[[1L]][["feedback"]],
    expected_last,
    tolerance = 0.08
  )
  nirsNeurofeedbackStop(fixture$controller)
  expect_identical(
    nirsNeurofeedbackState(fixture$controller)$lifecycle, "stopped"
  )
})

test_that("multi-target updates are atomic and preserve exact names", {
  contrast <- cbind(
    activation = c(left = 1, right = -1),
    global = c(left = 0.5, right = 0.5)
  )
  source <- make_feedback_source()
  controller <- make_feedback_controller(source, contrast = contrast)
  nirsNeurofeedbackStart(controller)
  time <- 0:19 / 10
  PhysioStream::loopbackFeed(source, feedback_values(time), time)
  repeat {
    result <- nirsNeurofeedbackStep(controller, 2L)
    if (!result$updated) break
  }
  state <- nirsNeurofeedbackState(controller)

  expect_gt(length(state$receipts), 0L)
  expect_true(all(vapply(
    state$receipts,
    function(x) identical(
      x$receipt$names, c("activation", "global")
    ),
    logical(1)
  )))
  expect_true(all(vapply(
    state$target_values, function(x) length(x) == 2L, logical(1)
  )))
})

test_that("neurofeedback validates source, region, and contrast identity", {
  source <- make_feedback_source()
  expect_error(
    nirsNeurofeedback(
      source,
      list(left = "missing", right = "S3_D3"),
      c(left = 1, right = -1),
      baseline_seconds = 1
    ),
    "exactly identify"
  )
  expect_error(
    nirsNeurofeedback(
      source,
      list(left = c("S1_D1", "S2_D2"), right = "S2_D2"),
      c(left = 1, right = -1),
      baseline_seconds = 1
    ),
    "at most one region"
  )
  expect_error(
    make_feedback_controller(source, contrast = c(right = -1, other = 1)),
    "every region exactly"
  )
  expect_error(
    make_feedback_controller(
      source,
      contrast = matrix(
        1:4, 2L, dimnames = list(c("left", "other"), c("a", "b"))
      )
    ),
    "exact region rows"
  )
  expect_error(
    make_feedback_controller(source, baseline_seconds = 0)
  )
  expect_error(
    make_feedback_controller(source, update_seconds = 0.001)
  )
  expect_error(
    make_feedback_controller(source, smoothing_seconds = Inf)
  )
})

test_that("neurofeedback rejects non-HbO and malformed metadata", {
  no_metadata <- make_feedback_source(metadata = list())
  expect_error(
    make_feedback_controller(no_metadata), "must declare governed"
  )
  wrong_units <- make_feedback_source(units = rep("mM", 3L))
  expect_error(
    make_feedback_controller(wrong_units), "unit `uM`"
  )
  metadata <- list(nirs = list(
    assay_name = "HbO",
    assay_kind = "optical_density",
    identity_kind = "pair",
    channel_id = c("S1_D1", "S2_D2", "S3_D3"),
    source_index = 1:3,
    detector_index = 1:3
  ))
  wrong_kind <- make_feedback_source(metadata = metadata)
  expect_error(
    make_feedback_controller(wrong_kind), "haemoglobin concentration"
  )

  metadata$nirs$assay_kind <- "haemoglobin_concentration"
  metadata$nirs$source_index <- matrix(1:3, nrow = 1L)
  shaped_identity <- make_feedback_source(metadata = metadata)
  expect_error(
    make_feedback_controller(shaped_identity), "identity is malformed"
  )
})

test_that("clock discontinuity stops processing without false delivery", {
  source <- make_feedback_source()
  controller <- make_feedback_controller(source)
  nirsNeurofeedbackStart(controller)
  time <- c(0, 0.1, 0.2, 0.5)
  PhysioStream::loopbackFeed(source, feedback_values(time), time)

  expect_error(
    nirsNeurofeedbackStep(controller, 4L), "operation.*failed"
  )
  state <- nirsNeurofeedbackState(controller)
  expect_identical(state$lifecycle, "error_stopped")
  expect_length(state$receipts, 0L)
  expect_null(state$last_successful_sequence)
})

test_that("external sink sequence interference stops without false claim", {
  source <- make_feedback_source()
  controller <- make_feedback_controller(source)
  nirsNeurofeedbackStart(controller)
  PhysioStream::biofeedbackUpdate(
    nirsNeurofeedbackScope(controller),
    values = c(feedback = 0),
    timestamp = 0,
    units = "uM",
    sequence = 0
  )
  time <- 0:19 / 10
  PhysioStream::loopbackFeed(source, feedback_values(time), time)

  expect_error(
    nirsNeurofeedbackStep(controller, 10L), "contiguous from zero"
  )
  state <- nirsNeurofeedbackState(controller)
  expect_identical(state$lifecycle, "error_stopped")
  expect_length(state$receipts, 0L)
  expect_null(state$last_successful_sequence)
})

test_that("controller lifecycle rejects re-entry and wrong states", {
  source <- make_feedback_source()
  controller <- make_feedback_controller(source)
  expect_error(nirsNeurofeedbackStep(controller))
  expect_error(nirsNeurofeedbackStop(controller))
  nirsNeurofeedbackStart(controller)
  expect_error(nirsNeurofeedbackStart(controller))
  nirsNeurofeedbackStop(controller)
  expect_error(nirsNeurofeedbackStep(controller))
  expect_error(nirsNeurofeedbackStop(controller))
})

test_that("created sources use owned lifecycle", {
  source <- make_feedback_source(open = FALSE)
  controller <- make_feedback_controller(source)
  nirsNeurofeedbackStart(controller)
  expect_identical(
    PhysioStream::biofeedbackState(
      nirsNeurofeedbackScope(controller)
    )$lifecycle,
    "running"
  )
  nirsNeurofeedbackStop(controller)
  expect_identical(
    PhysioStream::biofeedbackState(
      nirsNeurofeedbackScope(controller)
    )$lifecycle,
    "stopped"
  )
})
