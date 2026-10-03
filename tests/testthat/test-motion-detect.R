test_that("motion detector returns a typed empty mask", {
  x <- make_motion_nirs(values = matrix(0, 80, 4))
  mask <- motionArtifactDetect(
    x,
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 20,
    amplitude_threshold = 0.5
  )
  expect_s3_class(mask, "nirs_motion_mask")
  expect_identical(dim(mask$sample_by_measurement), c(80L, 4L))
  expect_false(any(mask$sample_by_measurement))
  expect_false(any(mask$global))
  expect_identical(
    names(mask$intervals),
    c(
      "measurement_index", "start_index", "end_index",
      "start_time", "end_time"
    )
  )
  expect_identical(nrow(mask$intervals), 0L)
  expect_identical(typeof(mask$intervals$measurement_index), "integer")
  expect_identical(
    mask$fingerprint,
    PhysioNIRS:::.nirs_motion_mask_fingerprint(mask)
  )
})

test_that("detector aligns diff to the later sample and groups wavelengths", {
  value <- matrix(0, 80, 4)
  value[40, 1] <- 1
  x <- make_motion_nirs(values = value)
  grouped <- motionArtifactDetect(
    x,
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 1e6,
    amplitude_threshold = 0.5,
    group_wavelengths = TRUE
  )
  expect_identical(which(grouped$sample_by_measurement[, 1]), 40:41)
  expect_identical(which(grouped$sample_by_measurement[, 2]), 40:41)
  expect_false(any(grouped$sample_by_measurement[, 3:4]))
  expect_identical(grouped$intervals$start_index, c(40L, 40L))
  expect_identical(grouped$intervals$end_index, c(41L, 41L))

  separate <- motionArtifactDetect(
    x,
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 1e6,
    amplitude_threshold = 0.5,
    group_wavelengths = FALSE
  )
  expect_false(any(separate$sample_by_measurement[, 2]))

  expanded <- motionArtifactDetect(
    x,
    t_motion = 0.1,
    t_mask = 0.2,
    sd_threshold = 1e6,
    amplitude_threshold = 0.5,
    group_wavelengths = FALSE
  )
  expect_identical(which(expanded$sample_by_measurement[, 1]), 38:43)
})

test_that("detector uses strict thresholds and clips boundaries", {
  equal <- matrix(0, 40, 4)
  equal[20:40, 1] <- 0.5
  x <- make_motion_nirs(values = equal)
  mask <- motionArtifactDetect(
    x,
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 1e6,
    amplitude_threshold = 0.5,
    group_wavelengths = FALSE
  )
  expect_false(any(mask$sample_by_measurement))

  edge <- matrix(0, 40, 4)
  edge[1, 1] <- 1
  edge[40, 3] <- 1
  edge_x <- make_motion_nirs(values = edge)
  edge_mask <- motionArtifactDetect(
    edge_x,
    t_motion = 0.1,
    t_mask = 1,
    sd_threshold = 1e6,
    amplitude_threshold = 0.5,
    group_wavelengths = FALSE
  )
  expect_true(all(which(edge_mask$sample_by_measurement[, 1]) >= 2L))
  expect_true(max(which(edge_mask$sample_by_measurement[, 1])) <= 12L)
  expect_true(all(which(edge_mask$sample_by_measurement[, 3]) <= 40L))
})

test_that("mask validation rejects stale and corrupted masks", {
  x <- make_motion_nirs()
  mask <- motionArtifactDetect(x)
  expect_silent(PhysioNIRS:::.nirs_validate_motion_mask(
    mask,
    PhysioNIRS:::.nirs_motion_context(x, "OD")
  ))

  changed <- x
  value <- motion_assay(changed, "OD")
  value[1, 1] <- value[1, 1] + 1
  SummarizedExperiment::assay(
    changed, "OD", withDimnames = FALSE
  ) <- value
  expect_error(
    splineMotionCorrect(changed, mask = mask),
    "stale"
  )

  corrupt <- mask
  corrupt$sample_by_measurement[1, 1] <-
    !corrupt$sample_by_measurement[1, 1]
  expect_error(splineMotionCorrect(x, mask = corrupt), "global|intervals")
  expect_error(
    splineMotionCorrect(x, mask = matrix(FALSE, nrow(x), ncol(x))),
    "nirs_motion_mask"
  )
})
