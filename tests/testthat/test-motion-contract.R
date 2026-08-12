test_that("motion context enforces governed finite uniform data", {
  x <- make_motion_nirs()
  context <- PhysioNIRS:::.nirs_motion_context(x, "OD")
  expect_equal(context$fs, 10)
  expect_identical(context$identity$measurement_index, 1:4)
  expect_length(context$fingerprint, 1L)

  irregular <- x
  SummarizedExperiment::rowData(irregular)$time_seconds[10] <- 0.905
  expect_error(
    motionArtifactDetect(irregular),
    "uniformly sampled"
  )

  nonfinite <- x
  value <- motion_assay(nonfinite, "OD")
  value[4, 2] <- Inf
  SummarizedExperiment::assay(
    nonfinite, "OD", withDimnames = FALSE
  ) <- value
  expect_error(motionArtifactDetect(nonfinite), "row 4, channel 2")
  expect_error(motionArtifactDetect(matrix(1, 4, 2)), "PhysioExperiment")
  expect_error(motionArtifactDetect(x, assay_name = "missing"), "available")
})

test_that("motion output helper preserves source object and contracts", {
  x <- make_motion_nirs()
  before <- serialize(x, NULL, version = 2)
  out <- tddr(x)
  expect_identical(serialize(x, NULL, version = 2), before)
  expect_true(all(c("OD", "OD_tddr") %in%
                    SummarizedExperiment::assayNames(out)))
  expect_identical(motion_assay(out, "OD"), motion_assay(x, "OD"))
  contract <- S4Vectors::metadata(out)$nirs$assays$OD_tddr
  expect_identical(contract$source_assay, "OD")
  expect_identical(contract$method, "tddr")
  expect_true(contract$motion_corrected)
  expect_error(tddr(out), "already exists")
})

test_that("motion methods accept governed MBLL pair identity", {
  collapsed <- mbll(make_motion_nirs(), pathlength_factor = 6)
  context <- PhysioNIRS:::.nirs_motion_context(collapsed, "HbO")
  expect_identical(context$identity$measurement_id, c("S1_D1", "S2_D2"))
  expect_identical(context$identity$source_index, c(1L, 2L))
  expect_identical(context$identity$detector_index, c(1L, 2L))

  mask <- motionArtifactDetect(
    collapsed,
    assay_name = "HbO",
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 1e6,
    amplitude_threshold = 1e12
  )
  expect_identical(mask$measurement_id, c("S1_D1", "S2_D2"))
  corrected <- tddr(collapsed, assay_name = "HbO")
  expect_true(all(is.finite(motion_assay(corrected, "HbO_tddr"))))
})

test_that("motion scalar and enum inputs reject coercion and partial values", {
  x <- make_motion_nirs()
  expect_error(motionArtifactDetect(x, t_motion = matrix(1)), "t_motion")
  expect_error(motionArtifactDetect(x, t_mask = -1), "t_mask")
  expect_error(motionArtifactDetect(x, sd_threshold = 0), "sd_threshold")
  expect_error(
    motionArtifactDetect(x, amplitude_threshold = Inf),
    "amplitude_threshold"
  )
  expect_error(
    waveletMotionCorrect(x, wavelet = "db"),
    "wavelet"
  )
  expect_error(tddr(x, max_iter = 2.5), "max_iter")
  expect_error(splineMotionCorrect(x, p = 1.1), "`p`")
})
