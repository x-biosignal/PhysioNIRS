test_that("csaps values match the fixed independent implementation", {
  x <- c(1, 2, 4, 6)
  y <- c(2, 4, 5, 7)
  expect_equal(
    PhysioNIRS:::.nirs_csaps_values(x, y, 0.5),
    c(
      2.3183673469387753, 3.4857142857142858,
      5.2326530612244895, 6.9632653061224490
    ),
    tolerance = 1e-14
  )
  expect_identical(
    PhysioNIRS:::.nirs_csaps_values(x, y, 1),
    y
  )
})

test_that("spline is exact identity for a clean mask", {
  x <- make_motion_nirs()
  mask <- make_manual_motion_mask(
    x, matrix(FALSE, nrow(x), ncol(x))
  )
  out <- splineMotionCorrect(x, mask = mask)
  expect_identical(
    motion_assay(out, "OD_spline"),
    motion_assay(x, "OD")
  )
})

test_that("spline handles internal, beginning, ending, and singleton intervals", {
  time <- seq(0, 19.9, by = 0.1)
  clean <- 0.02 * time + sin(2 * pi * 0.08 * time)
  corrupted <- clean
  corrupted[1:5] <- corrupted[1:5] + 2
  corrupted[80] <- corrupted[80] + 3
  corrupted[120:130] <- corrupted[120:130] - 2
  corrupted[195:200] <- corrupted[195:200] + 1.5
  x <- make_motion_nirs(values = cbind(
    corrupted, corrupted, clean, clean
  ))
  bad <- matrix(FALSE, nrow(x), ncol(x))
  bad[c(1:5, 80, 120:130, 195:200), 1:2] <- TRUE
  mask <- make_manual_motion_mask(x, bad)
  out <- splineMotionCorrect(x, mask = mask, p = 0.99)
  corrected <- motion_assay(out, "OD_spline")
  expect_true(all(is.finite(corrected)))
  expect_lt(
    sqrt(mean((corrected[, 1] - clean)^2)),
    sqrt(mean((corrupted - clean)^2))
  )
  expect_identical(corrected[, 3:4], motion_assay(x, "OD")[, 3:4])
})

test_that("spline rejects an all-artifact channel and malformed masks", {
  x <- make_motion_nirs()
  bad <- matrix(FALSE, nrow(x), ncol(x))
  bad[, 1] <- TRUE
  mask <- make_manual_motion_mask(x, bad)
  expect_error(
    splineMotionCorrect(x, mask = mask),
    "at least one clean sample"
  )
  expect_error(
    splineMotionCorrect(x, mask = list()),
    "nirs_motion_mask"
  )
})

test_that("spline can run detection once and preserves source metadata", {
  x <- make_motion_nirs()
  before_events <- PhysioCore::getEvents(x)
  before_identity <- S4Vectors::metadata(x)$snirf$measurement_list
  out <- splineMotionCorrect(
    x,
    t_motion = 0.1,
    t_mask = 0,
    sd_threshold = 1e6,
    amplitude_threshold = 10
  )
  expect_identical(PhysioCore::getEvents(out), before_events)
  expect_identical(
    S4Vectors::metadata(out)$snirf$measurement_list,
    before_identity
  )
})
