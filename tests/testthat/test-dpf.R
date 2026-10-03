scholkmann_oracle <- function(wavelength, age) {
  223.3 + 0.05624 * age^0.8493 -
    5.723e-7 * wavelength^3 +
    0.001245 * wavelength^2 -
    0.9025 * wavelength
}

test_that("Scholkmann DPF matches the published equation", {
  wavelength <- rep(c(690, 760, 807, 832, 850), each = 4L)
  age <- rep(c(0, 20, 50, 70), times = 5L)
  expected <- scholkmann_oracle(wavelength, age)
  actual <- dpf(wavelength, age)
  expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
  expect_identical(attr(actual, "model"), "scholkmann2013")
  expect_false(any(attr(actual, "extrapolated")))
  expect_false(isTRUE(all.equal(
    actual[[1L]],
    233.3 - 5.723e-7 * wavelength[[1L]]^3 +
      0.001245 * wavelength[[1L]]^2 - 0.9025 * wavelength[[1L]]
  )))
})

test_that("Duncan DPF matches all four published age equations", {
  wavelength <- rep(c(690, 744, 807, 832), each = 3L)
  age <- rep(c(0, 20, 50), times = 4L)
  expected <- c(
    5.38 + 0.049 * age[1:3]^0.877,
    5.11 + 0.106 * age[4:6]^0.723,
    4.99 + 0.067 * age[7:9]^0.814,
    4.67 + 0.062 * age[10:12]^0.819
  )
  actual <- dpf(wavelength, age, model = "duncan1996")
  expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
})

test_that("DPF vectorization and governed domains are exact", {
  expect_length(dpf(c(760, 850), 20), 2L)
  expect_length(dpf(760, c(20, 30)), 2L)
  expect_error(dpf(c(760, 800), c(10, 20, 30)), "one or equal")
  expect_error(dpf(0, 20), "positive")
  expect_error(dpf(760, -1), "non-negative")
  expect_error(dpf(760, Inf), "finite")
  expect_error(dpf(760, 20, model = "schol"), "exactly")
  expect_error(dpf(649, 20), "650-950")
  expect_warning(out <- dpf(649, 20, extrapolate = TRUE), "extrapolated")
  expect_true(attr(out, "extrapolated"))
  expect_error(dpf(760, 71), "0-70")
  expect_error(dpf(760, 20, model = "duncan1996"), "exact wavelengths")
  expect_error(
    dpf(760, 20, model = "duncan1996", extrapolate = TRUE),
    "exact wavelengths"
  )
  expect_error(dpf(690, 51, model = "duncan1996"), "0-50")
  expect_warning(
    dpf(690, 51, model = "duncan1996", extrapolate = TRUE),
    "extrapolated"
  )
})
