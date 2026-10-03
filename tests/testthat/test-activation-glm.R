# First-level fNIRS activation GLM.

# a minimal HbO PhysioExperiment with block events
.make_glm_pe <- function(Y, fs, onset, duration, value) {
  x <- PhysioCore::PhysioExperiment(
    assays = list(HbO = Y), samplingRate = fs,
    colData = S4Vectors::DataFrame(label = paste0("ch", seq_len(ncol(Y)))))
  PhysioCore::setEvents(
    x, PhysioCore::PhysioEvents(onset = onset, duration = duration,
                                type = "task", value = value))
}

test_that("nirsActivationGLM recovers planted activation betas", {
  fs <- 10; n <- 1400
  onsetA <- c(20, 60, 100); onsetB <- c(40, 80, 120); dur <- 12
  kern <- PhysioNIRS:::.nirs_hrf_kernel(fs, 6, 16, 6, deriv = 0L)
  rA <- PhysioNIRS:::.nirs_condition_regressor(onsetA, rep(dur, 3), n, fs, kern)
  rB <- PhysioNIRS:::.nirs_condition_regressor(onsetB, rep(dur, 3), n, fs, kern)
  set.seed(1)
  bA <- c(2.0, 0.4); bB <- c(0.3, 1.6)               # 2 channels
  Y <- cbind(bA[1] * rA + bB[1] * rB, bA[2] * rA + bB[2] * rB) +
    matrix(stats::rnorm(n * 2, 0, 0.03), n, 2)
  x <- .make_glm_pe(Y, fs, c(onsetA, onsetB), dur,
                    c("A", "A", "A", "B", "B", "B"))

  fit <- nirsActivationGLM(x, assay_name = "HbO", drift = 3)
  expect_s3_class(fit, "nirs_activation_glm")
  co <- fit$coefficients
  getb <- function(ch, cond) co$beta[co$channel == ch & co$condition == cond]
  expect_equal(getb("ch1", "A"), 2.0, tolerance = 0.1)
  expect_equal(getb("ch2", "B"), 1.6, tolerance = 0.1)
  expect_equal(getb("ch1", "B"), 0.3, tolerance = 0.1)
  # a real activation is highly significant at this SNR
  expect_lt(co$p[co$channel == "ch1" & co$condition == "A"], 1e-6)
})

test_that("nirsActivationContrast estimates A - B per channel", {
  fs <- 10; n <- 1200
  onsetA <- c(20, 60, 100); onsetB <- c(40, 80); dur <- 12
  kern <- PhysioNIRS:::.nirs_hrf_kernel(fs, 6, 16, 6)
  rA <- PhysioNIRS:::.nirs_condition_regressor(onsetA, rep(dur, 3), n, fs, kern)
  rB <- PhysioNIRS:::.nirs_condition_regressor(onsetB, rep(dur, 2), n, fs, kern)
  set.seed(2)
  Y <- cbind(2 * rA + 0.5 * rB, 0.4 * rA + 1.5 * rB) +
    matrix(stats::rnorm(n * 2, 0, 0.03), n, 2)
  x <- .make_glm_pe(Y, fs, c(onsetA, onsetB), dur,
                    c("A", "A", "A", "B", "B"))
  fit <- nirsActivationGLM(x, drift = 3)
  ctr <- nirsActivationContrast(fit, c(A = 1, B = -1))
  expect_equal(ctr$estimate[ctr$channel == "ch1"], 2 - 0.5, tolerance = 0.1)
  expect_equal(ctr$estimate[ctr$channel == "ch2"], 0.4 - 1.5, tolerance = 0.1)
  expect_true(all(ctr$se > 0))
  expect_error(nirsActivationContrast(fit, c(Z = 1)), "unknown condition")
})

test_that("AR(1) prewhitening estimates rho and still recovers betas", {
  fs <- 10; n <- 1600; dur <- 12
  onsetA <- c(20, 60, 100, 140)
  kern <- PhysioNIRS:::.nirs_hrf_kernel(fs, 6, 16, 6)
  rA <- PhysioNIRS:::.nirs_condition_regressor(onsetA, rep(dur, 4), n, fs, kern)
  set.seed(3)
  # AR(1) coloured noise, rho = 0.6
  e <- stats::filter(stats::rnorm(n, 0, 0.05), 0.6, method = "recursive")
  Y <- matrix(1.5 * rA + as.numeric(e), n, 1)
  x <- .make_glm_pe(Y, fs, onsetA, dur, rep("A", 4))
  fit <- nirsActivationGLM(x, drift = 3, prewhiten = "ar1")
  expect_equal(fit$prewhiten, "ar1")
  expect_gt(fit$rho[1], 0.3)                                    # detects the AR(1)
  expect_equal(fit$coefficients$beta[1], 1.5, tolerance = 0.2)  # beta still recovered
})
