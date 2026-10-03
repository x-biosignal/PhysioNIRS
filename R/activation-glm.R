# First-level cortical-activation GLM for fNIRS. Task regressors are boxcars over
# each condition's stimulus blocks convolved with a canonical (double-gamma) HRF;
# a discrete-cosine drift basis high-passes the slow trend and optional
# short-separation nuisance regressors (see shortSeparationDesign()) absorb
# systemic physiology. Every haemoglobin channel is fit at once by OLS, with
# optional AR(1) prewhitening because fNIRS residuals are strongly autocorrelated
# (naive OLS t-statistics are otherwise inflated). Mirrors the PsPM-style GLM in
# PhysioEDA (edaGLM), specialised to the HRF, block designs and many channels.

# SPM canonical double-gamma HRF sampled on the signal grid (unit peak). deriv=1
# returns the temporal-derivative basis.
.nirs_hrf_kernel <- function(fs, peak = 6, undershoot = 16, ratio = 6,
                             disp1 = 1, disp2 = 1, deriv = 0L, len_sec = 32) {
  dt <- 1 / fs
  t <- seq(0, len_sec, by = dt)
  canon <- function(pk) stats::dgamma(t, shape = pk / disp1, scale = disp1) -
    stats::dgamma(t, shape = undershoot / disp2, scale = disp2) / ratio
  h <- if (deriv == 1L) canon(peak) - canon(peak + 1) else canon(peak)
  mx <- max(abs(h)); if (mx > 0) h <- h / mx
  h
}

# Boxcar over each [onset, onset+duration] block, summed, convolved (causal) with
# the HRF kernel -> a length-n condition regressor.
.nirs_condition_regressor <- function(onsets, durations, n, fs, kern) {
  box <- numeric(n)
  for (i in seq_along(onsets)) {
    s <- as.integer(round(onsets[i] * fs)) + 1L
    e <- as.integer(round((onsets[i] + durations[i]) * fs))
    s <- max(1L, s); e <- min(n, e)
    if (e >= s) box[s:e] <- box[s:e] + 1
  }
  reg <- numeric(n + length(kern))
  kl <- length(kern)
  for (s in which(box != 0)) reg[s:(s + kl - 1L)] <- reg[s:(s + kl - 1L)] + box[s] * kern
  reg[seq_len(n)]
}

# Discrete-cosine (high-pass) drift basis: intercept + k low-frequency cosines.
.nirs_drift_basis <- function(n, k = 6) {
  tt <- seq(0, 1, length.out = n)
  B <- matrix(1, n, 1)
  if (k > 0) for (j in seq_len(k)) B <- cbind(B, cos(pi * j * tt))
  colnames(B) <- c("drift0", if (k > 0) paste0("drift", seq_len(k)))
  B
}

# SVD pseudo-inverse with a rank tolerance (no MASS dependency).
.nirs_pinv <- function(M) {
  sv <- svd(M)
  tol <- max(dim(M)) * .Machine$double.eps * max(sv$d)
  pos <- sv$d > tol
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  sv$v[, pos, drop = FALSE] %*% (t(sv$u[, pos, drop = FALSE]) / sv$d[pos])
}

# Core multivariate fit: OLS, or per-channel AR(1) prewhitened GLS.
.nirs_glm_fit <- function(Y, X, prewhiten = "none") {
  n <- nrow(X); p <- ncol(X); nchan <- ncol(Y)
  rank <- qr(X)$rank
  XtXinv <- .nirs_pinv(crossprod(X))
  B_ols <- XtXinv %*% crossprod(X, Y)
  beta <- matrix(NA_real_, p, nchan, dimnames = list(colnames(X), colnames(Y)))
  se <- matrix(NA_real_, p, nchan)
  sigma2 <- rep(NA_real_, nchan); rho <- rep(NA_real_, nchan)
  xtxinv_list <- vector("list", nchan); xtxinv_shared <- NULL
  resid <- Y - X %*% B_ols
  if (prewhiten == "none") {
    dfres <- n - rank
    beta <- B_ols
    sigma2 <- colSums(resid^2) / dfres
    dvar <- diag(XtXinv)
    for (ch in seq_len(nchan)) se[, ch] <- sqrt(sigma2[ch] * dvar)
    xtxinv_shared <- XtXinv
  } else {
    dfres <- (n - 1L) - rank
    for (ch in seq_len(nchan)) {
      r <- resid[, ch]
      rho_ch <- sum(r[-1] * r[-n]) / sum(r[-n]^2)
      rho_ch <- max(min(rho_ch, 0.99), -0.99); rho[ch] <- rho_ch
      Xw <- X[-1, , drop = FALSE] - rho_ch * X[-n, , drop = FALSE]
      yw <- Y[-1, ch] - rho_ch * Y[-n, ch]
      inv <- .nirs_pinv(crossprod(Xw))
      b <- inv %*% crossprod(Xw, yw)
      beta[, ch] <- b
      s2 <- sum((yw - Xw %*% b)^2) / dfres
      sigma2[ch] <- s2; se[, ch] <- sqrt(s2 * diag(inv))
      xtxinv_list[[ch]] <- inv
    }
  }
  list(beta = beta, se = se, sigma2 = sigma2, df = dfres, rank = rank,
       collinear = rank < p, rho = rho, prewhiten = prewhiten,
       xtxinv_shared = xtxinv_shared, xtxinv_list = xtxinv_list)
}

#' First-level fNIRS cortical-activation GLM
#'
#' Fits a general linear model of every haemoglobin channel against
#' HRF-convolved task regressors: for each condition, a boxcar over its stimulus
#' blocks convolved with a canonical double-gamma HRF, plus a discrete-cosine
#' drift basis and optional short-separation nuisance regressors. Returns the
#' per-channel activation amplitude (beta), standard error and t-statistic for
#' each condition. Run it separately on the `HbO` and `HbR` assays.
#'
#' @param x A `PhysioExperiment` of haemoglobin concentration (e.g. the `HbO`
#'   assay from [mbll()]).
#' @param events A `data.frame` with numeric `onset` and `duration` (seconds) and
#'   a `condition` column; if `NULL`, the object's events ([PhysioExperiment::getEvents]
#'   `onset`/`duration`/`value`) are used, with `value` as the condition.
#' @param assay_name Assay to model (default `"HbO"`).
#' @param conditions Optional character vector selecting/ordering conditions.
#' @param basis `"canonical"` (HRF only) or `"derivative"` (HRF + temporal
#'   derivative per condition).
#' @param drift Number of discrete-cosine drift terms (default 6).
#' @param nuisance Optional numeric matrix of extra nuisance regressors (e.g.
#'   [shortSeparationDesign()] output).
#' @param hrf List of HRF parameters (`peak`, `undershoot`, `ratio`).
#' @param prewhiten `"none"` or `"ar1"` (per-channel AR(1) prewhitening; strongly
#'   recommended for fNIRS to avoid inflated t-statistics).
#' @return An S3 `nirs_activation_glm` with `coefficients` (a `channel` x
#'   `condition` data.frame of `beta`, `se`, `t`, `df`, `p`), the `beta`/`se`
#'   matrices, per-channel `sigma2`/`rho`, and the `design` (for
#'   [nirsActivationContrast()]).
#' @references Friston KJ et al. (1994) Statistical parametric maps.
#'   *Hum Brain Mapp* 2:189-210. Huppert TJ et al. (2009) HomER. *Appl Opt*
#'   48:D280.
#' @seealso [nirsActivationContrast()], [shortSeparationDesign()], [mbll()]
#' @examples
#' hb <- mbll(intensityToOD(PhysioNIRS:::.nirs_demo_object()))
#' fit <- nirsActivationGLM(hb, assay_name = "HbO")
#' fit$coefficients
#' @export
nirsActivationGLM <- function(x, events = NULL, assay_name = "HbO",
                              conditions = NULL,
                              basis = c("canonical", "derivative"),
                              drift = 6, nuisance = NULL,
                              hrf = list(peak = 6, undershoot = 16, ratio = 6),
                              prewhiten = c("none", "ar1")) {
  stopifnot(inherits(x, "PhysioExperiment"))
  basis <- match.arg(basis); prewhiten <- match.arg(prewhiten)
  Y <- as.matrix(SummarizedExperiment::assay(x, assay_name))
  n <- nrow(Y); fs <- PhysioExperiment::samplingRate(x)
  if (is.null(events)) {
    events <- methods::slot(PhysioExperiment::getEvents(x), "events")
  }
  events <- as.data.frame(events)
  if (!all(c("onset", "duration") %in% names(events))) {
    stop("'events' needs 'onset' and 'duration' columns.", call. = FALSE)
  }
  if (!"condition" %in% names(events)) {
    events$condition <- if ("value" %in% names(events)) {
      as.character(events$value)
    } else "task"
  }
  if (is.null(conditions)) conditions <- unique(as.character(events$condition))
  hp <- utils::modifyList(list(peak = 6, undershoot = 16, ratio = 6), hrf)
  reg <- function(kern) vapply(conditions, function(cc) {
    sel <- as.character(events$condition) == cc
    .nirs_condition_regressor(events$onset[sel], events$duration[sel], n, fs, kern)
  }, numeric(n))
  Xc <- reg(.nirs_hrf_kernel(fs, hp$peak, hp$undershoot, hp$ratio, deriv = 0L))
  colnames(Xc) <- conditions
  Xd <- NULL
  if (basis == "derivative") {
    Xd <- reg(.nirs_hrf_kernel(fs, hp$peak, hp$undershoot, hp$ratio, deriv = 1L))
    colnames(Xd) <- paste0(conditions, "_d")
  }
  design <- cbind(Xc, Xd, .nirs_drift_basis(n, drift))
  if (!is.null(nuisance)) {
    nuisance <- as.matrix(nuisance)
    if (is.null(colnames(nuisance))) {
      colnames(nuisance) <- paste0("nuis", seq_len(ncol(nuisance)))
    }
    design <- cbind(design, nuisance)
  }
  chan <- colnames(Y)
  if (is.null(chan)) {
    cd <- SummarizedExperiment::colData(x)
    chan <- if (!is.null(cd$label)) as.character(cd$label) else
      paste0("ch", seq_len(ncol(Y)))
  }
  colnames(Y) <- chan
  cond_idx <- seq_along(conditions)
  fit <- .nirs_glm_fit(Y, design, prewhiten)
  if (fit$collinear) {
    warning(sprintf("nirsActivationGLM: design is rank-deficient (rank %d of %d).",
                    fit$rank, ncol(design)), call. = FALSE)
  }
  co <- do.call(rbind, lapply(seq_along(chan), function(ch) {
    b <- fit$beta[cond_idx, ch]; s <- fit$se[cond_idx, ch]
    data.frame(channel = chan[ch], condition = conditions, beta = b, se = s,
               t = b / s, df = fit$df, row.names = NULL, stringsAsFactors = FALSE)
  }))
  co$p <- 2 * stats::pt(-abs(co$t), co$df)
  structure(list(
    coefficients = co, beta = fit$beta, se = fit$se, sigma2 = fit$sigma2,
    df = fit$df, rank = fit$rank, collinear = fit$collinear, rho = fit$rho,
    prewhiten = prewhiten, basis = basis, design = design,
    conditions = conditions, cond_idx = cond_idx, channels = chan,
    assay = assay_name, xtxinv_shared = fit$xtxinv_shared,
    xtxinv_list = fit$xtxinv_list), class = "nirs_activation_glm")
}

#' Contrast of an fNIRS activation GLM
#'
#' Applies a linear contrast over the condition betas of a
#' [nirsActivationGLM()] fit, per channel: the estimate `c' beta`, its standard
#' error `sqrt(sigma^2 * c' (X'X)^-1 c)`, t-statistic, degrees of freedom and
#' p-value.
#'
#' @param fit A `nirs_activation_glm`.
#' @param contrast A numeric contrast over conditions: named by condition, or an
#'   unnamed vector with one weight per condition (e.g. `c(task = 1, rest = -1)`).
#' @return A per-channel `data.frame`: `channel`, `estimate`, `se`, `t`, `df`,
#'   `p`.
#' @seealso [nirsActivationGLM()]
#' @examples
#' hb <- mbll(intensityToOD(PhysioNIRS:::.nirs_demo_object()))
#' fit <- nirsActivationGLM(hb, assay_name = "HbO")
#' nirsActivationContrast(fit, c(task = 1))
#' @export
nirsActivationContrast <- function(fit, contrast) {
  stopifnot(inherits(fit, "nirs_activation_glm"))
  conds <- fit$conditions; ncond <- length(conds)
  cvec <- stats::setNames(numeric(ncond), conds)
  if (!is.null(names(contrast))) {
    unknown <- setdiff(names(contrast), conds)
    if (length(unknown)) {
      stop("unknown condition(s) in contrast: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    cvec[names(contrast)] <- contrast
  } else {
    if (length(contrast) != ncond) {
      stop("an unnamed contrast needs one weight per condition.", call. = FALSE)
    }
    cvec[] <- contrast
  }
  ci <- fit$cond_idx
  res <- do.call(rbind, lapply(seq_along(fit$channels), function(ch) {
    est <- sum(cvec * fit$beta[ci, ch])
    xtxinv <- if (fit$prewhiten == "ar1") fit$xtxinv_list[[ch]] else fit$xtxinv_shared
    v <- fit$sigma2[ch] * as.numeric(t(cvec) %*% xtxinv[ci, ci, drop = FALSE] %*% cvec)
    se <- sqrt(v)
    data.frame(channel = fit$channels[ch], estimate = est, se = se,
               t = est / se, df = fit$df, row.names = NULL,
               stringsAsFactors = FALSE)
  }))
  res$p <- 2 * stats::pt(-abs(res$t), res$df)
  res
}

#' @export
print.nirs_activation_glm <- function(x, ...) {
  cat(sprintf("fNIRS activation GLM (%s, %d channels, %d condition(s): %s)\n",
              x$assay, length(x$channels), length(x$conditions),
              paste(x$conditions, collapse = ", ")))
  cat(sprintf("  prewhiten = %s%s | drift + HRF%s design, rank %d, df %d\n",
              x$prewhiten,
              if (x$prewhiten == "ar1")
                sprintf(" (mean rho %.2f)", mean(x$rho, na.rm = TRUE)) else "",
              if (x$basis == "derivative") " + derivative" else "",
              x$rank, x$df))
  invisible(x)
}
