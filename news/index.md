# Changelog

## PhysioNIRS 0.6.0

First-level cortical-activation GLM.

- [`nirsActivationGLM()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationGLM.md):
  fits every haemoglobin channel against HRF-convolved task regressors
  (per-condition boxcars over the stimulus blocks convolved with a
  canonical double-gamma HRF), a discrete-cosine drift basis, and
  optional short-separation nuisance regressors, returning
  per-channel/per-condition activation beta, SE and t. Supports a
  temporal-derivative basis and per-channel AR(1) prewhitening (fNIRS
  residuals are strongly autocorrelated, so naive OLS t-statistics are
  inflated).
- [`nirsActivationContrast()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationContrast.md):
  linear contrasts over the condition betas (estimate, SE, t, p per
  channel).

## PhysioNIRS 0.5.1

- Fixed
  [`pruneChannels()`](https://x-biosignal.github.io/PhysioNIRS/reference/pruneChannels.md)
  / quality-contract validation on R \>= 4.4.0: R 4.4 changed
  `is.atomic(NULL)` from TRUE to FALSE, so `.nirs_quality_plain()`
  wrongly rejected governed quality results carrying NULL
  optional-parameter sentinels (e.g. `step_seconds = NULL`) as ‘stale or
  malformed’. NULL is now handled explicitly as a plain leaf (its
  pre-4.4 behaviour).
- Made the MBLL reference-fixture self-consistency check portable: the
  stored intensity is recomputed from optical density via exp() (a
  transcendental whose last bit is toolchain-dependent), so the exact
  (tolerance = 0) compare is now tolerance = 1e-12 (~2300x above the
  measured cross-platform noise). File-hash governance is unaffected.

## PhysioNIRS 0.5.0

- Add source-bound scalp-coupling, PHOEBE-style peak-power, and spectral
  SNR quality metrics with auditable mark/drop channel pruning.
- Add a synchronous governed HbO neurofeedback controller that publishes
  atomic multi-target updates through PhysioStream biofeedback scopes.

## PhysioNIRS 0.4.0

- Add source-bound short-channel identification, centred nearest/all
  short-separation regression, explicit physiology-band extraction, and
  GLM-ready nuisance-design matrices.

## PhysioNIRS 0.3.0

- Add governed motion-artifact detection, TDDR, shift-invariant db2
  wavelet correction, and cubic-spline correction.

## PhysioNIRS 0.2.0

- Add governed natural-log optical density, differential pathlength,
  extinction-coefficient, and modified Beer-Lambert conversions.

## PhysioNIRS 0.1.0

- Add strict SNIRF v1.1 continuous-wave read/write support.
- Add governed measurement-list and source-detector distance helpers.
