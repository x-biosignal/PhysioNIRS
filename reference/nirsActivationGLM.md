# First-level fNIRS cortical-activation GLM

Fits a general linear model of every haemoglobin channel against
HRF-convolved task regressors: for each condition, a boxcar over its
stimulus blocks convolved with a canonical double-gamma HRF, plus a
discrete-cosine drift basis and optional short-separation nuisance
regressors. Returns the per-channel activation amplitude (beta),
standard error and t-statistic for each condition. Run it separately on
the `HbO` and `HbR` assays.

## Usage

``` r
nirsActivationGLM(
  x,
  events = NULL,
  assay_name = "HbO",
  conditions = NULL,
  basis = c("canonical", "derivative"),
  drift = 6,
  nuisance = NULL,
  hrf = list(peak = 6, undershoot = 16, ratio = 6),
  prewhiten = c("none", "ar1")
)
```

## Arguments

- x:

  A `PhysioExperiment` of haemoglobin concentration (e.g. the `HbO`
  assay from
  [`mbll()`](https://x-biosignal.github.io/PhysioNIRS/reference/mbll.md)).

- events:

  A `data.frame` with numeric `onset` and `duration` (seconds) and a
  `condition` column; if `NULL`, the object's events
  ([PhysioCore::getEvents](https://x-biosignal.github.io/PhysioCore//reference/getEvents.html)
  `onset`/`duration`/`value`) are used, with `value` as the condition.

- assay_name:

  Assay to model (default `"HbO"`).

- conditions:

  Optional character vector selecting/ordering conditions.

- basis:

  `"canonical"` (HRF only) or `"derivative"` (HRF + temporal derivative
  per condition).

- drift:

  Number of discrete-cosine drift terms (default 6).

- nuisance:

  Optional numeric matrix of extra nuisance regressors (e.g.
  [`shortSeparationDesign()`](https://x-biosignal.github.io/PhysioNIRS/reference/shortSeparationDesign.md)
  output).

- hrf:

  List of HRF parameters (`peak`, `undershoot`, `ratio`).

- prewhiten:

  `"none"` or `"ar1"` (per-channel AR(1) prewhitening; strongly
  recommended for fNIRS to avoid inflated t-statistics).

## Value

An S3 `nirs_activation_glm` with `coefficients` (a `channel` x
`condition` data.frame of `beta`, `se`, `t`, `df`, `p`), the `beta`/`se`
matrices, per-channel `sigma2`/`rho`, and the `design` (for
[`nirsActivationContrast()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationContrast.md)).

## References

Friston KJ et al. (1994) Statistical parametric maps. *Hum Brain Mapp*
2:189-210. Huppert TJ et al. (2009) HomER. *Appl Opt* 48:D280.

## See also

[`nirsActivationContrast()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationContrast.md),
[`shortSeparationDesign()`](https://x-biosignal.github.io/PhysioNIRS/reference/shortSeparationDesign.md),
[`mbll()`](https://x-biosignal.github.io/PhysioNIRS/reference/mbll.md)
