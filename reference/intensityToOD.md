# Convert continuous-wave fNIRS intensity to optical density

Optical density is calculated with the natural-log convention
`-log(I / I0)`. The reference intensity and any replacement of
non-positive cells are explicit and recorded in provenance.

## Usage

``` r
intensityToOD(
  x,
  assay_name = NULL,
  baseline = c("mean", "median"),
  baseline_interval = NULL,
  reference = NULL,
  nonpositive = c("error", "replace"),
  replacement = NULL,
  output_assay = "OD"
)
```

## Arguments

- x:

  A governed continuous-wave `PhysioExperiment`.

- assay_name:

  Exact input assay name. `NULL` is allowed only when one eligible
  intensity assay exists.

- baseline:

  Exactly `"mean"` or `"median"`.

- baseline_interval:

  Optional closed interval in seconds.

- reference:

  Optional positive scalar or vector with one value per measurement.

- nonpositive:

  Exactly `"error"` or `"replace"`.

- replacement:

  Explicit positive scalar used only with `nonpositive = "replace"`.

- output_assay:

  Exact name for the new optical-density assay.

## Value

`x` with a unitless optical-density assay appended.
