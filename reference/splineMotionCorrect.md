# Correct detected fNIRS motion intervals with cubic smoothing splines

Uses the MATLAB `csaps` smoothing-parameter definition and stitches
corrected artifact and clean segments with short local means.

## Usage

``` r
splineMotionCorrect(
  x,
  mask = NULL,
  assay_name = "OD",
  output_assay = paste0(assay_name, "_spline"),
  p = 0.99,
  t_motion = 0.5,
  t_mask = 1,
  sd_threshold = 20,
  amplitude_threshold = 0.5,
  group_wavelengths = TRUE
)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- mask:

  A compatible `nirs_motion_mask`, or `NULL` to detect once.

- assay_name:

  Exact source assay name.

- output_assay:

  Exact new assay name.

- p:

  Cubic smoothing parameter in `[0, 1]`.

- t_motion:

  Positive detection lag window in seconds, used only when
  `mask = NULL`.

- t_mask:

  Non-negative detection-mask expansion in seconds, used only when
  `mask = NULL`.

- sd_threshold:

  Positive derivative-SD detection multiplier, used only when
  `mask = NULL`.

- amplitude_threshold:

  Positive absolute-change detection threshold, used only when
  `mask = NULL`.

- group_wavelengths:

  Whether detection unions wavelengths from each source-detector pair,
  used only when `mask = NULL`.

## Value

A clone of `x` with one motion-corrected assay.
