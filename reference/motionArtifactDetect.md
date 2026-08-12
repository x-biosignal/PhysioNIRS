# Detect fNIRS motion artifacts by channel

Uses strict amplitude and derivative-SD thresholds over a finite lag
window. Detections are expanded in time and can be unioned across
wavelengths from the same source-detector pair.

## Usage

``` r
motionArtifactDetect(
  x,
  assay_name = "OD",
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

- assay_name:

  Exact source assay name.

- t_motion:

  Positive lag window in seconds.

- t_mask:

  Non-negative expansion on each side in seconds.

- sd_threshold:

  Positive multiple of channel derivative SD.

- amplitude_threshold:

  Positive absolute-change threshold.

- group_wavelengths:

  Whether to union measurements sharing an exact source-detector pair.

## Value

A typed `nirs_motion_mask` aligned to assay rows and columns.
