# Repair fNIRS motion artifacts with TDDR

Implements temporal derivative distribution repair with a third-order
Butterworth low-pass and Tukey bisquare derivative weights.

## Usage

``` r
tddr(
  x,
  assay_name = "OD",
  output_assay = paste0(assay_name, "_tddr"),
  cutoff_hz = 0.5,
  tune = 4.685,
  max_iter = 50L
)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- assay_name:

  Exact source assay name.

- output_assay:

  Exact new assay name.

- cutoff_hz:

  Positive low-pass cutoff below Nyquist.

- tune:

  Positive Tukey bisquare tuning constant.

- max_iter:

  Positive integer maximum robust iterations.

## Value

A clone of `x` with one motion-corrected assay.
