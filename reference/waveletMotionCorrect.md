# Correct fNIRS motion with a shift-invariant db2 wavelet transform

Outlying detail coefficients are removed using a per-block IQR rule and
the signal is reconstructed at its original length.

## Usage

``` r
waveletMotionCorrect(
  x,
  assay_name = "OD",
  output_assay = paste0(assay_name, "_wavelet"),
  iqr = 1.5,
  wavelet = "db2",
  min_level = 4L
)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- assay_name:

  Exact source assay name.

- output_assay:

  Exact new assay name.

- iqr:

  Positive IQR multiplier.

- wavelet:

  Exact wavelet name; currently only `"db2"`.

- min_level:

  Positive integer lowest analysis scale.

## Value

A clone of `x` with one motion-corrected assay.
