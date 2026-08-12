# Regress governed short-separation nuisance signals

Fits a centred short-channel regression separately for each long
channel. Short channels are copied unchanged. Positive lag means that
the short signal at `t - lag_seconds` predicts the long signal at `t`.

## Usage

``` r
shortSeparationRegress(
  x,
  assay_name = "OD",
  short = NULL,
  method = c("nearest", "all"),
  lag_seconds = 0,
  output_assay = paste0(assay_name, "_ssr"),
  condition_limit = 1e+10
)
```

## Arguments

- x:

  A governed OD or pair-collapsed haemoglobin `PhysioExperiment`.

- assay_name:

  Exact source assay name.

- short:

  A compatible `nirs_short_channels` map returned by
  [`identifyShortChannels()`](https://x-biosignal.github.io/PhysioNIRS/reference/identifyShortChannels.md),
  or `NULL` to use the strict 10-mm default.

- method:

  Exactly `"nearest"` or `"all"`.

- lag_seconds:

  A finite lag that corresponds to an exact sample count.

- output_assay:

  Exact name for the added corrected assay.

- condition_limit:

  Finite regression condition-number limit greater than one.

## Value

A clone of `x` with one short-separation-corrected assay.

## Details

For a target `y` and selected short signals `X`, the fitted intercept is
retained and only `(X - colMeans(X)) beta` is subtracted. Regression
uses a deterministic SVD solve and rejects rank-deficient or over-limit
designs. Short columns and lag-edge rows remain unchanged.

## References

Saager and Berger (2005), DOI: 10.1364/JOSAA.22.001874.
