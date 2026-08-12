# Construct a governed short-separation nuisance design

Returns short-channel signals, optionally aggregated by compatibility
class, with requested physiology-band copies. No intercept, task, or
drift regressor is added.

## Usage

``` r
shortSeparationDesign(
  x,
  assay_name = "OD",
  short = NULL,
  aggregation = c("all", "mean"),
  physiology_bands = character(),
  band_ranges = NULL,
  standardize = TRUE
)
```

## Arguments

- x:

  A governed OD or pair-collapsed haemoglobin `PhysioExperiment`.

- assay_name:

  Exact source assay.

- short:

  A compatible `nirs_short_channels` map returned by
  [`identifyShortChannels()`](https://x-biosignal.github.io/PhysioNIRS/reference/identifyShortChannels.md),
  or `NULL` for the strict 10-mm default.

- aggregation:

  Exactly `"all"` or `"mean"`.

- physiology_bands:

  Unique exact band names to append.

- band_ranges:

  Optional named list providing one explicit range per requested band.

- standardize:

  Whether to centre and sample-SD standardize every column.

## Value

A time-aligned `nirs_nuisance_design` numeric matrix.
