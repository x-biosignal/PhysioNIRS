# Extract an explicit fNIRS physiology frequency band

Applies a zero-phase Butterworth band-pass independently to every source
assay column. Band names identify common nuisance ranges; they are not
physiological diagnoses.

## Usage

``` r
physiologyBandpass(
  x,
  assay_name = "OD",
  band = c("mayer", "respiration", "cardiac"),
  range_hz = NULL,
  output_assay = paste0(assay_name, "_", band[[1L]]),
  order = 3L
)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- assay_name:

  Exact source assay name.

- band:

  Exactly `"mayer"`, `"respiration"`, or `"cardiac"`.

- range_hz:

  Optional explicit positive passband below Nyquist.

- output_assay:

  Exact name for the added filtered assay.

- order:

  Positive integer Butterworth order.

## Value

A clone of `x` with one zero-phase filtered assay.

## Details

Default passbands in Hz are Mayer `0.07--0.13`, respiration
`0.15--0.40`, and cardiac `0.70--2.00`. Filtering uses an order-`order`
Butterworth band-pass through
[`signal::filtfilt()`](https://rdrr.io/pkg/signal/man/filtfilt.html).
Explicit ranges must remain strictly below the governed Nyquist
frequency.
