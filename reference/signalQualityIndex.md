# Calculate governed fNIRS signal quality

`method = "peak_power"` implements the pinned MNE-NIRS PHOEBE-style
cross-correlation peak-power metric. It is not the original PHOEBE
equation. `method = "snr"` returns a channel-specific cardiac-to-noise
spectral power ratio in decibels.

## Usage

``` r
signalQualityIndex(
  x,
  assay_name = "OD",
  method = c("peak_power", "snr"),
  window_seconds = 10,
  step_seconds = NULL,
  cardiac_range_hz = c(0.7, 1.5),
  noise_range_hz = NULL,
  threshold = NULL,
  order = 3L
)
```

## Arguments

- x:

  A governed measurement-level optical-density `PhysioExperiment`.

- assay_name:

  Exact optical-density assay name.

- method:

  Exactly `"peak_power"` or `"snr"`.

- window_seconds:

  Positive complete-window duration.

- step_seconds:

  Optional positive window step. `NULL` uses the window.

- cardiac_range_hz:

  Finite increasing cardiac band below Nyquist.

- noise_range_hz:

  For SNR, `NULL` or a disjoint two-column band matrix.

- threshold:

  Optional finite threshold. Defaults to `0.1` for peak power and `0` dB
  for SNR.

- order:

  Positive Butterworth filter order.

## Value

A source-bound `nirs_quality` result.

## References

Pollonini et al. (2016), DOI: 10.1364/BOE.7.005104.
