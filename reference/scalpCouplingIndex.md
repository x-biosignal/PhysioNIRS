# Calculate the governed fNIRS scalp-coupling index

Optical-density measurements are grouped by exact source-detector
identity, zero-phase filtered in the cardiac band, and scored by the
minimum pairwise wavelength correlation. Scores are copied to every
wavelength in the group. The default threshold is a project default, not
a clinical cutoff.

## Usage

``` r
scalpCouplingIndex(
  x,
  assay_name = "OD",
  l_freq = 0.7,
  h_freq = 1.5,
  window_seconds = NULL,
  step_seconds = NULL,
  threshold = 0.8,
  order = 3L
)
```

## Arguments

- x:

  A governed measurement-level optical-density `PhysioExperiment`.

- assay_name:

  Exact optical-density assay name.

- l_freq, h_freq:

  Cardiac-band limits in Hz.

- window_seconds:

  Optional positive complete-window duration.

- step_seconds:

  Optional positive window step. `NULL` uses the window.

- threshold:

  Finite pass threshold in `[-1, 1]`.

- order:

  Positive Butterworth filter order.

## Value

A source-bound `nirs_quality` result.

## References

Pollonini et al. (2014), DOI: 10.1117/1.JBO.19.8.086007.
