# Mark or drop channels using governed NIRS quality results

Window decisions are collapsed conservatively per metric.
Measurement-level optical-density drops propagate to the complete
source-detector pair. Inputs and quality results are fingerprint-bound
and are never modified.

## Usage

``` r
pruneChannels(x, quality, action = c("mark", "drop"), require_all = TRUE)
```

## Arguments

- x:

  The exact governed source `PhysioExperiment`.

- quality:

  One `nirs_quality` result or a non-empty named list of distinct
  results.

- action:

  Exactly `"mark"` or `"drop"`.

- require_all:

  Whether all metrics, rather than any metric, must pass.

## Value

A cloned marked or consistently subset `PhysioExperiment`.
