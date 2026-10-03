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

## Examples

``` r
od <- intensityToOD(PhysioNIRS:::.nirs_demo_object())
q <- signalQualityIndex(od)
pruneChannels(od, q, action = "mark")
#> class: PhysioExperiment
#> dim: 160 x 4 
#> assays(2): raw, OD
#> samplingRate: 10 Hz
#> channels(4): S1_D1_760, S1_D1_850, S2_D2_760, S2_D2_850
#> rowData names(1): time_seconds
#> colData names(20): measurement_index, source_index, detector_index, wavelength_index, wavelength_nm ...
#> events: 1 
#> provenance: 2 steps
```
