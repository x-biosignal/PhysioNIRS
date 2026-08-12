# Identify short source-detector channels

Calculates governed source-detector distances and midpoints in metres.
Channels strictly below `threshold_m` are classified as short. A
returned map is fingerprint-bound to the complete source object.

## Usage

``` r
identifyShortChannels(x, threshold_m = 0.01)
```

## Arguments

- x:

  A governed measurement-level or MBLL pair-level `PhysioExperiment`.

- threshold_m:

  Positive short-separation threshold in metres.

## Value

A source-bound `nirs_short_channels` data frame.
