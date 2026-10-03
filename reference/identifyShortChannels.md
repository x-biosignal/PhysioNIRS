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

## Examples

``` r
od <- intensityToOD(PhysioNIRS:::.nirs_demo_object())
identifyShortChannels(od, threshold_m = 0.015)
#>   channel_index  channel_id source_index detector_index distance_m midpoint_x_m
#> 1             1 S1_D1_wl760            1              1      0.030          0.0
#> 2             2 S1_D1_wl850            1              1      0.030          0.0
#> 3             3 S2_D2_wl760            2              2      0.008          0.1
#> 4             4 S2_D2_wl850            2              2      0.008          0.1
#>   midpoint_y_m midpoint_z_m is_short identity_kind
#> 1        0.015            0    FALSE   measurement
#> 2        0.015            0    FALSE   measurement
#> 3        0.004            0     TRUE   measurement
#> 4        0.004            0     TRUE   measurement
```
