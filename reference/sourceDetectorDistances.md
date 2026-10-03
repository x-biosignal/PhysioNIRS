# Calculate source-detector distances for each SNIRF measurement

Calculate source-detector distances for each SNIRF measurement

## Usage

``` r
sourceDetectorDistances(
  x,
  dimension = c("auto", "3d", "2d"),
  unit = c("m", "native")
)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- dimension:

  Geometry selection: exactly `"auto"`, `"3d"`, or `"2d"`.

- unit:

  Output unit: exactly `"m"` or `"native"`.

## Value

A `DataFrame` with one row per measurement and distance identity.

## Examples

``` r
pe <- PhysioNIRS:::.nirs_demo_object()
sourceDetectorDistances(pe)
#> DataFrame with 4 rows and 9 columns
#>   measurement_index source_index detector_index wavelength_index wavelength_nm
#>           <integer>    <integer>      <integer>        <integer>     <numeric>
#> 1                 1            1              1                1           760
#> 2                 2            1              1                2           850
#> 3                 3            2              2                1           760
#> 4                 4            2              2                2           850
#>   channel_label  distance distance_unit geometry_dimension
#>     <character> <numeric>   <character>        <character>
#> 1     S1_D1_760     0.030             m                 2d
#> 2     S1_D1_850     0.030             m                 2d
#> 3     S2_D2_760     0.008             m                 2d
#> 4     S2_D2_850     0.008             m                 2d
```
