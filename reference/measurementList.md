# Return the governed SNIRF measurement list

Return the governed SNIRF measurement list

## Usage

``` r
measurementList(x)
```

## Arguments

- x:

  A `PhysioExperiment` imported from SNIRF or carrying the same governed
  metadata contract.

## Value

A defensive copy of the normalized measurement `DataFrame`.

## Examples

``` r
pe <- PhysioNIRS:::.nirs_demo_object()
measurementList(pe)[, c("channel_label", "wavelength_nm")]
#> DataFrame with 4 rows and 2 columns
#>   channel_label wavelength_nm
#>     <character>     <numeric>
#> 1     S1_D1_760           760
#> 2     S1_D1_850           850
#> 3     S2_D2_760           760
#> 4     S2_D2_850           850
```
