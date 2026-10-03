# Display a governed NIRS quality result

Display a governed NIRS quality result

## Usage

``` r
# S3 method for class 'nirs_quality'
print(x, ...)
```

## Arguments

- x:

  A `nirs_quality` result.

- ...:

  Unused.

## Examples

``` r
od <- intensityToOD(PhysioNIRS:::.nirs_demo_object())
print(signalQualityIndex(od))
#> NIRS quality <phoebe_peak_power>: 4 channels x 1 windows; pass 4/4
```
