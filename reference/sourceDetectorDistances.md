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
