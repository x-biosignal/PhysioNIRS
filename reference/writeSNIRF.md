# Write a Shared Near Infrared Spectroscopy Format file

Write a Shared Near Infrared Spectroscopy Format file

## Usage

``` r
writeSNIRF(x, path, assay_name = NULL, overwrite = FALSE, compact_time = FALSE)
```

## Arguments

- x:

  A governed `PhysioExperiment`.

- path:

  Destination ending exactly in `.snirf`.

- assay_name:

  Assay to write. `NULL` uses the default assay.

- overwrite:

  Whether to replace an existing destination.

- compact_time:

  Whether to use `[start, spacing]` time encoding.

## Value

`path`, invisibly.
