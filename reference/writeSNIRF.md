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

## Examples

``` r
# \donttest{
# Write a governed object to a temporary SNIRF/HDF5 file
pe <- PhysioNIRS:::.nirs_demo_object()
pe <- PhysioExperiment::setEvents(pe, PhysioExperiment::PhysioEvents(
  onset = 6, duration = 4, type = "task", value = "1"))
path <- tempfile(fileext = ".snirf")
writeSNIRF(pe, path)
file.exists(path)
#> [1] TRUE
unlink(path)
# }
```
