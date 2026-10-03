# Read a Shared Near Infrared Spectroscopy Format file

Read a Shared Near Infrared Spectroscopy Format file

## Usage

``` r
readSNIRF(path, nirs_index = 1L, data_index = 1L)
```

## Arguments

- path:

  Path to one readable `.snirf` file.

- nirs_index:

  Positive index of the NIRS root.

- data_index:

  Positive index of the data block.

## Value

A
[`PhysioExperiment::PhysioExperiment()`](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/PhysioExperiment.html)
with a time-by-measurement `raw` assay.

## Examples

``` r
# \donttest{
# Round-trip through a temporary SNIRF/HDF5 file (needs a .snirf on disk)
pe <- PhysioNIRS:::.nirs_demo_object()
pe <- PhysioExperiment::setEvents(pe, PhysioExperiment::PhysioEvents(
  onset = 6, duration = 4, type = "task", value = "1"))
path <- tempfile(fileext = ".snirf")
writeSNIRF(pe, path)
back <- readSNIRF(path)
dim(back)
#> [1] 160   4
unlink(path)
# }
```
