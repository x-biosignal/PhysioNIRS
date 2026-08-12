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
[`PhysioCore::PhysioExperiment()`](https://x-biosignal.r-universe.dev/PhysioCore/reference/PhysioExperiment.html)
with a time-by-measurement `raw` assay.
