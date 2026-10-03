# Display a governed NIRS neurofeedback runtime

Display a governed NIRS neurofeedback runtime

## Usage

``` r
# S3 method for class 'NIRSNeurofeedback'
print(x, ...)
```

## Arguments

- x:

  A `NIRSNeurofeedback` runtime.

- ...:

  Unused.

## Examples

``` r
# \donttest{
info <- PhysioStream::streamInfo(
  "nirs-demo", type = "NIRS", channel_names = c("S1_D1", "S2_D2"),
  nominal_srate = 10, dtype = "float64", source_id = "nirs-demo",
  clock_domain = "demo-clock", channel_units = c("uM", "uM"),
  metadata = list(nirs = list(
    assay_name = "HbO", assay_kind = "haemoglobin_concentration",
    identity_kind = "pair", channel_id = c("S1_D1", "S2_D2"),
    source_index = 1:2, detector_index = 1:2)))
stream <- PhysioStream::streamOpen(PhysioStream::loopbackSource(info, 8192L))
controller <- nirsNeurofeedback(
  stream, regions = list(motor = "S1_D1", control = "S2_D2"),
  contrast = c(motor = 1, control = -1), baseline_seconds = 1)
print(controller)
#> NIRS neurofeedback <created>: 2 regions -> 1 targets; delivered 0
# }
```
