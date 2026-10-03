# Construct a governed live fNIRS neurofeedback controller

The controller builds a public PhysioStream pipeline and
`BiofeedbackScope`. It computes frozen-baseline regional HbO contrasts
in a committed pipeline operation and publishes each eligible
multi-target value through one atomic `biofeedbackUpdate()` call.
Construction is side-effect-free.

## Usage

``` r
nirsNeurofeedback(
  stream,
  regions,
  contrast,
  assay_name = "HbO",
  baseline_seconds = 20,
  update_seconds = 0.1,
  smoothing_seconds = 2
)
```

## Arguments

- stream:

  A created or open regular-rate numeric `StreamSource`.

- regions:

  A named list of non-overlapping exact HbO channel IDs.

- contrast:

  A region-named numeric vector or region-by-target matrix.

- assay_name:

  Exact live assay name, currently `"HbO"`.

- baseline_seconds:

  Positive frozen-baseline duration.

- update_seconds:

  Positive scheduled update interval.

- smoothing_seconds:

  Positive trailing time-weighted mean duration.

## Value

A side-effect-free `NIRSNeurofeedback` runtime.

## Details

The stream must declare `metadata$nirs` fields `assay_name`,
`assay_kind = "haemoglobin_concentration"`, `identity_kind = "pair"`,
`channel_id`, `source_index`, and `detector_index`; ordered channel
names and units must be the exact HbO identities and `"uM"`.

## References

Kober et al. (2017), DOI: 10.3389/fnhum.2017.00081.

## Examples

``` r
# \donttest{
# Offline demonstration on a deterministic loopback stream. In practice
# `stream` is a live regular-rate HbO inlet (e.g. an LSL inlet from a device).
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
  contrast = c(motor = 1, control = -1),
  baseline_seconds = 1, update_seconds = 0.1, smoothing_seconds = 0.5)
nirsNeurofeedbackStart(controller)
time <- (0:19) / 10
PhysioStream::loopbackFeed(stream, cbind(
  S1_D1 = 1 + 0.2 * sin(2 * pi * 0.5 * time),
  S2_D2 = 1 + 0.1 * cos(time)), time)
while (nirsNeurofeedbackStep(controller, 4L)$updated) {}
nirsNeurofeedbackState(controller)$lifecycle
#> [1] "running"
nirsNeurofeedbackStop(controller)
# }
```
