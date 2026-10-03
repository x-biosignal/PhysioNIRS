# Control and inspect live NIRS neurofeedback

`nirsNeurofeedbackStep()` processes at most `max_chunks` synchronously,
flushing committed target values after each chunk so the external update
timestamp remains the latest scope signal time.

## Usage

``` r
nirsNeurofeedbackStart(controller)

nirsNeurofeedbackStep(controller, max_chunks = 16L)

nirsNeurofeedbackState(controller)

nirsNeurofeedbackStop(controller)

nirsNeurofeedbackScope(controller)
```

## Arguments

- controller:

  A `NIRSNeurofeedback` runtime.

- max_chunks:

  Positive exact per-call work bound.

## Value

Lifecycle functions return `controller` invisibly. Step and state
functions return portable plain lists. `nirsNeurofeedbackScope()`
returns the governed `BiofeedbackScope` for a viewer.

## Examples

``` r
# \donttest{
# A controller is driven by an open regular-rate HbO stream; here a loopback
# stream stands in for a live device inlet.
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
nirsNeurofeedbackStart(controller)
time <- (0:19) / 10
PhysioStream::loopbackFeed(stream, cbind(
  S1_D1 = 1 + 0.2 * sin(2 * pi * 0.5 * time), S2_D2 = 1), time)
while (nirsNeurofeedbackStep(controller, 4L)$updated) {}
nirsNeurofeedbackScope(controller)
#> <BiofeedbackScope: state=running, traces=3, frames=20, time=1.9>
nirsNeurofeedbackState(controller)$lifecycle
#> [1] "running"
nirsNeurofeedbackStop(controller)
# }
```
