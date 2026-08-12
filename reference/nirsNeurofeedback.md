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
