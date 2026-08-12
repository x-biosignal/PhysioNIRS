# Convert optical density to haemoglobin concentration change

Solves the modified Beer-Lambert system independently for each exact
source-detector pair. The returned `PhysioExperiment` has one column per
pair and assays `HbO`, `HbR`, and `HbT` in micromolar.

## Usage

``` r
mbll(
  x,
  assay_name = "OD",
  pathlength_factor = 6,
  age_years = NULL,
  dpf_model = "scholkmann2013",
  distance_m = NULL,
  output_unit = "uM"
)
```

## Arguments

- x:

  A governed `PhysioExperiment` carrying natural-log optical density.

- assay_name:

  Exact optical-density assay name.

- pathlength_factor:

  Positive scalar/vector, or `NULL` to derive DPF from `age_years`.

- age_years:

  Age used only when `pathlength_factor = NULL`.

- dpf_model:

  Exact DPF model name passed to
  [`dpf()`](https://x-biosignal.github.io/PhysioNIRS/reference/dpf.md).

- distance_m:

  Optional positive scalar, per-measurement vector, or per-pair vector
  overriding probe geometry.

- output_unit:

  Exactly `"uM"`.

## Value

A pair-collapsed `PhysioExperiment` with concentration assays.
