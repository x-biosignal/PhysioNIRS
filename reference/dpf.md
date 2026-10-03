# Estimate a differential pathlength factor

Implements published population equations. These values are population
estimates, not subject-specific pathlength measurements.

## Usage

``` r
dpf(
  wavelength_nm,
  age_years,
  model = c("scholkmann2013", "duncan1996"),
  extrapolate = FALSE
)
```

## Arguments

- wavelength_nm:

  Positive wavelength values in nanometres.

- age_years:

  Non-negative age values in years.

- model:

  Exactly `"scholkmann2013"` or `"duncan1996"`.

- extrapolate:

  Whether to allow age/wavelength extrapolation where the selected model
  defines an unambiguous equation.

## Value

A positive numeric vector with model and extrapolation attributes.

## Examples

``` r
# Age-dependent differential pathlength factors for a 30-year-old adult
dpf(c(760, 850), age_years = 30)
#> [1] 6.296601 5.234328
#> attr(,"model")
#> [1] "scholkmann2013"
#> attr(,"extrapolated")
#> [1] FALSE FALSE
```
