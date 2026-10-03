# Get governed haemoglobin extinction coefficients

Coefficients are linearly interpolated from the package-owned table.
Stored coefficients use the base-10 molar convention in `cm^-1 M^-1`;
natural-log coefficients in `m^-1 M^-1` include the explicit
`100 * log(10)` conversion.

## Usage

``` r
extinctionCoefficients(
  wavelength_nm,
  unit = c("m-1 M-1", "cm-1 M-1"),
  extrapolate = FALSE
)
```

## Arguments

- wavelength_nm:

  Positive wavelengths in nanometres.

- unit:

  Exactly `"m-1 M-1"` or `"cm-1 M-1"`.

- extrapolate:

  Whether to permit terminal-slope linear extrapolation.

## Value

A numeric matrix with exact columns `HbO` and `HbR`.

## Examples

``` r
# Molar extinction (m^-1 M^-1) at two typical fNIRS wavelengths
extinctionCoefficients(c(760, 850))
#>           HbO      HbR
#> [1,] 134931.5 356559.9
#> [2,] 243613.5 159182.3
#> attr(,"wavelength_nm")
#> [1] 760 850
#> attr(,"unit")
#> [1] "m-1 M-1"
#> attr(,"extrapolated")
#> [1] FALSE FALSE
# The same coefficients expressed per centimetre
extinctionCoefficients(c(690, 830), unit = "cm-1 M-1")
#>      HbO     HbR
#> [1,] 276 2051.96
#> [2,] 974  693.04
#> attr(,"wavelength_nm")
#> [1] 690 830
#> attr(,"unit")
#> [1] "cm-1 M-1"
#> attr(,"extrapolated")
#> [1] FALSE FALSE
```
