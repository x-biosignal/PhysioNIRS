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
