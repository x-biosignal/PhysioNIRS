# Contrast of an fNIRS activation GLM

Applies a linear contrast over the condition betas of a
[`nirsActivationGLM()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationGLM.md)
fit, per channel: the estimate `c' beta`, its standard error
`sqrt(sigma^2 * c' (X'X)^-1 c)`, t-statistic, degrees of freedom and
p-value.

## Usage

``` r
nirsActivationContrast(fit, contrast)
```

## Arguments

- fit:

  A `nirs_activation_glm`.

- contrast:

  A numeric contrast over conditions: named by condition, or an unnamed
  vector with one weight per condition (e.g. `c(task = 1, rest = -1)`).

## Value

A per-channel `data.frame`: `channel`, `estimate`, `se`, `t`, `df`, `p`.

## See also

[`nirsActivationGLM()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsActivationGLM.md)

## Examples

``` r
hb <- mbll(intensityToOD(PhysioNIRS:::.nirs_demo_object()))
fit <- nirsActivationGLM(hb, assay_name = "HbO")
nirsActivationContrast(fit, c(task = 1))
#>   channel   estimate        se         t  df         p
#> 1   S1_D1 0.03532496 0.1592257 0.2218547 152 0.8247248
#> 2   S2_D2 0.13246861 0.5970962 0.2218547 152 0.8247248
```
