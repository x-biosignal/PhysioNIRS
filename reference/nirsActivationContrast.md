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
