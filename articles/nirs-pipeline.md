# A continuous-wave fNIRS pipeline with synthetic data

`PhysioNIRS` carries continuous-wave functional near-infrared
spectroscopy (fNIRS) data in a governed `PhysioExperiment` and preserves
measurement order, probe geometry, the exact time base, and processing
provenance through every step. This vignette walks the standard
preprocessing-to-analysis path

> raw intensity → optical density → haemoglobin concentration (MBLL) →
> motion correction → quality control → short-separation regression →
> activation GLM

entirely offline on a tiny synthetic recording, so no SNIRF file or
device is required.

``` r

library(PhysioNIRS)
#> Loading required package: PhysioExperiment
```

## A small synthetic recording

A governed continuous-wave object needs an intensity assay, a
`time_seconds` column, a SNIRF measurement list, and a probe. The
recording below has two source-detector pairs – one long (3 cm) and one
short (0.8 cm) – each sampled at 760 nm and 850 nm for 16 s at 10 Hz,
with a shared cardiac pulsation, a slow task response, and one motion
spike.

``` r

fs <- 10
n_time <- 160L
t <- seq_len(n_time) / fs
cardiac <- 0.06 * sin(2 * pi * 1.0 * t)
task <- 0.08 * exp(-((t - 8)^2) / (2 * 2^2))
base_od <- c(1.2, 1.0, 1.1, 0.9)
wl_gain <- c(1.0, 0.8, 1.0, 0.8)
od_true <- sapply(1:4, function(j) base_od[j] + cardiac + wl_gain[j] * task)
od_true[round(n_time / 2), ] <- od_true[round(n_time / 2), ] + 0.7   # motion
intensity <- sweep(exp(-od_true), 2L, c(5000, 4000, 5200, 4100), "*")
channel_label <- c("S1_D1_760", "S1_D1_850", "S2_D2_760", "S2_D2_850")
colnames(intensity) <- channel_label

measurement <- S4Vectors::DataFrame(
  measurement_index = 1:4,
  source_index = c(1L, 1L, 2L, 2L), detector_index = c(1L, 1L, 2L, 2L),
  wavelength_index = c(1L, 2L, 1L, 2L), wavelength_nm = c(760, 850, 760, 850),
  wavelength_actual_nm = NA_real_, data_type = 1L, data_type_index = 1L,
  data_type_label = NA_character_, data_unit = NA_character_,
  source_power = NA_real_, detector_gain = NA_real_,
  source_label = c("S1", "S1", "S2", "S2"),
  detector_label = c("D1", "D1", "D2", "D2"), channel_label = channel_label)

probe <- list(
  wavelengths = c(760, 850), sourceLabels = c("S1", "S2"),
  detectorLabels = c("D1", "D2"),
  sourcePos2D = rbind(c(0, 0), c(0.1, 0)),
  detectorPos2D = rbind(c(0, 0.03), c(0.1, 0.008)), LengthUnit = "m")

pe <- PhysioExperiment::PhysioExperiment(
  assays = list(raw = intensity),
  rowData = S4Vectors::DataFrame(time_seconds = t),
  colData = cbind(measurement, label = channel_label),
  metadata = list(snirf = list(measurement_list = measurement, probe = probe,
                               stim = list(), time_sampling = "uniform")),
  samplingRate = fs)
pe <- PhysioExperiment::setEvents(pe, PhysioExperiment::PhysioEvents(
  onset = 6, duration = 4, type = "task", value = "task"))
dim(pe)
#> [1] 160   4
```

The probe geometry yields the source-detector distance of each
measurement:

``` r

sourceDetectorDistances(pe)
#> DataFrame with 4 rows and 9 columns
#>   measurement_index source_index detector_index wavelength_index wavelength_nm
#>           <integer>    <integer>      <integer>        <integer>     <numeric>
#> 1                 1            1              1                1           760
#> 2                 2            1              1                2           850
#> 3                 3            2              2                1           760
#> 4                 4            2              2                2           850
#>   channel_label  distance distance_unit geometry_dimension
#>     <character> <numeric>   <character>        <character>
#> 1     S1_D1_760     0.030             m                 2d
#> 2     S1_D1_850     0.030             m                 2d
#> 3     S2_D2_760     0.008             m                 2d
#> 4     S2_D2_850     0.008             m                 2d
```

## Optical density

[`intensityToOD()`](https://x-biosignal.github.io/PhysioNIRS/reference/intensityToOD.md)
applies the natural-log convention `-log(I / I0)` and records the
reference it used in provenance.

``` r

od <- intensityToOD(pe)
SummarizedExperiment::assayNames(od)
#> [1] "raw" "OD"
```

## Haemoglobin concentration (MBLL)

[`mbll()`](https://x-biosignal.github.io/PhysioNIRS/reference/mbll.md)
collapses each source-detector pair’s two wavelengths into HbO, HbR, and
HbT concentrations using the modified Beer-Lambert law, deriving the
differential pathlength from the probe geometry.

``` r

hb <- mbll(od)
SummarizedExperiment::assayNames(hb)
#> [1] "HbO" "HbR" "HbT"
dim(hb)
#> [1] 160   2
```

## Motion correction

Motion shows up as abrupt shifts in optical density. Detect it, then
correct it with temporal-derivative distribution repair
([`tddr()`](https://x-biosignal.github.io/PhysioNIRS/reference/tddr.md));
spline- and wavelet-based correctors are also available.

``` r

mask <- motionArtifactDetect(od)
corrected <- tddr(od)
SummarizedExperiment::assayNames(corrected)
#> [1] "raw"     "OD"      "OD_tddr"
```

## Signal quality and channel pruning

The scalp-coupling index and a peak-power signal-quality index both
score whether a channel carries a clean cardiac pulsation.
[`pruneChannels()`](https://x-biosignal.github.io/PhysioNIRS/reference/pruneChannels.md)
then marks or drops channels that fail.

``` r

sqi <- signalQualityIndex(od)
sqi
#> NIRS quality <phoebe_peak_power>: 4 channels x 1 windows; pass 4/4
pruned <- pruneChannels(od, sqi, action = "mark")
dim(pruned)
#> [1] 160   4
```

## Short-separation regression

Short channels mostly see scalp haemodynamics. Identify them by
distance, then regress their signal out of the long channels as a
systemic-physiology nuisance.

``` r

short <- identifyShortChannels(od, threshold_m = 0.015)
short[, c("channel_id", "distance_m", "is_short")]
#>    channel_id distance_m is_short
#> 1 S1_D1_wl760      0.030    FALSE
#> 2 S1_D1_wl850      0.030    FALSE
#> 3 S2_D2_wl760      0.008     TRUE
#> 4 S2_D2_wl850      0.008     TRUE
od_ssr <- shortSeparationRegress(od, short = short)
SummarizedExperiment::assayNames(od_ssr)
#> [1] "raw"    "OD"     "OD_ssr"
```

## Activation GLM

Finally, fit a general linear model of the HbO concentration against an
HRF regressor built from the object’s stimulus events, and read off a
per-channel task contrast.

``` r

fit <- nirsActivationGLM(hb, assay_name = "HbO")
fit$coefficients
#>   channel condition       beta        se         t  df        p
#> 1   S1_D1      task 0.02519858 0.1536975 0.1639491 152 0.869989
#> 2   S2_D2      task 0.09449466 0.5763657 0.1639491 152 0.869989
nirsActivationContrast(fit, c(task = 1))
#>   channel   estimate        se         t  df        p
#> 1   S1_D1 0.02519858 0.1536975 0.1639491 152 0.869989
#> 2   S2_D2 0.09449466 0.5763657 0.1639491 152 0.869989
```

## Where to go next

- **Real data.**
  [`readSNIRF()`](https://x-biosignal.github.io/PhysioNIRS/reference/readSNIRF.md)
  and
  [`writeSNIRF()`](https://x-biosignal.github.io/PhysioNIRS/reference/writeSNIRF.md)
  move governed objects to and from SNIRF/HDF5 files, preserving the
  measurement list, probe, stimulus tables, and provenance.
- **Live neurofeedback.**
  [`nirsNeurofeedback()`](https://x-biosignal.github.io/PhysioNIRS/reference/nirsNeurofeedback.md)
  and its lifecycle functions run a closed loop over a regular-rate HbO
  stream; they require a live stream (for example an LSL inlet provided
  by `PhysioStream`) and so are not exercised here.
- **Upstream and downstream packages.** `PhysioExperiment` defines the
  shared data model, and `PhysioStream` provides the streaming
  transports used by the neurofeedback loop. \`\`\`
