# PhysioNIRS

`PhysioNIRS` provides governed SNIRF input/output, optical-density,
haemoglobin, motion-correction, and short-separation contracts for
near-infrared spectroscopy data represented as `PhysioExperiment` objects.

```r
library(PhysioNIRS)

x <- readSNIRF("recording.snirf")
measurementList(x)
sourceDetectorDistances(x, unit = "m")
writeSNIRF(x, "roundtrip.snirf")
```

For a governed optical-density assay, short-separation processing is explicit:

```r
x_od <- intensityToOD(x, assay_name = "raw")
short <- identifyShortChannels(x_od, threshold_m = 0.01)
corrected <- shortSeparationRegress(
  x_od, assay_name = "OD", short = short
)
nuisance <- shortSeparationDesign(
  x_od,
  assay_name = "OD",
  short = short,
  aggregation = "mean",
  physiology_bands = c("mayer", "respiration")
)
```

Short-channel regression is a signal-processing choice, not evidence that a
channel contains only extracerebral physiology. Named physiology bands are
nuisance ranges rather than diagnoses, and population-appropriate explicit
ranges should be used when the adult defaults are unsuitable.

SNIRF metadata can contain subject identifiers and acquisition dates. Remove or
replace identifying metadata before sharing files.
