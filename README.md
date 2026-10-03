# PhysioNIRS

`PhysioNIRS` provides governed SNIRF input/output, optical-density,
haemoglobin, motion-correction, and short-separation contracts for
near-infrared spectroscopy data represented as `PhysioExperiment` objects.

```r
library(PhysioNIRS)

# A continuous-wave SNIRF file bundled with the package (no download needed).
snirf <- system.file(
  "extdata", "snirf_official_simple_probe.snirf", package = "PhysioNIRS"
)
x <- readSNIRF(snirf)
measurementList(x)
sourceDetectorDistances(x, unit = "m")

x_od <- intensityToOD(x, assay_name = "raw")   # raw intensity -> optical density
x_hb <- mbll(x_od, assay_name = "OD")          # optical density -> haemoglobin

writeSNIRF(x, file.path(tempdir(), "roundtrip.snirf"))
```

For a governed optical-density assay, short-separation processing is explicit:

```r
short <- identifyShortChannels(x_od, threshold_m = 0.01)

# When the probe includes short channels, regress them out and build a
# physiology nuisance design. (This symmetric example probe has only
# standard-separation channels, so the step below is skipped here.)
if (any(short$is_short) && !all(short$is_short)) {
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
}
```

Short-channel regression is a signal-processing choice, not evidence that a
channel contains only extracerebral physiology. Named physiology bands are
nuisance ranges rather than diagnoses, and population-appropriate explicit
ranges should be used when the adult defaults are unsuitable.

SNIRF metadata can contain subject identifiers and acquisition dates. Remove or
replace identifying metadata before sharing files.

## Installation

The ecosystem builds on Bioconductor, so its repositories have to be on the
list as well -- without them the install stops at `SummarizedExperiment`.

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages(
  "PhysioNIRS",
  repos = c("https://x-biosignal.r-universe.dev", BiocManager::repositories())
)
```

From GitHub instead:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("x-biosignal/PhysioNIRS")
```

