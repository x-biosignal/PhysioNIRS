#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PhysioNIRS)
  library(PhysioCore)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script <- normalizePath(file_arg[[1L]], mustWork = TRUE)
validation_dir <- dirname(script)
package_root <- normalizePath(
  file.path(validation_dir, "..", ".."), mustWork = TRUE
)

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE)
}

make_experiment <- function(seed, n_seconds = 180) {
  set.seed(seed)
  fs <- c(8, 10, 12, 20)[(seed %% 4L) + 1L]
  time <- seq.int(0, n_seconds * fs - 1L) / fs
  n <- length(time)
  target_x <- if (seed %% 2L) 0.018 else 0.062
  distances <- c(
    runif(1, 0.005, 0.0075),
    runif(1, 0.008, 0.0095),
    0.01,
    runif(1, 0.025, 0.04)
  )
  centers <- c(0, 0.08, 0.04, target_x)
  source_position <- cbind(centers, 0)
  detector_position <- cbind(centers, distances)
  midpoint <- (source_position + detector_position) / 2
  task <- as.numeric(time >= 60 & time < 120)
  activation <- 0.5 * task
  superficial_1 <- 0.3 + sin(2 * pi * (0.08 + runif(1, 0, 0.025)) * time) +
    0.15 * sin(2 * pi * 0.27 * time + runif(1, 0, 2 * pi))
  superficial_2 <- -0.2 +
    cos(2 * pi * (0.075 + runif(1, 0, 0.025)) * time) +
    0.13 * sin(2 * pi * 0.31 * time + runif(1, 0, 2 * pi))
  nearest_pair <- which.min(sqrt(rowSums(sweep(
    midpoint[1:2, , drop = FALSE], 2L, midpoint[4L, ], "-"
  )^2)))
  nearest_surface <- if (nearest_pair == 1L) superficial_1 else superficial_2
  noise <- outer(
    time,
    seq_len(8L),
    function(t, j) {
      0.006 * sin(2 * pi * (0.43 + 0.011 * j) * t + j + seed / 100)
    }
  )
  values <- cbind(
    superficial_1,
    0.82 * superficial_1 + 0.01 * sin(2 * pi * 0.37 * time),
    superficial_2,
    0.76 * superficial_2 + 0.01 * cos(2 * pi * 0.39 * time),
    0.5 * activation + 1.1 * superficial_1,
    0.4 * activation + 0.9 * 0.82 * superficial_1,
    activation + 1.8 * nearest_surface,
    0.7 * activation + 1.4 * if (nearest_pair == 1L) {
      0.82 * superficial_1 + 0.01 * sin(2 * pi * 0.37 * time)
    } else {
      0.76 * superficial_2 + 0.01 * cos(2 * pi * 0.39 * time)
    }
  ) + noise
  source_index <- rep(1:4, each = 2L)
  detector_index <- source_index
  wavelength_nm <- rep(c(760, 850), 4L)
  labels <- paste0(
    "S", source_index, "_D", detector_index, "_", wavelength_nm
  )
  colnames(values) <- labels
  measurement <- S4Vectors::DataFrame(
    measurement_index = seq_len(8L),
    source_index = source_index,
    detector_index = detector_index,
    wavelength_index = rep(1:2, 4L),
    wavelength_nm = wavelength_nm,
    wavelength_actual_nm = rep(NA_real_, 8L),
    data_type = rep(1L, 8L),
    data_type_index = rep(1L, 8L),
    data_type_label = rep(NA_character_, 8L),
    data_unit = rep(NA_character_, 8L),
    source_power = rep(NA_real_, 8L),
    detector_gain = rep(NA_real_, 8L),
    source_label = paste0("S", source_index),
    detector_label = paste0("D", detector_index),
    channel_label = labels
  )
  metadata <- list(
    snirf = list(
      measurement_list = measurement,
      probe = list(
        wavelengths = c(760, 850),
        sourceLabels = paste0("S", 1:4),
        detectorLabels = paste0("D", 1:4),
        sourcePos2D = source_position,
        detectorPos2D = detector_position,
        LengthUnit = "m"
      ),
      metadata_tags = list(
        SubjectID = "synthetic-validation",
        MeasurementDate = "unknown",
        MeasurementTime = "unknown",
        LengthUnit = "m",
        TimeUnit = "s",
        FrequencyUnit = "Hz"
      ),
      stim = list(),
      time_sampling = "uniform"
    ),
    nirs = list(assays = list(OD = list(
      kind = "optical_density",
      unit = "1",
      log_convention = "natural",
      source_assay = "synthetic"
    )))
  )
  x <- PhysioCore::PhysioExperiment(
    assays = list(OD = values),
    rowData = S4Vectors::DataFrame(time_seconds = time),
    colData = cbind(measurement, label = labels),
    metadata = metadata,
    samplingRate = fs
  )
  x <- PhysioCore::setEvents(
    x,
    PhysioCore::PhysioEvents(
      onset = 60, duration = 60, type = "task", value = "1"
    )
  )
  list(
    x = x,
    fs = fs,
    time = time,
    task = task,
    activation = activation,
    distances = distances,
    midpoint = midpoint,
    nearest_pair = nearest_pair
  )
}

oracle_map <- function(case, threshold = 0.01) {
  identity <- measurementList(case$x)
  midpoint <- case$midpoint[identity$source_index, , drop = FALSE]
  distance <- case$distances[identity$source_index]
  list(
    distance = distance,
    midpoint = midpoint,
    short = distance < threshold
  )
}

nearest_index <- function(midpoint, candidates, target) {
  distance <- sqrt(rowSums(sweep(
    midpoint[candidates, , drop = FALSE],
    2L, midpoint[target, ], "-"
  )^2))
  candidates[[which.min(distance)]]
}

oracle_regress <- function(case, map) {
  source <- SummarizedExperiment::assay(case$x, "OD")
  corrected <- source
  short <- which(map$short)
  long <- which(!map$short)
  wavelength <- measurementList(case$x)$wavelength_nm
  selected <- integer(length(long))
  for (i in seq_along(long)) {
    target <- long[[i]]
    compatible <- short[wavelength[short] == wavelength[[target]]]
    predictor <- nearest_index(map$midpoint, compatible, target)
    selected[[i]] <- predictor
    value <- source[, predictor]
    beta <- stats::cov(source[, target], value) / stats::var(value)
    corrected[, target] <- source[, target] -
      (value - mean(value)) * beta
  }
  list(corrected = corrected, long = long, selected = selected)
}

contrast <- function(value, task) {
  mean(value[task == 1]) - mean(value[task == 0])
}

cnr <- function(value, task) {
  abs(contrast(value, task)) /
    stats::sd(value[task == 0] - mean(value[task == 0]))
}

cases <- vector("list", 100L)
for (i in seq_len(100L)) {
  case <- make_experiment(15000L + i)
  input_before <- serialize(case$x, NULL, version = 3)
  expected_map <- oracle_map(case)
  actual_map <- identifyShortChannels(case$x)
  expected_regression <- oracle_regress(case, expected_map)
  result <- shortSeparationRegress(case$x)
  actual <- SummarizedExperiment::assay(result, "OD_ssr")
  source <- SummarizedExperiment::assay(case$x, "OD")
  target <- 7L
  before_cnr <- cnr(source[, target], case$task)
  after_cnr <- cnr(actual[, target], case$task)
  true_amplitude <- contrast(case$activation, case$task)
  corrected_amplitude <- contrast(actual[, target], case$task)
  design <- shortSeparationDesign(case$x, standardize = FALSE)
  expected_design <- source[, which(expected_map$short), drop = FALSE]
  cases[[i]] <- data.frame(
    case = i,
    seed = 15000L + i,
    fs_hz = case$fs,
    geometry_exact = identical(actual_map$is_short, expected_map$short) &&
      max(abs(actual_map$distance_m - expected_map$distance)) <= 1e-14,
    nearest_exact = identical(
      expected_regression$selected,
      vapply(expected_regression$long, function(target_index) {
        compatible <- which(expected_map$short &
          measurementList(case$x)$wavelength_nm ==
            measurementList(case$x)$wavelength_nm[[target_index]])
        nearest_index(expected_map$midpoint, compatible, target_index)
      }, integer(1))
    ),
    corrected_max_abs = max(abs(actual - expected_regression$corrected)),
    mean_error = max(abs(colMeans(actual) - colMeans(source))),
    cnr_before = before_cnr,
    cnr_after = after_cnr,
    cnr_factor = after_cnr / before_cnr,
    cnr_improved = after_cnr > before_cnr,
    activation_ratio = corrected_amplitude / true_amplitude,
    design_max_abs = max(abs(unname(design) - expected_design)),
    short_unchanged = identical(actual[, which(expected_map$short)],
                                source[, which(expected_map$short)]),
    input_immutable = identical(
      serialize(case$x, NULL, version = 3), input_before
    ),
    finite = all(is.finite(actual)) && all(is.finite(design)),
    stringsAsFactors = FALSE
  )
}
case_table <- do.call(rbind, cases)

filter_case <- make_experiment(15101L, n_seconds = 400)
time <- filter_case$time
pass <- sin(2 * pi * 0.1 * time)
stop_component <- sin(2 * pi * 1.0 * time)
filter_values <- matrix(
  rep(pass + stop_component, 8L), nrow = length(time), ncol = 8L
)
dimnames(filter_values) <- dimnames(
  SummarizedExperiment::assay(filter_case$x, "OD")
)
SummarizedExperiment::assay(
  filter_case$x, "OD", withDimnames = FALSE
) <- filter_values
filtered_object <- physiologyBandpass(
  filter_case$x, band = "mayer", range_hz = c(0.07, 0.13)
)
filtered <- SummarizedExperiment::assay(
  filtered_object, "OD_mayer"
)[, 1L]
middle <- seq.int(filter_case$fs * 50L + 1L, length(time) - 50L * filter_case$fs)
pass_gain <- abs(sum(filtered[middle] * pass[middle]) / sum(pass[middle]^2))
stop_gain <- abs(
  sum(filtered[middle] * stop_component[middle]) /
    sum(stop_component[middle]^2)
)
stop_attenuation_db <- -20 * log10(stop_gain)

caught <- function(expr) {
  tryCatch({
    force(expr)
    FALSE
  }, error = function(e) TRUE)
}

base <- make_experiment(15001L)
x <- base$x
map <- identifyShortChannels(x)
source <- SummarizedExperiment::assay(x, "OD")
out <- shortSeparationRegress(x)
corrected <- SummarizedExperiment::assay(out, "OD_ssr")
oracle <- oracle_regress(base, oracle_map(base))

mm <- x
mm_meta <- S4Vectors::metadata(mm)
mm_meta$snirf$probe$sourcePos2D <- mm_meta$snirf$probe$sourcePos2D * 1000
mm_meta$snirf$probe$detectorPos2D <- mm_meta$snirf$probe$detectorPos2D * 1000
mm_meta$snirf$probe$LengthUnit <- "mm"
S4Vectors::metadata(mm) <- mm_meta
mm_map <- identifyShortChannels(mm)

stale <- x
stale_value <- SummarizedExperiment::assay(stale, "OD")
stale_value[1L, 1L] <- stale_value[1L, 1L] + 1
SummarizedExperiment::assay(stale, "OD", withDimnames = FALSE) <- stale_value

rank_bad <- x
rank_value <- SummarizedExperiment::assay(rank_bad, "OD")
rank_value[, 3L] <- rank_value[, 1L]
rank_value[, 4L] <- rank_value[, 2L]
SummarizedExperiment::assay(rank_bad, "OD", withDimnames = FALSE) <- rank_value

nonfinite <- x
nonfinite_value <- SummarizedExperiment::assay(nonfinite, "OD")
nonfinite_value[1L, 1L] <- Inf
SummarizedExperiment::assay(
  nonfinite, "OD", withDimnames = FALSE
) <- nonfinite_value

population <- sweep(source[, 1:4, drop = FALSE], 2L,
                    colMeans(source[, 1:4, drop = FALSE]), "-")
population <- sweep(
  population, 2L,
  sqrt(colMeans(population^2)), "/"
)
sample_design <- shortSeparationDesign(x)

lag_n <- 300L
lag_short_760 <- sin(seq_len(lag_n) / 7) + 0.2
lag_short_850 <- cos(seq_len(lag_n) / 9) - 0.1
lag_values <- cbind(
  lag_short_760, lag_short_850,
  c(9, lag_short_760[-lag_n]),
  c(-9, lag_short_850[-lag_n]),
  c(7, lag_short_760[-lag_n]),
  c(-7, lag_short_850[-lag_n]),
  c(5, lag_short_760[-lag_n]),
  c(-5, lag_short_850[-lag_n])
)
lag_case <- make_experiment(15003L)
lag_case$x <- lag_case$x[seq_len(lag_n), ]
SummarizedExperiment::rowData(lag_case$x)$time_seconds <-
  seq.int(0, lag_n - 1L) / lag_case$fs
SummarizedExperiment::assay(
  lag_case$x, "OD", withDimnames = FALSE
) <- lag_values
lag_out <- shortSeparationRegress(
  lag_case$x, lag_seconds = 1 / lag_case$fs
)
lag_corrected <- SummarizedExperiment::assay(lag_out, "OD_ssr")

band_a <- SummarizedExperiment::assay(
  physiologyBandpass(x, range_hz = c(0.07, 0.13)), "OD_mayer"
)
band_b <- SummarizedExperiment::assay(
  physiologyBandpass(
    x, range_hz = c(0.08, 0.14), output_assay = "OD_other"
  ),
  "OD_other"
)
band_order_2 <- SummarizedExperiment::assay(
  physiologyBandpass(
    x, range_hz = c(0.07, 0.13), order = 2L,
    output_assay = "OD_order2"
  ),
  "OD_order2"
)

mutations <- data.frame(
  mutation = c(
    "threshold_less_equal", "millimetres_as_metres",
    "stored_distance_trusted", "missing_midpoint_z",
    "separation_used_as_proximity", "tie_not_first",
    "display_name_identity", "cross_wavelength_regression",
    "cross_chromophore_regression", "stale_map_accepted",
    "short_channels_corrected", "intercept_subtracted",
    "target_mean_not_preserved", "regression_direction_reversed",
    "zero_intercept_fit", "positive_lag_reversed",
    "lag_silently_rounded", "lag_edges_zero_filled",
    "rank_deficiency_accepted", "condition_gate_omitted",
    "causal_filter", "wrong_band_edge", "nyquist_equality",
    "wrong_filter_order", "population_sd_standardization",
    "wavelengths_averaged_together", "duplicate_design_names",
    "output_assay_overwritten", "metadata_lost", "nonfinite_propagated",
    "partial_enum_accepted"
  ),
  detected = c(
    !identical(map$is_short, map$distance_m <= 0.01),
    max(abs(mm_map$distance_m - map$distance_m * 1000)) > 1,
    {
      hb <- suppressWarnings(mbll(x, pathlength_factor = 6))
      SummarizedExperiment::colData(hb)$source_detector_distance_m[[1L]] <-
        0.02
      caught(identifyShortChannels(hb))
    },
    {
      bad_map <- map
      bad_map$midpoint_z_m[[1L]] <- NA_real_
      caught(shortSeparationRegress(x, short = bad_map))
    },
    {
      target <- 7L
      compatible <- c(1L, 3L)
      mutant <- compatible[[which.min(abs(
        map$distance_m[compatible] - map$distance_m[[target]]
      ))]]
      mutant != oracle$selected[which(oracle$long == target)]
    },
    {
      tie <- base
      tie_meta <- S4Vectors::metadata(tie$x)
      tie_meta$snirf$probe$sourcePos2D[1L, ] <- c(0, 0.0005)
      tie_meta$snirf$probe$detectorPos2D[1L, ] <- c(0, 0.0085)
      tie_meta$snirf$probe$sourcePos2D[2L, ] <- c(0.08, 0)
      tie_meta$snirf$probe$detectorPos2D[2L, ] <- c(0.08, 0.009)
      tie_meta$snirf$probe$sourcePos2D[4L, ] <- c(0.04, 0)
      tie_meta$snirf$probe$detectorPos2D[4L, ] <- c(0.04, 0.03)
      S4Vectors::metadata(tie$x) <- tie_meta
      tie_result <- SummarizedExperiment::assay(
        shortSeparationRegress(tie$x), "OD_ssr"
      )
      tie_source <- SummarizedExperiment::assay(tie$x, "OD")
      beta_first <- stats::cov(tie_source[, 7L], tie_source[, 1L]) /
        stats::var(tie_source[, 1L])
      first <- tie_source[, 7L] -
        (tie_source[, 1L] - mean(tie_source[, 1L])) * beta_first
      beta_last <- stats::cov(tie_source[, 7L], tie_source[, 3L]) /
        stats::var(tie_source[, 3L])
      last <- tie_source[, 7L] -
        (tie_source[, 3L] - mean(tie_source[, 3L])) * beta_last
      max(abs(tie_result[, 7L] - first)) < 1e-10 &&
        max(abs(tie_result[, 7L] - last)) > 1e-5
    },
    !identical(
      map$channel_id,
      as.character(SummarizedExperiment::colData(x)$label)
    ),
    {
      wrong <- source[, 7L] -
        (source[, 2L] - mean(source[, 2L])) *
        (stats::cov(source[, 7L], source[, 2L]) / stats::var(source[, 2L]))
      max(abs(wrong - corrected[, 7L])) > 1e-5
    },
    {
      hb <- suppressWarnings(mbll(x, pathlength_factor = 6))
      hb_result <- shortSeparationRegress(hb, assay_name = "HbO")
      !identical(
        SummarizedExperiment::assay(hb_result, "HbO_ssr"),
        SummarizedExperiment::assay(hb_result, "HbR")
      )
    },
    caught(shortSeparationRegress(stale, short = map)),
    {
      mutant <- corrected
      mutant[, map$is_short] <- mutant[, map$is_short] + 1
      !identical(mutant[, map$is_short], source[, map$is_short])
    },
    {
      predictor <- source[, oracle$selected[[1L]]]
      beta <- stats::cov(source[, oracle$long[[1L]]], predictor) /
        stats::var(predictor)
      mutant <- source[, oracle$long[[1L]]] -
        (mean(source[, oracle$long[[1L]]]) - beta * mean(predictor)) -
        predictor * beta
      abs(mean(mutant) - mean(source[, oracle$long[[1L]]])) > 1e-4
    },
    max(abs(colMeans(corrected) - colMeans(source))) < 1e-12,
    {
      y <- source[, oracle$long[[1L]]]
      predictor <- source[, oracle$selected[[1L]]]
      reverse_beta <- stats::cov(y, predictor) / stats::var(y)
      mutant <- y - (predictor - mean(predictor)) * reverse_beta
      max(abs(mutant - corrected[, oracle$long[[1L]]])) > 1e-5
    },
    {
      y <- source[, oracle$long[[1L]]]
      predictor <- source[, oracle$selected[[1L]]]
      beta <- sum(predictor * y) / sum(predictor^2)
      mutant <- y - predictor * beta
      max(abs(mutant - corrected[, oracle$long[[1L]]])) > 1e-5
    },
    stats::sd(lag_corrected[-1L, 7L]) < 1e-10 &&
      stats::sd(lag_corrected[-1L, 7L] -
        lag_values[-lag_n, 1L]) > 1e-3,
    caught(shortSeparationRegress(x, lag_seconds = 0.15)),
    lag_corrected[[1L, 7L]] == lag_values[[1L, 7L]] &&
      lag_corrected[[1L, 7L]] != 0,
    caught(shortSeparationRegress(rank_bad, method = "all")),
    caught(shortSeparationRegress(x, condition_limit = 1.01)),
    {
      causal <- as.numeric(signal::filter(
        signal::butter(3L, 2 * c(0.07, 0.13) / base$fs, "pass"),
        source[, 1L]
      ))
      max(abs(causal - band_a[, 1L])) > 1e-4
    },
    max(abs(band_a - band_b)) > 1e-4,
    caught(physiologyBandpass(x, range_hz = c(0.1, base$fs / 2))),
    max(abs(band_a - band_order_2)) > 1e-4,
    max(abs(unname(sample_design) - population)) > 1e-6,
    ncol(shortSeparationDesign(
      x, aggregation = "mean", standardize = FALSE
    )) == 2L,
    anyDuplicated(rep("duplicate", ncol(sample_design))) > 0,
    caught(physiologyBandpass(x, output_assay = "OD")),
    identical(
      S4Vectors::metadata(out)$snirf$probe,
      S4Vectors::metadata(x)$snirf$probe
    ) && !identical(list(), S4Vectors::metadata(x)$snirf$probe),
    caught(shortSeparationRegress(nonfinite)),
    caught(shortSeparationRegress(x, method = "near"))
  ),
  stringsAsFactors = FALSE
)

gates <- data.frame(
  gate = c(
    "cases_100", "geometry_exact_100", "nearest_exact_100",
    "cnr_improved_95", "cnr_median_factor_1_5",
    "activation_median_within_10pct", "corrected_max_abs_1e_10",
    "mne_centered_equation_1e_12", "mean_preserved",
    "passband_gain_5pct", "stopband_attenuation_20db",
    "design_exact", "short_unchanged", "input_immutable", "all_finite",
    "mutations_all_detected", "fixture_manifest_valid"
  ),
  pass = c(
    nrow(case_table) == 100L,
    all(case_table$geometry_exact),
    all(case_table$nearest_exact),
    sum(case_table$cnr_improved) >= 95L,
    stats::median(case_table$cnr_factor) >= 1.5,
    abs(stats::median(case_table$activation_ratio) - 1) <= 0.1,
    max(case_table$corrected_max_abs) <= 1e-10,
    max(case_table$corrected_max_abs) <= 1e-12,
    max(case_table$mean_error) <= 1e-12,
    abs(pass_gain - 1) <= 0.05,
    stop_attenuation_db >= 20,
    max(case_table$design_max_abs) <= 1e-12,
    all(case_table$short_unchanged),
    all(case_table$input_immutable),
    all(case_table$finite),
    nrow(mutations) >= 25L && all(mutations$detected),
    {
      manifest_path <- file.path(
        package_root, "inst", "extdata", "short_separation_reference.sha256"
      )
      fixture_path <- file.path(
        package_root, "inst", "extdata", "short_separation_reference.rds"
      )
      manifest <- readLines(manifest_path, warn = FALSE)
      length(manifest) >= 2L &&
        identical(strsplit(manifest[[1L]], " ", fixed = TRUE)[[1L]][[1L]],
                  sha256_file(fixture_path))
    }
  ),
  stringsAsFactors = FALSE
)

options(digits = 17)
utils::write.csv(
  case_table,
  file.path(validation_dir, "ws10-15-validation.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  mutations,
  file.path(validation_dir, "ws10-15-mutations.csv"),
  row.names = FALSE,
  na = ""
)

report <- c(
  "# WS10-15 independent validation",
  "",
  paste("- Cases:", sum(case_table$geometry_exact), "/ 100 geometry exact;"),
  paste(" ", sum(case_table$nearest_exact), "/ 100 nearest exact"),
  paste("- CNR improved:", sum(case_table$cnr_improved), "/ 100"),
  paste("- Median CNR factor:", format(
    stats::median(case_table$cnr_factor), digits = 17
  )),
  paste("- Median activation ratio:", format(
    stats::median(case_table$activation_ratio), digits = 17
  )),
  paste("- Corrected max absolute error:", format(
    max(case_table$corrected_max_abs), digits = 17
  )),
  paste("- Mean-preservation max error:", format(
    max(case_table$mean_error), digits = 17
  )),
  paste("- Passband gain:", format(pass_gain, digits = 17)),
  paste("- Stopband attenuation dB:", format(
    stop_attenuation_db, digits = 17
  )),
  paste("- Mutations:", sum(mutations$detected), "/", nrow(mutations)),
  paste("- Gates:", sum(gates$pass), "/", nrow(gates)),
  "",
  "CNR improvement is specific to the governed synthetic injection model.",
  "",
  "## Gate table",
  "",
  "| gate | pass |",
  "|---|---|",
  paste0("| ", gates$gate, " | ", gates$pass, " |")
)
writeLines(
  report,
  file.path(validation_dir, "ws10-15-validation.md"),
  useBytes = TRUE
)

if (!all(gates$pass)) {
  print(gates)
  stop("WS10-15 validation gates failed", call. = FALSE)
}
cat(
  "WS10-15 validation:",
  nrow(case_table), "cases;",
  nrow(mutations), "mutations;",
  sum(gates$pass), "/", nrow(gates), "gates PASS\n"
)
