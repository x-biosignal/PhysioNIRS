.snirf_required_tags <- c(
  "SubjectID", "MeasurementDate", "MeasurementTime",
  "LengthUnit", "TimeUnit", "FrequencyUnit"
)

.snirf_measurement_fields <- c(
  "sourceIndex", "detectorIndex", "wavelengthIndex", "dataType",
  "dataTypeIndex", "dataTypeLabel", "dataUnit", "wavelengthActual",
  "sourcePower", "detectorGain"
)

.snirf_probe_fields <- c(
  "wavelengths", "wavelengthsEmission",
  "sourcePos2D", "detectorPos2D", "sourcePos3D", "detectorPos3D",
  "sourceLabels", "detectorLabels",
  "landmarkPos2D", "landmarkPos3D", "landmarkLabels",
  "coordinateSystem", "coordinateSystemDescription"
)

.snirf_enum <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !x %in% choices) {
    stop("`", arg, "` must be exactly one of: ",
         paste(choices, collapse = ", "), call. = FALSE)
  }
  x
}

.snirf_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be one non-missing logical value", call. = FALSE)
  }
  x
}

.snirf_index <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 1 || x > .Machine$integer.max ||
      x != floor(x)) {
    stop("`", arg, "` must be one finite positive whole number",
         call. = FALSE)
  }
  as.integer(x)
}

.snirf_tree <- function(path) {
  tree <- rhdf5::h5ls(path, recursive = TRUE, all = TRUE)
  tree$.path <- ifelse(
    tree$group == "/", paste0("/", tree$name),
    paste0(tree$group, "/", tree$name)
  )
  tree
}

.snirf_tree_row <- function(tree, path, type = NULL) {
  row <- tree[tree$.path == path, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Required SNIRF object is missing or ambiguous: ", path,
         call. = FALSE)
  }
  if (!is.null(type) && !identical(as.character(row$otype), type)) {
    stop("SNIRF object has the wrong HDF5 type: ", path, call. = FALSE)
  }
  row
}

.snirf_has <- function(tree, path, type = NULL) {
  row <- tree[tree$.path == path, , drop = FALSE]
  nrow(row) == 1L && (is.null(type) || identical(as.character(row$otype), type))
}

.snirf_rank <- function(tree, path) {
  as.integer(.snirf_tree_row(tree, path, "H5I_DATASET")$rank)
}

.snirf_read_scalar <- function(path, dataset, tree, kind = c("any", "string",
                                                              "numeric")) {
  kind <- .snirf_enum(kind[[1L]], c("any", "string", "numeric"), "kind")
  if (.snirf_rank(tree, dataset) != 0L) {
    stop("SNIRF field must use a scalar HDF5 dataspace: ", dataset,
         call. = FALSE)
  }
  value <- rhdf5::h5read(path, dataset)
  if (length(value) != 1L || is.list(value)) {
    stop("SNIRF scalar is unreadable: ", dataset, call. = FALSE)
  }
  if (kind == "string") {
    value <- as.character(value)
    if (is.na(value) || !nzchar(value)) {
      stop("SNIRF string scalar must be non-empty: ", dataset, call. = FALSE)
    }
  } else if (kind == "numeric") {
    if (!is.numeric(value) || is.na(value) || !is.finite(value)) {
      stop("SNIRF numeric scalar must be finite: ", dataset, call. = FALSE)
    }
    value <- as.numeric(value)
  }
  value
}

.snirf_read_vector <- function(path, dataset, tree, kind = c("numeric",
                                                              "string"),
                               allow_empty = FALSE, allow_missing = FALSE) {
  kind <- .snirf_enum(kind[[1L]], c("numeric", "string"), "kind")
  if (.snirf_rank(tree, dataset) != 1L) {
    stop("SNIRF field must be a rank-1 array: ", dataset, call. = FALSE)
  }
  value <- as.vector(rhdf5::h5read(path, dataset))
  if (!allow_empty && !length(value)) {
    stop("SNIRF array must not be empty: ", dataset, call. = FALSE)
  }
  if (kind == "numeric") {
    invalid <- if (allow_missing) {
      !is.na(value) & !is.finite(value)
    } else {
      is.na(value) | !is.finite(value)
    }
    if (!is.numeric(value) || any(invalid)) {
      stop("SNIRF numeric array must be finite: ", dataset, call. = FALSE)
    }
    value <- as.numeric(value)
  } else {
    value <- as.character(value)
    if (allow_missing) {
      value[is.na(value) | !nzchar(value)] <- NA_character_
    } else if (anyNA(value) || any(!nzchar(value))) {
      stop("SNIRF string array must contain non-empty values: ", dataset,
           call. = FALSE)
    }
  }
  value
}

# rhdf5 maps R's column-major dimensions in reverse HDF5 order.
.snirf_read_matrix <- function(path, dataset, tree) {
  if (.snirf_rank(tree, dataset) != 2L) {
    stop("SNIRF field must be a rank-2 matrix: ", dataset, call. = FALSE)
  }
  value <- rhdf5::h5read(path, dataset)
  if (!is.numeric(value) || length(dim(value)) != 2L ||
      any(dim(value) == 0L) || anyNA(value) || any(!is.finite(value))) {
    stop("SNIRF matrix must be non-empty and finite: ", dataset,
         call. = FALSE)
  }
  unname(t(as.matrix(value)))
}

.snirf_direct <- function(tree, group, type = NULL) {
  out <- tree[tree$group == group, , drop = FALSE]
  if (!is.null(type)) {
    out <- out[out$otype == type, , drop = FALSE]
  }
  out
}

.snirf_indexed_groups <- function(tree, parent, stem, allow_empty = FALSE) {
  rows <- .snirf_direct(tree, parent, "H5I_GROUP")
  pattern <- paste0("^", stem, "([0-9]+)$")
  keep <- grepl(pattern, rows$name)
  rows <- rows[keep, , drop = FALSE]
  if (!nrow(rows)) {
    if (allow_empty) return(character())
    stop("Required indexed SNIRF group is missing: ", parent, "/", stem,
         "1", call. = FALSE)
  }
  indices <- as.integer(sub(pattern, "\\1", rows$name))
  expected <- seq_len(max(indices))
  if (anyDuplicated(indices) || !identical(sort(indices), expected)) {
    paths <- paste(rows$.path, collapse = ", ")
    stop("Indexed SNIRF groups must be contiguous from 1: ", paths,
         call. = FALSE)
  }
  rows$.path[order(indices)]
}

.snirf_select_root <- function(tree, nirs_index) {
  rows <- .snirf_direct(tree, "/", "H5I_GROUP")
  names <- rows$name[grepl("^nirs([0-9]+)?$", rows$name)]
  if (!length(names)) {
    stop("SNIRF contains no /nirs root", call. = FALSE)
  }
  if ("nirs" %in% names && any(grepl("^nirs[0-9]+$", names))) {
    stop("Conflicting SNIRF roots: ",
         paste0("/", sort(names), collapse = ", "), call. = FALSE)
  }
  if (identical(names, "nirs")) {
    if (nirs_index != 1L) {
      stop("Requested nirs_index does not exist: ", nirs_index, call. = FALSE)
    }
    return("/nirs")
  }
  indexed <- as.integer(sub("^nirs", "", names))
  expected <- seq_len(max(indexed))
  if (anyDuplicated(indexed) || !identical(sort(indexed), expected)) {
    stop("Indexed SNIRF roots must be contiguous from /nirs1: ",
         paste0("/", sort(names), collapse = ", "), call. = FALSE)
  }
  requested <- paste0("/nirs", nirs_index)
  if (!requested %in% rows$.path) {
    stop("Requested SNIRF root does not exist: ", requested, call. = FALSE)
  }
  requested
}

.snirf_select_data <- function(tree, root, data_index) {
  groups <- .snirf_indexed_groups(tree, root, "data")
  requested <- paste0(root, "/data", data_index)
  if (!requested %in% groups) {
    stop("Requested SNIRF data block does not exist: ", requested,
         call. = FALSE)
  }
  requested
}

.snirf_time_factor <- function(unit) {
  factors <- c(s = 1, ms = 1e-3, us = 1e-6)
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) ||
      !unit %in% names(factors)) {
    stop("Unsupported SNIRF TimeUnit: ", paste(unit, collapse = ", "),
         call. = FALSE)
  }
  unname(factors[[unit]])
}

.snirf_length_factor <- function(unit) {
  factors <- c(m = 1, cm = 1e-2, mm = 1e-3)
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) ||
      !unit %in% names(factors)) {
    stop("Unsupported SNIRF LengthUnit: ", paste(unit, collapse = ", "),
         call. = FALSE)
  }
  unname(factors[[unit]])
}

.snirf_uniform_step <- function(time) {
  if (length(time) < 2L) return(NA_real_)
  delta <- diff(time)
  step <- stats::median(delta)
  if (!is.finite(step) || step <= 0) return(NA_real_)
  tolerance <- max(1e-12, abs(step) * 1e-9)
  if (max(abs(delta - step)) <= tolerance) step else NA_real_
}

.snirf_positive_whole <- function(value, path) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value > .Machine$integer.max ||
      value != floor(value)) {
    stop("SNIRF index must be a finite positive whole number: ", path,
         call. = FALSE)
  }
  as.integer(value)
}

.snirf_write_scalar <- function(path, dataset, value) {
  mode <- if (is.character(value)) "character" else if (is.integer(value)) {
    "integer"
  } else {
    "double"
  }
  args <- list(
    file = path, dataset = dataset, dims = NULL, storage.mode = mode
  )
  if (mode == "character") args$size <- NULL
  do.call(rhdf5::h5createDataset, args)
  rhdf5::h5write(value, path, dataset)
}

.snirf_write_vector <- function(path, dataset, value) {
  if (!length(value)) return(invisible(NULL))
  mode <- if (is.character(value)) "character" else if (is.integer(value)) {
    "integer"
  } else {
    "double"
  }
  args <- list(
    file = path, dataset = dataset, dims = length(value),
    storage.mode = mode
  )
  if (mode == "character") args$size <- NULL
  do.call(rhdf5::h5createDataset, args)
  rhdf5::h5write(value, path, dataset)
  invisible(NULL)
}

.snirf_write_matrix <- function(path, dataset, value) {
  value <- as.matrix(value)
  rhdf5::h5createDataset(path, dataset, dims = rev(dim(value)),
                         storage.mode = "double")
  rhdf5::h5write(t(value), path, dataset)
  invisible(NULL)
}
