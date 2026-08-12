#' PhysioNIRS: governed NIRS data and SNIRF input/output
#'
#' `PhysioNIRS` maps continuous-wave SNIRF data into the
#' [PhysioCore::PhysioExperiment()] contract without changing measurement
#' order, probe geometry, stimulus tables, or the exact time base. Governed
#' optical-density, modified Beer-Lambert, motion-correction,
#' short-separation nuisance, quality/pruning, and live-neurofeedback
#' operations retain their log, distance, pathlength, geometry, filter,
#' identity, and reference conventions.
#'
#' SNIRF metadata may contain identifying subject and acquisition fields.
#' Callers are responsible for de-identification before sharing a file.
#'
#' @references
#' Shared Near Infrared Spectroscopy Format specification:
#' \url{https://fnirs.github.io/snirf/}
#'
#' @importFrom PhysioCore PhysioEvents PhysioExperiment appendProvenance
#' @importFrom PhysioCore defaultAssay getEvents samplingRate setEvents
"_PACKAGE"
