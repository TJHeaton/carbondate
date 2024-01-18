#' Calibrate a single radiocarbon determination
#'
#' Uses the supplied calibration curve to take a single radiocarbon
#' determination and uncertainty (expressed as the radiocarbon age BP or F14C concentration) and
#' calculate the calendar age probability density for it.
#'
#' @param rc_determination A single observed radiocarbon determination
#' (\eqn{{}^{14}}C BP age or F14C concentration)
#' @param rc_sigma The uncertainty of the radiocarbon determination in the same units
#' @param calibration_curve A dataframe which must contain one column `calendar_age_BP`, and also
#' columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both sets).
#' This format matches the curves supplied with this package e.g. [carbondate::intcal20],
#' which contain all 5 columns.
#' @param F14C_inputs `TRUE` if the provided rc_determinations are F14C concentrations and `FALSE`
#' if they are radiocarbon age BP. Defaults to `FALSE`.
#'
#' @export
#'
#' @return A data frame with one column `calendar_age_BP` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age.
#'
#' @examples
#' # Calibration of a single determination expressed as 14C age BP
#' calib = CalibrateSingleDetermination(31020, 35, intcal20)
#' plot(calib, type = "l", xlim = c(36000, 34000))
#'
#' # Calibration of a single determination expressed as F14C concentration
#' calib <- CalibrateSingleDetermination(
#'     0.02103493, 9.164975e-05, intcal20, F14C_inputs = TRUE)
#' plot(calib, type = "l", xlim = c(36000, 34000))
CalibrateSingleDetermination <- function(
    rc_determination, rc_sigma, calibration_curve, F14C_inputs = FALSE) {

  arg_check <- .InitializeErrorList()
  .CheckNumber(arg_check, rc_determination)
  .CheckNumber(arg_check, rc_sigma)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, F14C_inputs)
  .ReportErrors(arg_check)

  probabilities <- .ProbabilitiesForSingleDetermination(rc_determination, rc_sigma, F14C_inputs, calibration_curve)

  return(
    data.frame(
      calendar_age_BP=calibration_curve$calendar_age, probability=probabilities))
}


.ProbabilitiesForSingleDetermination <- function(
    rc_determination, rc_sigma, F14C_inputs, calibration_curve) {

  if (F14C_inputs) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$f14c
    calcurve_rc_sigs <- calibration_curve$f14c_sig
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$c14_age
    calcurve_rc_sigs <- calibration_curve$c14_sig
  }

  probabilities <- stats::dnorm(rc_determination, mean=calcurve_rc_ages, sd=sqrt(calcurve_rc_sigs^2 + rc_sigma^2))
  probabilities <- probabilities / sum(probabilities)
  return(probabilities)
}

