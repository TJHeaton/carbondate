#' Find the summed probability distribution (SPD) for a set of observations
#'
#' Takes a set of radiocarbon determinations and uncertainties and independently
#' calibrates each one, and then averages them to give the SPD estimate.
#'
#' @inheritParams CalibrateSingleDetermination
#' @param calendar_age_range_BP A vector of length 2 with the start and end
#' calendar age BP to calculate the SPD over
#' @param rc_determinations A vector of observed radiocarbon determinations
#' @param rc_sigmas A vector of the radiocarbon determinations
#' uncertainties (1-sigma). Must be the same length as `rc_determinations`.
#' @param calibration_curve A dataframe which must contain one column `calendar_age_BP`, and also
#' columns `c14_age` and `c14_sig` or `f14c` and `f14c_sig` (or both sets).
#' This format matches the curves supplied with this package e.g. [carbondate::intcal20], [carbondate::intcal13],
#' which contain all 5 columns.
#'
#' @return A data frame with one column `calendar_age_BP` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age
#' @export
#'
#' @examples
#' # An example using 14C age BP and the IntCal 20 curve
#' FindSummedProbabilityDistribution(
#'    calendar_age_range_BP=c(600, 1600),
#'    rc_determinations=c(602, 805, 1554),
#'    rc_sigmas=c(35, 34, 45),
#'    calibration_curve=intcal20)
#'
#' # An example using F14C concentrations and the IntCal 13 curve
#' FindSummedProbabilityDistribution(
#'    calendar_age_range_BP=c(600, 1600),
#'    rc_determinations=c(0.8, 0.85, 0.9),
#'    rc_sigmas=c(0.01, 0.015, 0.012),
#'    F14C_inputs=TRUE,
#'    calibration_curve=intcal13)
FindSummedProbabilityDistribution <- function(
    calendar_age_range_BP,
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE) {

  arg_check <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(calendar_age_range_BP, len = 2, any.missing = FALSE, add = arg_check)
  checkmate::assert_flag(F14C_inputs, add=arg_check)
  .CheckInputData(arg_check, rc_determinations, rc_sigmas, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  checkmate::reportAssertions(arg_check)

  calibration_data_range = range(calibration_curve$calendar_age)
  range_margin = 400

  new_calendar_ages <- seq(
    max(calibration_data_range[1], calendar_age_range_BP[1] - range_margin),
    min(calibration_data_range[2], calendar_age_range_BP[2] + range_margin),
    by = 1)

  interpolated_calibration_data <- InterpolateCalibrationCurve(
    new_calendar_ages, calibration_curve, F14C_inputs)

  # Each column of matrix gives probabilities per calendar age for a given determination
  individual_probabilities <- mapply(
    .ProbabilitiesForSingleDetermination,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(F14C_inputs = F14C_inputs, calibration_curve = interpolated_calibration_data))

  probabilities_per_calendar_age = apply(individual_probabilities, 1, sum) /
    dim(individual_probabilities)[2]

  SPD <- data.frame(
    calendar_age_BP=new_calendar_ages,
    probability=probabilities_per_calendar_age)
  return(SPD)
}
