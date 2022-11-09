#' Find the summed probability distribution (SPD) for a set of observations
#'
#' Takes a set of radiocarbon determinations and uncertainties and independently
#' calibrates each one, and then averages them to give the SPD estimate.
#'
#' @inheritParams CalibrateSingleDetermination
#' @param calendar_age_range A vector of length 2 with the start and end
#' calendar age to calculate the SPD over
#' @param c14_determinations A vector of observed radiocarbon determinations
#' @param c14_sigmas A vector of the radiocarbon determinations
#' uncertainties (1-sigma). Must be the same length as `c14_determinations`.
#'
#' @return A data frame with one column `calendar_age` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age
#' @export
#'
#' @examples
#' FindSummedProbabilityDistribution(
#'    calendar_age_range=c(600, 1600),
#'    c14_determinations=c(602, 805, 1554),
#'    c14_sigmas=c(35, 34, 45),
#'    calibration_curve=intcal20)
FindSummedProbabilityDistribution <- function(
    calendar_age_range,
    c14_determinations,
    c14_sigmas,
    calibration_curve) {

  arg_check <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(
    calendar_age_range, len = 2, any.missing = FALSE, add = arg_check)
  .CheckInputData(
    arg_check, c14_determinations, c14_sigmas, calibration_curve)
  checkmate::reportAssertions(arg_check)

  calibration_data_range = range(calibration_curve$calendar_age)
  range_margin = 400

  new_calendar_ages <- seq(
    max(calibration_data_range[1], calendar_age_range[1] - range_margin),
    min(calibration_data_range[2], calendar_age_range[2] + range_margin),
    by = 1)

  interpolated_calibration_data <- InterpolateCalibrationCurve(
    new_calendar_ages, calibration_curve)

  # Each column of matrix gives probabilities per calendar age for a given determination
  individual_probabilities <- mapply(
    .ProbabilitiesForSingleDetermination,
    c14_determinations,
    c14_sigmas,
    MoreArgs = list(calibration_curve = interpolated_calibration_data))

  probabilities_per_calendar_age = apply(individual_probabilities, 1, sum) /
    dim(individual_probabilities)[2]

  SPD <- data.frame(
    calendar_age=new_calendar_ages,
    probability=probabilities_per_calendar_age)
  return(SPD)
}
