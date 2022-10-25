#' Calibrate a single radiocarbon determination
#'
#' Uses the supplied calibration curve to take a single radiocardbon
#' determination and uncertainty and calculate the calendar age probability
#' density for it.
#'
#' @param c14_determination A single observed radiocarbon determination (c14 age)
#' @param c14_uncertainty The uncertainty of the radiocarbon determination
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#'
#' @return A vector containing the probability for each calendar age in the
#' calibration data
#' @export
#'
#' @examples
#' CalibrateSingleDetermination(51020, 35, intcal20)
CalibrateSingleDetermination <- function(
    c14_determination,
    c14_uncertainty,
    calibration_curve) {

  c14_ages = calibration_curve$c14_age
  c14_sigs = calibration_curve$c14_sig

  probabilities <- stats::dnorm(
    c14_determination, mean=c14_ages, sd=sqrt(c14_sigs^2 + c14_uncertainty^2))
  probabilities <- probabilities / sum(probabilities)
  return(probabilities)
}


#' Find the summed probability function (SPD) for a set of observations
#'
#' Takes a set of radiocarbon determinations and uncertainties and independently
#' calibrates each one, and then averages them to give the SPD estimate.
#'
#' @param calendar_age_range An array of length 2 with the start and end
#' calendar age to calculate the SPD over
#' @param c14_determinations An array of observed radiocarbon determinations
#' @param c14_uncertainties An array of the radiocarbon determinations
#' uncertainties. Must be the same length as `c14_determinations`.
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#'
#' @return A data frame with one column `calendar_age` containing the calendar
#' ages, and the other column `probability` containing the probability at that
#' calendar age
#' @export
#'
#' @examples
#' FindSPD(
#'    calendar_age_range=c(600, 1600),
#'    c14_determinations=c(602, 805, 1554),
#'    c14_uncertainties=c(35, 34, 45),
#'    calibration_curve=intcal20)
FindSPD <- function(
    calendar_age_range,
    c14_determinations,
    c14_uncertainties,
    calibration_curve) {

  # TODO(error-handling): check determinations and uncertainties are same length

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
    CalibrateSingleDetermination,
    c14_determinations,
    c14_uncertainties,
    MoreArgs = list(calibration_curve = interpolated_calibration_data))

  probabilities_per_calendar_age = apply(individual_probabilities, 1, sum) /
    dim(individual_probabilities)[2]

  SPD <- data.frame(
    calendar_age=new_calendar_ages,
    probability=probabilities_per_calendar_age)
  return(SPD)
}
