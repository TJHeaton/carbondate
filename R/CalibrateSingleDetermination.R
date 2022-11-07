#' Calibrate a single radiocarbon determination
#'
#' Uses the supplied calibration curve to take a single radiocardbon
#' determination and uncertainty and calculate the calendar age probability
#' density for it.
#'
#' @param c14_determination A single observed radiocarbon determination (c14 age)
#' @param c14_sigma The uncertainty of the radiocarbon determination
#' @param calibration_curve A dataframe which should contain at least 3 columns
#' entitled calendar_age, c14_age and c14_sig.
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
    c14_sigma,
    calibration_curve) {

  arg_check <- checkmate::makeAssertCollection()
  checkmate::assertNumber(c14_determination, add = arg_check)
  checkmate::assertNumber(c14_sigma, add = arg_check)
  .CheckCalibrationCurve(arg_check, calibration_curve)
  checkmate::reportAssertions(arg_check)

  c14_ages = calibration_curve$c14_age
  c14_sigs = calibration_curve$c14_sig

  probabilities <- stats::dnorm(
    c14_determination, mean=c14_ages, sd=sqrt(c14_sigs^2 + c14_sigma^2))
  probabilities <- probabilities / sum(probabilities)
  return(probabilities)
}
