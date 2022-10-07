#' Calibrates a single radiocarbon determination
#'
#' @param rc_determination A single observed radiocarbon determination (c14 age)
#' @param rc_uncertainty The uncertainty of the radiocarbon determination
#' @param calibration_data A dataframe which should contain one column entitled
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
    rc_determination,
    rc_uncertainty,
    calibration_data) {

  c14_ages = calibration_data$c14_age
  c14_sigs = calibration_data$c14_sig

  probabilities <- stats::dnorm(
    rc_determination, mean=c14_ages, sd=sqrt(c14_sigs^2 + rc_uncertainty^2))
  probabilities <- probabilities/sum(probabilities)
  return(probabilities)
}
