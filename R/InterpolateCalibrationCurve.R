#' Interpolate calendar ages for a calibration curve
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param new_calendar_ages A scalar or vector containing calendar ages (in yr BP) to
#' interpolate the calibration curve to. If not provided (and `NA` is given) uses the range
#' from the minimum calendar age to the maximum calendar age of the original calibration curve
#' spaced by 1.
#' @param F14C_outputs `TRUE` if only F14C concentrations are required, `FALSE`
#' if only radiocarbon age BP are required and `NA` if both are required for the new curve.
#'
#' @return A new dataframe with entries for the interpolated `c14_age`, and
#' `c14_sig`, `f14c` and `f14c_sig` values at the `calendar_age_BP` values given
#' in `new_calendar_ages`
#' @export
#'
#' @examples
#' InterpolateCalibrationCurve(51020.5, intcal20)
#' InterpolateCalibrationCurve(c(51020.5, 51021.5), intcal20, TRUE)
#' InterpolateCalibrationCurve(c(51020.5, 51021.5), intcal20, FALSE)
#' InterpolateCalibrationCurve(NA, intcal20)
InterpolateCalibrationCurve <- function(new_calendar_ages, calibration_curve, F14C_outputs = NA) {

  .CheckCalibrationCurve(NULL, calibration_curve, NA)
  if (!any(is.na(new_calendar_ages))) {
    checkmate::assertNumeric(new_calendar_ages)
  } else {
    start_age <- floor(min(calibration_curve$calendar_age_BP))
    end_age <- ceiling(max(calibration_curve$calendar_age_BP))
    diff <- 1
    new_calendar_ages <- seq(start_age, end_age, by=diff)
  }

  new_calibration_curve <- data.frame(calendar_age_BP=new_calendar_ages)

  calendar_ages <- calibration_curve$calendar_age_BP

  if (is.na(F14C_outputs) || F14C_outputs == FALSE) {
    calibration_curve = .AddC14ageColumns(calibration_curve)
    new_calibration_curve$c14_age <- stats::approx(
      calendar_ages, calibration_curve$c14_age, new_calendar_ages, rule=2)$y
    new_calibration_curve$c14_sig <- stats::approx(
      calendar_ages, calibration_curve$c14_sig, new_calendar_ages, rule=2)$y
  }
  if (is.na(F14C_outputs) || F14C_outputs == TRUE) {
    calibration_curve = .AddF14cColumns(calibration_curve)
    new_calibration_curve$f14c <- stats::approx(
      calendar_ages, calibration_curve$f14c, new_calendar_ages, rule=2)$y
    new_calibration_curve$f14c_sig <- stats::approx(
      calendar_ages, calibration_curve$f14c_sig, new_calendar_ages, rule=2)$y
  }

  return(new_calibration_curve)
}
