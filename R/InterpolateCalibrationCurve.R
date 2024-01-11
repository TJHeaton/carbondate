#' Interpolate calendar ages for a calibration curve
#'
#' @inheritParams FindSummedProbabilityDistribution
#' @param new_calendar_ages_BP A scalar or vector containing calendar ages (in yr BP) to
#' interpolate the calibration curve to. If not provided (and `NA` is given) uses the range
#' from the minimum calendar age to the maximum calendar age of the original calibration curve
#' spaced by 1.
#' @param F14C_outputs `TRUE` if only F14C concentrations are required, `FALSE`
#' if only radiocarbon age BP are required and `NA` if both are required for the new curve.
#'
#' @return A new dataframe with entries for the interpolated `c14_age`, and
#' `c14_sig`, `f14c` and `f14c_sig` values at the `calendar_age_BP` values given
#' in `new_calendar_ages_BP`
#' @export
#'
#' @examples
#'
#' # Interpolate a single calendar age. Generates both radiocarbon scales.
#' InterpolateCalibrationCurve(51020.5, intcal20)
#'
#' # Interpolate at two calendar ages. Generates F14C ages only.
#' InterpolateCalibrationCurve(c(51020.5, 51021.5), intcal20, TRUE)
#'
#' # Interpolate at two calendar ages. Generate 14C age BP only.
#' InterpolateCalibrationCurve(c(51020.5, 51021.5), intcal20, FALSE)
#'
#' # Interpolate at every integer calendar age within the range of dates (for intcal20 this is 0
#' # to 55000 BP), and add for both radiocarbon scales.
#' cal_curve = InterpolateCalibrationCurve(NA, intcal20)
InterpolateCalibrationCurve <- function(new_calendar_ages_BP, calibration_curve, F14C_outputs = NA) {

  arg_check <- .InitializeErrorList()
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  .CheckFlag(arg_check, F14C_outputs)
  if (!any(is.na(new_calendar_ages_BP))) {
    .CheckNumberVector(arg_check, new_calendar_ages_BP)
  } else {
    start_age <- floor(min(calibration_curve$calendar_age_BP))
    end_age <- ceiling(max(calibration_curve$calendar_age_BP))
    diff <- 1
    new_calendar_ages_BP <- seq(start_age, end_age, by=diff)
  }
  .ReportErrors(arg_check)

  new_calibration_curve <- data.frame(calendar_age_BP=new_calendar_ages_BP)

  calendar_ages <- calibration_curve$calendar_age_BP

  if (is.na(F14C_outputs) || F14C_outputs == FALSE) {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    new_calibration_curve$c14_age <- stats::approx(
      calendar_ages, calibration_curve$c14_age, new_calendar_ages_BP, rule=2)$y
    new_calibration_curve$c14_sig <- stats::approx(
      calendar_ages, calibration_curve$c14_sig, new_calendar_ages_BP, rule=2)$y
  }
  if (is.na(F14C_outputs) || F14C_outputs == TRUE) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    new_calibration_curve$f14c <- stats::approx(
      calendar_ages, calibration_curve$f14c, new_calendar_ages_BP, rule=2)$y
    new_calibration_curve$f14c_sig <- stats::approx(
      calendar_ages, calibration_curve$f14c_sig, new_calendar_ages_BP, rule=2)$y
  }

  return(new_calibration_curve)
}
