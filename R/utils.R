#' Interpolate calendar ages for a calibration curve
#'
#' @param new_calendar_ages A scalar or vector containing calendar ages to
#' interpolate the calibration curve to
#' @param calibration_curve A dataframe which should contain one column entitled
#' c14_age and one column entitled c14_sig.
#' This format matches [carbondate::intcal20].
#'
#' @return A new dataframe with entries for the interpolated `c14_age` and
#' `c14_sig` values at the `calendar_age` values given in `new_calendar_ages`
#' @export
#'
#' @examples
#' InterpolateCalibrationCurve(51020.5, intcal20)
#' InterpolateCalibrationCurve(c(51020.5, 51021.5), intcal20)
InterpolateCalibrationCurve <- function(
    new_calendar_ages,
    calibration_curve) {

  calendar_ages =  calibration_curve$calendar_age
  c14_ages = calibration_curve$c14_age
  c14_sigs = calibration_curve$c14_sig

  new_c14_ages = stats::approx(
    calendar_ages, c14_ages, new_calendar_ages, rule=2)$y
  new_c14_sigs = stats::approx(
    calendar_ages, c14_sigs, new_calendar_ages, rule=2)$y

  new_calibration_curve = data.frame(
    calendar_age=new_calendar_ages,
    c14_age=new_c14_ages,
    c14_sig=new_c14_sigs)
  return(new_calibration_curve)
}


# The as.row function which we will use to count number of different clusters
as.row <- function(x) {
  y <- NULL
  x <- as.matrix(x)
  if (nrow(x) > ncol(x)) {
    y <- t(x)
  } else {
    y <- x
  }
  y
}
