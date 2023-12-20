.FindIntegral <- function(rate_s, rate_h) {
  nchangepoints <- length(rate_s)
  nheights <- length(rate_h)

  if(nheights != (nchangepoints - 1)) stop("Incompatible information on rate s and h")

  sum(rate_h * diff(rate_s))
}



# Find likelihood (given calibration curve) of:
# vector of calendar ages theta for a single 14C observation
.CalendarAgeLikelihoodGivenCurve <- function(
    rc_determination,
    rc_sigma,
    theta,
    F14C_inputs,
    calibration_curve)
{
  if (F14C_inputs) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$f14c
    calcurve_rc_sigs <- calibration_curve$f14c_sig
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$c14_age
    calcurve_rc_sigs <- calibration_curve$c14_sig
  }

  cal_curve_mu <- stats::approx(x = calibration_curve$calendar_age,
                         y = calcurve_rc_ages,
                         xout = theta)$y
  cal_curve_sigma <- stats::approx(x = calibration_curve$calendar_age,
                            y = calcurve_rc_sigs,
                            xout = theta)$y

  likelihood <- stats::dnorm(rc_determination, cal_curve_mu, sqrt(cal_curve_sigma^2 + rc_sigma^2))

  return(likelihood)
}

.FindCalendarRangeForSingleDetermination <- function(
    rc_determination, rc_sigma, F14C_inputs, calibration_curve, prob_cutoff) {

  if (F14C_inputs) {
    calibration_curve <- .AddF14cColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$f14c
    calcurve_rc_sigs <- calibration_curve$f14c_sig
  } else {
    calibration_curve <- .AddC14ageColumns(calibration_curve)
    calcurve_rc_ages <- calibration_curve$c14_age
    calcurve_rc_sigs <- calibration_curve$c14_sig
  }

  probabilities <- stats::dnorm(
    rc_determination, mean=calcurve_rc_ages, sd=sqrt(calcurve_rc_sigs^2 + rc_sigma^2))
  probabilities <- probabilities / sum(probabilities)
  cumulativeprobabilities <- cumsum(probabilities)

  min_potential_calendar_age <- calibration_curve$calendar_age_BP[
    min(which(cumulativeprobabilities > prob_cutoff))]
  max_potential_calendar_age <- calibration_curve$calendar_age_BP[
    min(which(cumulativeprobabilities > (1 - prob_cutoff)))]

  calendar_range <- c(min_potential_calendar_age, max_potential_calendar_age)

  return(calendar_range)
}

.FindBoundingCalendarRange <- function(
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs,
    prob_cutoff = 0.005)
{
  ##############################################################################
  ## Interpolate cal curve onto single year (regular) grid
  ## Must be regular calendar grid for individual_possible_calendar_ranges
  ## as this works with normalised vector of probabilities
  integer_cal_year_curve <- InterpolateCalibrationCurve(
    NA, calibration_curve, F14C_inputs)

  # Need to ensure rc_determinations and rc_sigmas are in form of F14C_inputs
  # This is not checked internally in this function
  individual_possible_calendar_ranges <- mapply(
    .FindCalendarRangeForSingleDetermination,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(
      F14C_inputs = F14C_inputs,
      calibration_curve = integer_cal_year_curve,
      prob_cutoff = prob_cutoff))

  min_potential_calendar_age <- min(
    individual_possible_calendar_ranges[1,])
  max_potential_calendar_age <- max(
    individual_possible_calendar_ranges[2,])

  cal_age_range <- c(min_potential_calendar_age,
                      max_potential_calendar_age)

  return(cal_age_range)
}


.FindTrimmedVectorAndIndices <- function(
    vector,
    prob_cutoff = 0.005) {

  cumulative_sum <- cumsum(vector)
  total_sum <- cumulative_sum[length(cumulative_sum)]

  min_index <- min(which(cumulative_sum > prob_cutoff * total_sum))
  max_index <- min(which(cumulative_sum > total_sum * (1 - prob_cutoff)))

  return (list(values = vector[min_index:max_index], start_index = min_index, end_index = max_index))
}


# Need care with using the sample command as sometimes we pass a single integer j.
# If use sample() then will draw from 1:j which is not what we want
# This resample function will stop this happening
.resample <- function(x, size, ...)
{
  if(length(x) <= 1) {
    if(!missing(size) && size == 0) x[FALSE]
    else x
  }
  else sample(x, size, ...)
}

