# Find likelihood (given calibration curve) of:
# vector of calendar ages theta for a single 14C observation
# when you specify the source of the sample
.CalendarAgeLikelihoodGivenMultipleCurves <- function(
    rc_determination,
    rc_sigma,
    sample_source,
    theta,
    F14C_inputs,
    cal_curves)
{
  calibration_curve <- cal_curves[[sample_source]]

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


# Find bounding calendar ages of a single 14C observation
# when you specify the source of the sample
.FindBoundingCalendarRangebyCurve <- function(
    rc_determinations,
    rc_sigmas,
    sample_source,
    F14C_inputs,
    prob_cutoff = 0.005,
    cal_curves)
{
  calibration_curve <- cal_curves[[sample_source]]

  ##############################################################################
  ## Interpolate cal curve onto single year (regular) grid
  ## Must be regular calendar grid for individual_possible_calendar_ranges
  ## as this works with normalised vector of probabilities
  integer_cal_year_curve <- InterpolateCalibrationCurve(
    NA, calibration_curve, F14C_inputs)

  # Need to ensure rc_determinations and rc_sigmas are in form of F14C_inputs
  # This is not checked internally in this function
  .FindCalendarRangeForSingleDetermination(
    rc_determinations,
    rc_sigmas,
    F14C_inputs = F14C_inputs,
    calibration_curve = integer_cal_year_curve,
    prob_cutoff = prob_cutoff)
}



# Find bounding calendar ages of a vector of 14C observation
# when you specify the source of the sample
.FindBoundingCalendarRangeMultipleCurves <- function(
    rc_determinations,
    rc_sigmas,
    sample_source,
    F14C_inputs,
    prob_cutoff = 0.005,
    cal_curves)
{
  individual_possible_calendar_ranges <- mapply(
    .FindBoundingCalendarRangebyCurve,
    rc_determinations,
    rc_sigmas,
    sample_source,
    MoreArgs = list(
      F14C_inputs = F14C_inputs,
      cal_curves = cal_curves,
      prob_cutoff = prob_cutoff))

  min_potential_calendar_age <- min(
    individual_possible_calendar_ranges[1,])
  max_potential_calendar_age <- max(
    individual_possible_calendar_ranges[2,])

  cal_age_range <- c(min_potential_calendar_age,
                     max_potential_calendar_age)

  return(cal_age_range)
}


# Add offset function
.AddOffset <- function(
    rc_determination,
    rc_sigma,
    delta_r,
    delta_r_sig,
    F14C_inputs)
{
  if(!F14C_inputs) {
    rc_determination_adj <- rc_determination - delta_r
    rc_sigma_adj <- sqrt(rc_sigma^2 + delta_r_sig^2)
  } else{
    rc_age_info <- .ConvertF14cTo14Cage(rc_determination, rc_sigmas)
    rc_age_info$c14_age <- rc_age_info$c14age - delta_r
    rc_age_info$c14_sig <- sqrt(rc_age_info$c14_sig^2 - delta_r_sig^2)
    f14c_info_adj <- .Convert14CageToF14c(rc_age_info$c14_age, rc_age_info$c14_sig)
    rc_determination_adj <- f14c_info_adj$f14c
    rc_sigma_adj <- f14c_info_adj$f14c_sig
  }

  adjusted_values <- list(rc_determination = rc_determination_adj,
                          rc_sigma = rc_sigma_adj)

  return(adjusted_values)
}

# Create list containing all curves for specific year
.CreateCalCurveList <- function(
    curve_year)
{
  if(curve_year != 2009) {
  cal_curves <- list(NH = get(paste("intcal", substr(curve_year, 3, 4), sep = "")),
                     SH = get(paste("shcal", substr(curve_year, 3, 4), sep = "")),
                     Marine = get(paste("marine", substr(curve_year, 3, 4), sep = "")))
  } else { # No new SHCal curve in 2009 so use SHCal04
    cal_curves <- list(NH = get(paste("intcal", substr(curve_year, 3, 4), sep = "")),
                       SH = shcal04,
                       Marine = get(paste("marine", substr(curve_year, 3, 4), sep = "")))
  }
  return(cal_curves)
}



