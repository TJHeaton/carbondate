# Testing mixed calibration

# Read in PP_calibrate_mixed function
load_all()

# Test marine calibration
CalibrateSingleDetermination(8060,
                             35,
                             marine20,
                             delta_r = 400,
                             delta_r_sig = 20,
                             plot_output = TRUE)


# Create some variables
x <- pp_uniform_phase$c14_age
x_sig <- pp_uniform_phase$c14_sig
n_obs <- length(x)

sample_source <- rep("NH", n_obs)
delta_r <- rep(0, n_obs)
delta_r_sig <- rep(0, n_obs)

sample_source[1] <- "Marine"
delta_r[1:40] <- 100
delta_r_sig[1:40] <- 50


# Test calibration
test_pp_mixed <- PPcalibrateMixedCurves(rc_determinations = x,
                  rc_sigmas = x_sig,
                  delta_r = delta_r,
                  delta_r_sig = delta_r_sig,
                  sample_source = sample_source)

# Default plot with 2 sigma interval
PlotPosteriorMeanRate(test_pp_mixed)




# Manually test some individual functions
bounding_range_prob_cutoff <- 0.001
calendar_age_grid <- seq(443, 642, by = 1)
use_F14C_space <- FALSE


# Reduce x by correct amount
x_adj <- x - delta_r
x_sig_adj <- sqrt(x_sig^2 + delta_r_sig^2)

rc_determinations <- x_adj
rc_sigmas <- x_sig_adj

sample_source <- factor(sample_source,
                        levels = c("NH",
                                   "SH",
                                   "Marine"))

cal_curves <- list(NH = intcal20,
                   SH = shcal20,
                   Marine = marine20)

min(cal_curves[[1]]$calendar_age_BP)

max(sapply(cal_curves, function(x) min(x$calendar_age_BP)))

new_likelihood_calendar_ages_from_calibration_curve <- mapply(
  .CalendarAgeLikelihoodGivenMultipleCurves,
  rc_determinations,
  rc_sigmas,
  sample_source,
  MoreArgs = list(
    theta = calendar_age_grid,
    F14C_inputs = use_F14C_space,
    cal_curves = cal_curves))

old_likelihood_calendar_ages_from_calibration_curve <- mapply(
  .CalendarAgeLikelihoodGivenCurve,
  rc_determinations,
  rc_sigmas,
  MoreArgs = list(
    theta = calendar_age_grid,
    F14C_inputs = use_F14C_space,
    calibration_curve = intcal20))

identical(new_likelihood_calendar_ages_from_calibration_curve,
          old_likelihood_calendar_ages_from_calibration_curve)



plot(calendar_age_grid,
     new_likelihood_calendar_ages_from_calibration_curve[,30],
     type = "l" )



new_bounds <- .FindBoundingCalendarRangeMultipleCurves(
  rc_determinations,
  rc_sigmas,
  sample_source,
  F14C_inputs = use_F14C_space,
  cal_curves = cal_curves,
  prob_cutoff = bounding_range_prob_cutoff)

old_bounds <- .FindBoundingCalendarRange(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = intcal20,
  F14C_inputs = use_F14C_space,
  prob_cutoff = bounding_range_prob_cutoff)

identical(new_bounds,
          old_bounds)
