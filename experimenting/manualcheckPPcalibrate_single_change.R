# Create artificial dataset with a single change points in it
set.seed(15)
n_observed <- 200
n_iter <- 1000
n_thin <- 10
calibration_curve <- intcal20
calendar_age_range <- c(0,4000)
test_rc_sigma <- 15

# Create some true thetas
true_changepoint <- 2000
observed_age_range <- c(calendar_age_range[1],
                        true_changepoint)
true_theta <- seq(from = observed_age_range[1],
                  to = observed_age_range[2],
                  length = n_observed)

# Create rc_determinations
rc_sigmas <- rep(15, n_observed)
rc_determinations <- stats::rnorm(n = n_observed,
                                  mean = approx(
                                    x = calibration_curve$calendar_age_BP,
                                    y = calibration_curve$c14_age,
                                    xout = true_theta)$y,
                                  sd = rc_sigmas)


# Find the plausible calendar age range for plotting
min_cal_age <- max(0, calibration_curve$calendar_age_BP[which.min(
  abs(calibration_curve$c14_age - calendar_age_range[1]))
  ] - 250)
max_cal_age <- min(55000, calibration_curve$calendar_age_BP[which.min(
  abs(calibration_curve$c14_age - calendar_age_range[2]))
] + 250)
plot_calendar_age_range <- c(min_cal_age, max_cal_age)

# Find required
plot_radiocarbon_age_range <- approx(x = calibration_curve$calendar_age_BP,
                                     y = calibration_curve$c14_age,
                                     xout = calendar_age_range)$y + c(-100,100)

plot(intcal20$calendar_age_BP,
     intcal20$c14_age,
     type = "l",
     xlim = rev(plot_calendar_age_range),
     ylim = plot_radiocarbon_age_range)
rug(rc_determinations, side = 2)

Test_Output <- PP_fit_output <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = intcal20,
  calendar_age_range = calendar_age_range,
  calendar_grid_resolution = 10,
  n_iter = 100000
)

PlotPosteriorMeanRate(Test_Output)
PlotPosteriorChangePoints(Test_Output)
PlotPosteriorHeights(Test_Output)
PlotNumberOfInternalChanges(Test_Output)

