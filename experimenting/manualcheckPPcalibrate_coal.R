library(boot)
set.seed(13)


true_theta <- as.vector(coal[,1])
true_theta <- 365 * (true_theta - min(true_theta))
n_observed <- length(true_theta)

rc_sigmas <- rep(30, n_observed)
rc_determinations <- rnorm(
  n = n_observed,
  mean = approx(
    x = intcal20$calendar_age_BP,
    y = intcal20$c14_age,
    xout = true_theta)$y,
  sd = rc_sigmas
)

calendar_grid_resolution <- 100
calendar_age_range <- c(0,50000)
rate_s <- NA
rate_h <- NA
prior_n_internal_changepoints_lambda <- 5
k_max_internal_changepoints <- 30
rescale_factor_rev_jump <- 0.9
n_iter <- 1000
F14C_inputs <- FALSE
use_F14C_space <- FALSE
calibration_curve <- intcal20
calendar_grid_resolution <- 100

Test_Output <- PPcalibrate(
  rc_determinations = rc_determinations,
  rc_sigmas = rc_sigmas,
  calibration_curve = calibration_curve,
  prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
  prior_h_rate = 0.1,
  F14C_inputs = F14C_inputs,
  use_F14C_space = use_F14C_space,
  k_max_internal_changepoints = k_max_internal_changepoints,
  rescale_factor_rev_jump = rescale_factor_rev_jump,
  calendar_grid_resolution = calendar_grid_resolution,
  n_iter = n_iter)

rm(true_theta, n_observed)
