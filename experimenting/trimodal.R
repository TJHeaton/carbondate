# 10th October 2023

# Study to see how well the densities are reconstructed
library(carbondate)
###############################################################################
# Create some observed data from the clusters according to probabilities

num_observations <- 200
c14_sig <- rep(50, num_observations)

# Choose base function of calendar age
weights_true <- c(0.4, 0.3, 0.3)
phi_true <- c(14400, 13800, 12900)
tau_true <- 1 / c(200, 150, 150)^2
n_clust_true <- length(weights_true)

true_density <- data.frame(x = intcal20$calendar_age, y = 0)
true_density$y <- 0
for(i in seq_along(weights_true)) {
  true_density$y <- true_density$y +
    weights_true[i] * dnorm(
      true_density$x, mean = phi_true[i], sd = 1/sqrt(tau_true[i]))
}

cluster_identifiers_true <- sample(
  1:n_clust_true, num_observations, replace = TRUE, prob = weights_true)
calendar_ages_true <- rnorm(
  num_observations,
  mean = phi_true[cluster_identifiers_true],
  sd = 1 / sqrt(tau_true[cluster_identifiers_true]))
hist(calendar_ages_true, breaks = 20, freq = FALSE)
lines(true_density, col="red")

# Create some radiocarbon determinations
interpolated_calibration_curve <- InterpolateCalibrationCurve(
  new_calendar_ages_BP = calendar_ages_true, calibration_curve = intcal20)
interpolated_c14_age <- interpolated_calibration_curve$c14_age
interpolated_c14_sig <- interpolated_calibration_curve$c14_sig

# Sample some calibration curve values
xcalcurve <- rnorm(num_observations, interpolated_c14_age, interpolated_c14_sig)
c14_ages <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)

###############################################################################
# Implement the Neal version of the DPMM
polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = c14_ages,
  rc_sigmas = c14_sig,
  F14C_inputs = FALSE,
  calibration_curve = intcal20,
  n_iter = 100000,
  n_thin = 50)

PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  n_posterior_samples = 5000,
  show_SPD = TRUE)


par(mfrow = c(1,1))
