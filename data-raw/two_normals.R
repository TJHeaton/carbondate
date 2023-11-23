library(usethis)

set.seed(2)


# Create some observed data from the clusters according to probabilities
num_observations <- 50

# Choose base function of calendar age
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2
n_clust_true <- length(weights_true)

true_density <- data.frame(x = intcal20$calendar_age_BP, y = 0)
true_density$y <- 0
for(i in 1:length(weights_true)) {
  true_density$y <- true_density$y +
    weights_true[i] * dnorm(
      true_density$x, mean = cluster_means_true_calBP[i], sd = 1/sqrt(cluster_precisions_true[i]))
}

cluster_identifiers_true <- sample(
  1:n_clust_true, num_observations, replace = TRUE, prob = weights_true)
calendar_ages_true <- rnorm(
  num_observations,
  mean = cluster_means_true_calBP[cluster_identifiers_true],
  sd = 1 / sqrt(cluster_precisions_true[cluster_identifiers_true]))
hist(calendar_ages_true, breaks = 20, freq = FALSE)
lines(true_density, col="red")

# Create some radiocarbon determinations
calibration_curve_at_true_calendar_ages <- InterpolateCalibrationCurve(
  new_calendar_ages = calendar_ages_true, calibration_curve = intcal20)
c14_age_for_true_calendar_ages <- calibration_curve_at_true_calendar_ages$c14_age
c14_sig_for_true_calendar_ages <- calibration_curve_at_true_calendar_ages$c14_sig

# Sample from the true C14 age and C14 sig values
xcalcurve <- rnorm(num_observations, c14_age_for_true_calendar_ages, c14_sig_for_true_calendar_ages)
c14_sig <- rep(25, num_observations)
c14_ages <- rnorm(num_observations, mean = xcalcurve, sd = c14_sig)

two_normals <- data.frame(
  c14_age = rnorm(num_observations, mean = xcalcurve, sd = c14_sig), c14_sig = c14_sig)
two_normals$f14c <- exp(-two_normals$c14_age / 8033)
two_normals$f14c_sig <- two_normals$f14c * two_normals$c14_sig / 8033

usethis::use_data(two_normals, overwrite = TRUE)
